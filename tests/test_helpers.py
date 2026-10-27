"""Small synthetic checks; no research data, Scanpy, LIANA, or SCOTIA needed."""

import ast
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
from scipy import sparse

from pdac_macrophages.preprocessing import replace_nan_with_zero, sample_expression_zscores
from pdac_macrophages.spatial import exponential_model
from pdac_macrophages.statistics import outliers_iqr, erase_outliers
from pdac_macrophages.ligand_receptor import (
    combine_top_hits, subset_lr_data, load_scotia_results, summarize_lr,
    add_between_group_pvals,
)


ROOT = Path(__file__).resolve().parents[1]
CLEANED = ROOT / "MacrophagesAnalysisValidation.cleaned.ipynb"


class ObservationData:
    """Minimal obs-only AnnData stand-in for spatial and subsetting tests."""

    def __init__(self, obs):
        self.obs = obs

    def __len__(self):
        return len(self.obs)

    def __getitem__(self, mask):
        return ObservationData(self.obs.loc[mask])

    def copy(self):
        return ObservationData(self.obs.copy())


class NotebookTests(unittest.TestCase):
    def test_cleaned_notebook_is_unexecuted_and_syntax_valid(self):
        notebook = json.loads(CLEANED.read_text())
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                self.assertEqual(cell["outputs"], [])
                self.assertIsNone(cell["execution_count"])
                ast.parse("".join(cell["source"]))

    def test_opt_in_guards_stop_before_data_loading_and_lr_loop(self):
        notebook = json.loads(CLEANED.read_text())
        for index, flag in [(5, "RUN_ANALYSIS"), (119, "RUN_LR_ANALYSIS")]:
            with self.assertRaisesRegex(RuntimeError, "disabled"):
                exec("".join(notebook["cells"][index]["source"]), {flag: False})


class PreprocessingTests(unittest.TestCase):
    def test_sparse_nan_replacement_matches_dense_rule_without_densifying(self):
        matrix = sparse.csr_matrix([[1.0, np.nan], [0.0, np.inf]])
        expected = np.array([[1.0, 0.0], [0.0, np.inf]])
        # A sparse input whose toarray raises catches accidental densification.
        matrix.toarray = lambda: (_ for _ in ()).throw(AssertionError("densified"))
        result = replace_nan_with_zero(matrix)
        np.testing.assert_array_equal(result.toarray(), expected)
        self.assertTrue(np.isnan(matrix.data).any())

    def test_dense_nan_replacement_copies_input_and_keeps_infinities(self):
        matrix = np.array([[np.nan, -np.inf], [0.0, 4.0]])
        result = replace_nan_with_zero(matrix)
        self.assertTrue(sparse.isspmatrix_csr(result))
        np.testing.assert_array_equal(result.toarray(), [[0.0, -np.inf], [0.0, 4.0]])
        self.assertTrue(np.isnan(matrix[0, 0]))

    def test_sample_zscores_use_sample_means_and_ddof_one(self):
        matrix = np.array([[1., 2.], [3., 4.], [6., 9.]])
        expected = pd.DataFrame(
            [[-1 / np.sqrt(2), 1 / np.sqrt(2)]] * 2,
            index=["G1", "G2"], columns=["A", "B"],
        )
        for values in [matrix, sparse.csr_matrix(matrix)]:
            actual = sample_expression_zscores(values, ["A", "A", "B"], ["G1", "G2"])
            assert_frame_equal(actual, expected)

    def test_constant_gene_retains_undefined_zscores(self):
        result = sample_expression_zscores([[2.], [2.]], ["A", "B"], ["G"])
        self.assertTrue(result.isna().all().all())


class SpatialTests(unittest.TestCase):
    def setUp(self):
        self.data = ObservationData(pd.DataFrame({
            "sample_id": ["U1", "T1", "T1", "T1", "T1", "U1"],
            "fov": [1, 1, 1, 2, 1, 3],
            "annotation_majortypes": ["B", "A", "B", "A", "B", "A"],
            "CenterX_local_px": [0., 0., 3., 0., 0., 0.],
            "CenterY_local_px": [0., 0., 4., 0., 0., 0.],
        }))

    def test_distance_units_sums_and_sample_fov_isolation(self):
        first, second = exponential_model(self.data, "A", "B", radius=50)
        weight = np.exp(-5 * 0.18 / 50)
        np.testing.assert_allclose(first, [1 + weight, 0, 0])
        np.testing.assert_allclose(second, [weight, 1, 0])

    def test_identical_populations_include_self_distances(self):
        data = ObservationData(pd.DataFrame({
            "sample_id": ["T1", "T1"],
            "fov": [1, 1],
            "annotation_majortypes": ["A", "A"],
            "CenterX_local_px": [0., 3.],
            "CenterY_local_px": [0., 4.],
        }))
        first, second = exponential_model(data, "A", "A", radius=50)
        expected = [1 + np.exp(-5 * 0.18 / 50)] * 2
        np.testing.assert_allclose(first, expected)
        np.testing.assert_allclose(second, expected)

    def test_absent_populations_return_empty_arrays(self):
        first, second = exponential_model(self.data, "missing", "also missing")
        self.assertEqual((first.size, second.size), (0, 0))

    def test_nonpositive_or_nonfinite_radius_is_rejected(self):
        for radius in [0, -1, np.nan, np.inf]:
            with self.assertRaises(ValueError):
                exponential_model(self.data, "A", "B", radius=radius)


class OutlierTests(unittest.TestCase):
    def test_iqr_identifies_extreme_value(self):
        values = pd.Series([1., 2., 2., 3., 100.])
        assert_series_equal(outliers_iqr(values), pd.Series([100.], index=[4]))

    def test_filtering_is_group_specific_and_preserves_input(self):
        frame = pd.DataFrame({
            "interacting_cell": ["CAF"] * 10,
            "Condition": ["Treated"] * 5 + ["Untreated"] * 5,
            "mean_proximity": [1., 2., 2., 3., 100.] + [90., 100., 100., 110., 120.],
        })
        before = frame.copy(deep=True)
        result = erase_outliers(frame, "mean")
        self.assertNotIn(4, result.index)
        self.assertIn(6, result.index)
        assert_frame_equal(result, before.drop(index=4))
        assert_frame_equal(frame, before)

    def test_duplicate_indices_are_rejected(self):
        frame = pd.DataFrame(index=[0, 0])
        with self.assertRaisesRegex(ValueError, "unique"):
            erase_outliers(frame, "mean")


class LigandReceptorTests(unittest.TestCase):
    def test_hit_thresholds_union_and_deduplication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "L-R connectome full"
            second = root / "L-R cellphonedb full"
            first.mkdir()
            second.mkdir()
            filename = "liana_A_to_B_full.csv"
            pd.DataFrame({
                "ligand": ["L1", "L2", "L3"], "receptor": ["R1", "R2", "R3"],
                "scaled_weight": [1., .99, 2.],
            }).to_csv(first / filename, index=False)
            pd.DataFrame({
                "ligand": ["L1", "L2", "L4"], "receptor": ["R1", "R2", "R4"],
                "cellphone_pvals": [.01, .049, .05],
            }).to_csv(second / filename, index=False)
            actual = combine_top_hits("A", "B", root)
            expected = pd.DataFrame({
                "l_gene": ["L1", "L3", "L2"],
                "r_gene": ["R1", "R3", "R2"],
            })
            assert_frame_equal(actual, expected)

    def test_subsetting_uses_explicit_input_and_returns_copy(self):
        data = ObservationData(pd.DataFrame({"mac_type": ["A", "B", "C"]}))
        result = subset_lr_data(data, "A", "B")
        self.assertEqual(result.obs["mac_type"].tolist(), ["A", "B"])
        result.obs.iloc[0, 0] = "changed"
        self.assertEqual(data.obs.iloc[0, 0], "A")

    def test_scotia_tsv_loader_handles_upstream_typo_and_empty_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "empty.csv").touch()
            (root / "T1_fov1.csv").write_text("ligand_recptor\tlikelihood\nL_R\t0.5\n")
            result = load_scotia_results(root)
            self.assertEqual(result["ligand"].tolist(), ["L"])
            self.assertEqual(result["receptor"].tolist(), ["R"])
            self.assertEqual(result["fov"].tolist(), ["T1_fov1"])
            self.assertEqual(result["likelihood"].tolist(), [0.5])

    def test_scotia_loader_reports_missing_folders_and_bad_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FileNotFoundError):
                load_scotia_results(root / "missing")
            (root / "bad.csv").write_text("wrong\tlikelihood\nx\t0.5\n")
            with self.assertRaisesRegex(ValueError, "bad.csv.*missing columns"):
                load_scotia_results(root)

    def test_scotia_loader_reports_malformed_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "bad.csv").write_text("ligand_receptor\tlikelihood\ninvalid\t0.5\n")
            with self.assertRaisesRegex(ValueError, "expected ligand_receptor"):
                load_scotia_results(tmp)

    @staticmethod
    def interaction_frame(offset):
        return pd.DataFrame({
            "ligand": ["L1"] * 4 + ["L2"] * 3,
            "receptor": ["R1"] * 4 + ["R2"] * 3,
            "fov": ["f1", "f1", "f2", "f3", "f1", "f2", "f3"],
            "likelihood": np.array([1., 3., 5., 7., 0., 0., 0.]) + offset,
        })

    def test_summary_counts_unique_fovs_and_averages_rows(self):
        frame = self.interaction_frame(0)
        actual = summarize_lr(frame)
        expected = pd.DataFrame({
            "ligand": ["L1", "L2"], "receptor": ["R1", "R2"],
            "num_fovs_detected": [3, 3], "mean_likelihood": [4., 0.],
            "max_likelihood": [7., 0.],
        })
        assert_frame_equal(actual, expected)

    def test_comparison_pvals_qvals_fold_changes_and_csvs(self):
        treated, untreated = self.interaction_frame(2), self.interaction_frame(0)
        t_summary, u_summary = summarize_lr(treated), summarize_lr(untreated)
        t_before, u_before = treated.copy(deep=True), untreated.copy(deep=True)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = (treated, untreated, t_summary, u_summary)
            actual = add_between_group_pvals(*args, root / "t.csv", root / "u.csv")
            # Fixed reference values from the original helper on these synthetic
            # inputs, recorded before removing the private-notebook dependency.
            expected = pd.DataFrame({
                "ligand": ["L1", "L2"], "receptor": ["R1", "R2"],
                "n_fov_treated": [3, 3], "n_fov_untreated": [3, 3],
                "log2FC_T_over_U": [0.514573172737014, 30.897352854707609],
                "pval": [0.506555169049040, 0.046854177603874],
                "qval": [0.506555169049040, 0.093708355207748],
            })
            assert_frame_equal(actual, expected, check_exact=False, rtol=1e-12, atol=1e-14)
            columns = ["ligand", "receptor", "log2FC_T_over_U", "pval", "qval"]
            for name, summary in [("t", t_summary), ("u", u_summary)]:
                expected_csv = summary.merge(
                    expected[columns], on=["ligand", "receptor"], how="left"
                )
                assert_frame_equal(
                    pd.read_csv(root / f"{name}.csv"), expected_csv,
                    check_exact=False, rtol=1e-12, atol=1e-14,
                )
            self.assertTrue((actual["n_fov_treated"] == 3).all())
            self.assertTrue(np.isfinite(actual["log2FC_T_over_U"]).all())
        assert_frame_equal(treated, t_before)
        assert_frame_equal(untreated, u_before)

    def test_insufficient_fovs_preserve_empty_result_and_no_write(self):
        frame = self.interaction_frame(0).query("fov == 'f1'")
        summary = summarize_lr(frame)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "summary.csv"
            result = add_between_group_pvals(frame, frame, summary, summary, output, output)
            self.assertTrue(result.empty)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
