"""LIANA/SCOTIA result handling (original cells 108–113).

These helpers never run LIANA or SCOTIA. Statistical choices are preserved.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


DEFAULT_OUTPUT_ROOT = Path("LigandReceptor_Analysis_Output")


def combine_top_hits(ligand, receptor, out_root=DEFAULT_OUTPUT_ROOT):
    """Union Connectome scaled_weight >= 1 and CellPhoneDB p < 0.05 hits.

    The Connectome threshold is a score cutoff, not a significance test.
    Preserve input ordering and remove duplicate ligand/receptor gene pairs.
    """
    root = Path(out_root)
    filename = f"liana_{ligand}_to_{receptor}_full.csv"
    connectome = pd.read_csv(root / "L-R connectome full" / filename)
    cpdb = pd.read_csv(root / "L-R cellphonedb full" / filename)
    first = connectome.loc[connectome["scaled_weight"] >= 1.0, ["ligand", "receptor"]]
    second = cpdb.loc[cpdb["cellphone_pvals"] < 0.05, ["ligand", "receptor"]]
    combined = pd.concat([first, second], ignore_index=True)
    return combined.drop_duplicates().reset_index(drop=True).rename(
        columns={"ligand": "l_gene", "receptor": "r_gene"}
    )


def subset_lr_data(adata, ligand, receptor):
    """Copy the two requested cell populations without using notebook globals."""
    return adata[adata.obs["mac_type"].isin([ligand, receptor])].copy()


def load_scotia_results(folder_path):
    """Load tab-separated SCOTIA ``*.csv`` files, retaining legacy parsing.

    Empty files are skipped. Missing folders and malformed nonempty files
    raise useful errors rather than silently dropping results. FOV identity
    remains the filename stem; gene pairs still split at the first underscore.
    """
    folder = Path(folder_path)
    if not folder.is_dir():
        raise FileNotFoundError(f"SCOTIA result directory not found: {folder}")
    frames = []
    for file in sorted(folder.glob("*.csv")):
        try:
            frame = pd.read_csv(file, sep="\t", header=0)
        except pd.errors.EmptyDataError:
            continue
        if frame.empty:
            continue
        frame = frame.rename(columns={"ligand_recptor": "ligand_receptor"})
        required = {"ligand_receptor", "likelihood"}
        if missing := required.difference(frame.columns):
            raise ValueError(f"{file}: missing columns {sorted(missing)}")
        pairs = frame["ligand_receptor"].str.split("_", n=1, expand=True)
        if pairs.shape[1] != 2 or pairs.isna().any().any() or (pairs == "").any().any():
            raise ValueError(f"{file}: expected ligand_receptor values like LIGAND_RECEPTOR")
        frame["fov"] = file.stem
        frame[["ligand", "receptor"]] = pairs
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize_lr(df):
    """Summarize distinct FOV count and row-level mean/max likelihood per pair."""
    return (
        df.groupby(["ligand", "receptor"])
        .agg(
            num_fovs_detected=("fov", "nunique"),
            mean_likelihood=("likelihood", "mean"),
            max_likelihood=("likelihood", "max"),
        )
        .reset_index()
        .sort_values(["num_fovs_detected", "mean_likelihood"], ascending=False)
    )


def add_between_group_pvals(
    treated_df, untreated_df, treated_summary, untreated_summary,
    treated_csv_path, untreated_csv_path,
):
    """Preserve the original per-FOV Mann–Whitney tests and BH adjustment.

    Average likelihood per pair/FOV; require two detected FOVs per group.
    Use a 1e-9 pseudocount for log2 fold change. Write merged summaries to the
    supplied CSV paths and return comparisons. Return an empty frame without
    writing if no pairs qualify, matching the notebook's behavior.

    FOV independence, missing detections, and pooled directionality require
    scientific review; see docs/REVIEW.md. This does not validate those choices.
    """
    keys = ["ligand", "receptor", "fov"]
    treated = treated_df.groupby(keys, as_index=False)["likelihood"].mean()
    untreated = untreated_df.groupby(keys, as_index=False)["likelihood"].mean()
    pairs = pd.concat([
        treated[["ligand", "receptor"]], untreated[["ligand", "receptor"]]
    ]).drop_duplicates()
    records = []
    for _, row in pairs.iterrows():
        ligand, receptor = row["ligand"], row["receptor"]
        t_values = treated.loc[
            (treated["ligand"] == ligand) & (treated["receptor"] == receptor),
            "likelihood",
        ].values
        u_values = untreated.loc[
            (untreated["ligand"] == ligand) & (untreated["receptor"] == receptor),
            "likelihood",
        ].values
        if len(t_values) < 2 or len(u_values) < 2:
            continue
        _, p_value = mannwhitneyu(t_values, u_values, alternative="two-sided")
        fold_change = (t_values.mean() + 1e-9) / (u_values.mean() + 1e-9)
        records.append([
            ligand, receptor, len(t_values), len(u_values),
            np.log2(fold_change), p_value,
        ])
    if not records:
        return pd.DataFrame()
    results = pd.DataFrame(records, columns=[
        "ligand", "receptor", "n_fov_treated", "n_fov_untreated",
        "log2FC_T_over_U", "pval",
    ])
    results["qval"] = multipletests(results["pval"], method="fdr_bh")[1]
    comparison_columns = ["ligand", "receptor", "log2FC_T_over_U", "pval", "qval"]
    for summary, path in [
        (treated_summary, treated_csv_path), (untreated_summary, untreated_csv_path)
    ]:
        summary.merge(
            results[comparison_columns], on=["ligand", "receptor"], how="left"
        ).to_csv(path, index=False)
    return results
