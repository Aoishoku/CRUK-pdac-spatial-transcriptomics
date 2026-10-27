# Technical notes

These notes describe the existing numerical implementation and its main limitations. Analytical choices and known bugs are retained; the results should be treated as exploratory.

## Numerical implementation

- **Expression:** NaNs are replaced with zero, cells are normalized to 10,000 counts, and `log1p` is applied. Sparse matrices stay sparse; infinities remain unchanged. Heatmaps use sample-average expression and gene-wise z-scores with sample standard deviation (`ddof=1`). Constant genes retain undefined scores.
- **Spatial proximity:** Pairwise Euclidean distances are calculated within each sample and field of view (FOV). The main score sums `exp(-distance_px × 0.18 / 50)` over target cells. The 50 µm parameter is a decay length, not a cutoff. Missing targets score zero; identical populations include self-distances. Distance matrices are dense and can be memory-intensive. The 0.18 µm/px calibration needs confirmation from acquisition metadata.
- **Interactions:** Candidates combine Connectome scores ≥ 1 and CellPhoneDB p-values < 0.05. SCOTIA likelihoods are averaged per gene pair/FOV before two-sided Mann–Whitney tests, requiring two detected FOVs per condition. Log2 fold changes use a `1e-9` pseudocount; Benjamini–Hochberg adjustment is applied within each cell-type comparison. Descriptive summaries average interaction rows, so their weighting differs from the tests.

## Scientific limitations

- **Independent observations:** DGE treats cells as observations, although cells from one sample or donor are correlated. Proximity tests aggregate by `sample_id`; sample independence and any pairing need confirmation. SCOTIA treats FOVs as independent, but multiple FOVs from one donor are not biological replicates. These dependencies can overstate evidence ([Scanpy guidance](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html)). [Mann–Whitney tests](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) compare distributions, not necessarily medians.
- **Signature selection:** Genes are selected and scored on the same data. Separate 90th-percentile thresholds define relative enrichment within each treatment, not a common absolute identity; ties can select more than 10%. The “Antigen signature presenting” label lacks independent validation. g:Profiler is mentioned, but enrichment results and a background gene list were not supplied.
- **Preprocessing:** Prior normalization and the meaning of missing values were not verified. NaN-to-zero replacement and possible repeated normalization remain assumptions. DGE and signature scoring omit `use_raw`; if `.raw` exists, Scanpy may use a different matrix from normalized `.X`. Its contents were not inspected.
- **Filtering and multiplicity:** Group-wise IQR removal can change effects and uncertainty in small cohorts. Shapiro/Levene pretesting adds a data-dependent choice; failure to reject normality does not establish it. Proximity p-values are unadjusted across cell types and summary measures. Connectome's cutoff is a score filter; CellPhoneDB uses raw p-values. Later BH adjustment does not control error across the entire selection pipeline.
- **Spatial interpretation:** Scores reflect partner abundance and geometry as well as proximity, without normalization by abundance or tissue area. Cross-FOV neighbors are omitted, and tissue boundaries can affect comparisons.
- **Interaction interpretation:** Combining gene pairs may pool opposite source/target directions. Only detected interactions enter tests; missing detections are not filled with zero, and low-FOV pairs are skipped. Detection differences can influence comparisons, while near-zero likelihoods produce pseudocount-sensitive fold changes. LIANA/SCOTIA predictions do not demonstrate signaling or causal treatment effects.

## Known result-affecting issues

- **Mean-proximity tests reuse median-based assumption checks**, potentially selecting an inappropriate test. Groups with fewer than three values are skipped; constant/missing data are not handled explicitly, and NaN test p-values can be treated as passing.
- **Every non-`T` sample identifier is labeled untreated**, without validating missing or unexpected identifiers. Check against study metadata, including `treatment_status`.
- **Empty or degenerate inputs can interrupt processing:** empty proximity records cannot be concatenated; empty gene lists, constant-gene z-scores or fewer than two samples per treatment can break clustering. Categorical groupings may introduce unobserved combinations with NaN summaries; pandas defaults are version-sensitive.
- **Parsing and incomplete runs need attention:** SCOTIA parsing assumes the first underscore separates genes and filename stems uniquely identify FOVs across samples. The loader reports malformed nonempty files, but notebook exception handlers can still skip cell-type pairs. If no pairs qualify, no summary is written, so older output files can remain stale.

## Reproducibility

The cleaned notebook has not been executed after refactoring. The original notebook remains local and excluded from Git; selected figures preserve its saved outputs. The dataset is also excluded; its published source is documented in the [README](../README.md#data-and-optional-notebook-setup).

The original kernel recorded Python 3.11.5, but exact package versions, random-state settings and the L-R resource snapshot are unavailable. Dependency ranges are not a historical lockfile; SCOTIA is unpinned and full notebook compatibility is unverified. The notebook retains the `np.Inf = np.inf` workaround. SciPy's distance API is deprecated, with its announced removal excluded by the dependency range. LIANA repeats full-data analysis per target cell type, with 1,000 permutations per call, so full execution can be expensive.

The 21 public tests use synthetic inputs and fixed reference values without the original notebook or dataset. They check selected numerical behavior, output handling and execution guards, not full-data equivalence or biological validity. No analysis was rerun to prepare these notes.
