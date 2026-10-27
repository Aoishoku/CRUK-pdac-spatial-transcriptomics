"""Original IQR filtering rules (cells 86–87); no new statistical methods."""


def outliers_iqr(proximity):
    """Return values strictly outside Q1 − 1.5 IQR and Q3 + 1.5 IQR."""
    lower_quartile = proximity.quantile(0.25)
    upper_quartile = proximity.quantile(0.75)
    iqr = upper_quartile - lower_quartile
    lower = lower_quartile - 1.5 * iqr
    upper = upper_quartile + 1.5 * iqr
    return proximity[(proximity < lower) | (proximity > upper)]


def erase_outliers(proximity_dataframe, mean_or_median):
    """Return a copy filtered within each cell type and treatment group.

    This preserves the original analysis choice; it is not a recommendation
    to remove biological observations. Input row indices must be unique.
    """
    if mean_or_median not in {"mean", "median"}:
        raise ValueError("mean_or_median must be 'mean' or 'median'")
    if not proximity_dataframe.index.is_unique:
        raise ValueError("row indices must be unique before outlier filtering")
    outlier_indices = []
    for cell_type in proximity_dataframe["interacting_cell"].unique():
        group = proximity_dataframe.loc[
            proximity_dataframe["interacting_cell"] == cell_type
        ]
        for condition in ["Treated", "Untreated"]:
            values = group.loc[
                group["Condition"] == condition, f"{mean_or_median}_proximity"
            ]
            if not values.empty:
                outlier_indices.extend(outliers_iqr(values).index.tolist())
    return proximity_dataframe.drop(index=outlier_indices)
