"""The notebook's exponential spatial proximity score (original cell 62)."""

import numpy as np
from scipy.spatial import distance_matrix


def exponential_model(
    adata, ct1, ct2, ct1_obs="annotation_majortypes",
    ct2_obs="annotation_majortypes", radius=20,
):
    """Sum exp(-distance_px * 0.18 / radius) within each sample/FOV.

    Accept an object with a pandas ``obs`` table, including AnnData. Return
    one array per population, in sorted sample/FOV order and original row
    order within each FOV, matching the notebook. Missing partners score zero.
    The radius is a decay length in microns, not a hard distance cutoff.
    The dataset-specific calibration remains 0.18 microns per pixel.

    This still builds a full pairwise distance matrix for each FOV. Identical
    populations still include self-distances, as in the original calculation.
    """
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be finite and positive")
    obs = adata.obs
    coordinates = ["CenterX_local_px", "CenterY_local_px"]
    density1, density2 = [], []
    for sample_id in np.unique(obs["sample_id"]):
        sample = obs.loc[obs["sample_id"] == sample_id]
        for fov in np.unique(sample["fov"]):
            subset = sample.loc[sample["fov"] == fov]
            first = subset.loc[subset[ct1_obs] == ct1, coordinates]
            second = subset.loc[subset[ct2_obs] == ct2, coordinates]
            if len(first) and len(second):
                distances = distance_matrix(first, second)
                weights = np.exp(-(distances * 0.18) / radius)
                density1.append(np.sum(weights, axis=1))
                density2.append(np.sum(weights, axis=0))
            elif len(first):
                density1.append(np.zeros(len(first)))
            elif len(second):
                density2.append(np.zeros(len(second)))
    return (
        np.concatenate(density1) if density1 else np.array([]),
        np.concatenate(density2) if density2 else np.array([]),
    )
