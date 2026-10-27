"""Matrix preparation and sample-level heatmap values (original cells 24, 42, 46).

Cell numbers in this project are zero-based JSON notebook indices.
"""

import numpy as np
import pandas as pd
from scipy import sparse


def replace_nan_with_zero(matrix):
    """Return a CSR copy with NaNs replaced by zero, preserving infinities.

    This retains the notebook's imputation rule without densifying sparse data.
    It does not normalize or log-transform expression.
    """
    if sparse.issparse(matrix):
        result = sparse.csr_matrix(matrix, copy=True)
        result.data[np.isnan(result.data)] = 0
        return result
    values = np.array(matrix, copy=True)
    values[np.isnan(values)] = 0
    return sparse.csr_matrix(values)


def sample_expression_zscores(matrix, sample_ids, gene_names):
    """Average cells per sample, z-score each gene, and return genes × samples.

    Preserve pandas sample standard deviation (ddof=1). Constant genes and
    single-sample inputs retain undefined (NaN) z-scores; no values are filled.
    """
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    means = pd.DataFrame(
        values, index=pd.Index(sample_ids).astype(str), columns=gene_names
    ).groupby(level=0).mean()
    return means.subtract(means.mean(0), axis=1).divide(means.std(0), axis=1).T
