from __future__ import annotations

import hashlib
from typing import Any

import numpy as np
import scipy.sparse as sp


def _parse_count_frac(count_frac: float) -> float:
    count_frac = float(count_frac)
    if count_frac < 0.0 or count_frac > 1.0:
        raise ValueError(f"count_frac must be in [0, 1], got {count_frac}")
    return count_frac


def _validate_nonneg_integer_counts(values: np.ndarray) -> np.ndarray:
    """Return rounded int64 counts, raising if values are negative or non-integer."""
    if values.size == 0:
        return np.asarray(values, dtype=np.int64)
    if np.any(values < 0):
        raise ValueError("Counts must be non-negative for binomial downsampling")
    rounded = np.round(values)
    if not np.allclose(values, rounded, rtol=0.0, atol=1e-6):
        raise ValueError("Counts must be integers for binomial downsampling")
    return rounded.astype(np.int64)


def binomial_downsample_matrix(
    X: Any,
    count_frac: float,
    *,
    seed_key: str = "",
) -> Any:
    """
    Thin each matrix entry independently as Binomial(x_ij, count_frac).

    Sparse matrices are thinned in-place on their stored non-zeros (CSR).
    Dense matrices are thinned elementwise. Seeding is deterministic given
    seed_key and count_frac.
    """
    count_frac = _parse_count_frac(count_frac)
    if count_frac >= 1.0:
        return X.copy() if sp.issparse(X) else np.array(X, copy=True)

    seed_hex = hashlib.md5(
        f"{seed_key}|count_frac={count_frac}".encode("utf-8")
    ).hexdigest()[:8]
    rng = np.random.default_rng(int(seed_hex, 16))

    if sp.issparse(X):
        out = X.tocsr(copy=True)
        if out.data.size:
            n = _validate_nonneg_integer_counts(out.data)
            thinned = rng.binomial(n, count_frac)
            out.data = thinned.astype(out.dtype, copy=False)
            out.eliminate_zeros()
        return out

    arr = np.asarray(X)
    n = _validate_nonneg_integer_counts(arr.ravel()).reshape(arr.shape)
    thinned = rng.binomial(n, count_frac)
    return thinned.astype(arr.dtype, copy=False)


def binomial_downsample_counts(
    adata: Any,
    count_frac: float,
    *,
    seed_key: str = "",
    copy: bool = True,
):
    """
    Deterministically thin per-cell gene counts in `adata.X`.

    Each entry is sampled as Binomial(x_ij, count_frac). If count_frac >= 1,
    returns adata unchanged (or a copy when copy=True).
    """
    count_frac = _parse_count_frac(count_frac)
    if count_frac >= 1.0:
        return adata.copy() if copy else adata

    out = adata.copy() if copy else adata
    out.X = binomial_downsample_matrix(out.X, count_frac, seed_key=seed_key)
    return out
