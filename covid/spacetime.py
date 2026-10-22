"""Spatial and space-time statistics: global Moran's I, Getis-Ord Gi*,
Benjamini-Hochberg FDR, Mann-Kendall trend test and an emerging hot spot
classification modelled on ArcGIS Pro's Emerging Hot Spot Analysis.

Gi* follows Getis & Ord (1992) and Ord & Getis (1995). The space-time version
uses the same formula over a space-time cube: a bin's neighbourhood is its own
county plus its queen-contiguous neighbours, in the current time step and the
previous `time_window` steps, and the mean and standard deviation come from
every bin in the cube (as described in Esri's "How Emerging Hot Spot Analysis
works").
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse, stats
from esda.moran import Moran
from libpysal.weights import Queen


# ---------------------------------------------------------------- weights
def queen_weights(gdf):
    """Binary queen contiguity weights (libpysal), in gdf row order."""
    w = Queen.from_dataframe(gdf, use_index=False, silence_warnings=True)
    return w


def binary_matrix_with_self(w) -> sparse.csr_matrix:
    """Sparse binary neighbour matrix including the focal unit (for Gi*)."""
    m, _ = w.full()
    m = (m > 0).astype(float)
    np.fill_diagonal(m, 1.0)
    return sparse.csr_matrix(m)


# ---------------------------------------------------------------- Gi*
def gi_star(x: np.ndarray, W: sparse.csr_matrix) -> np.ndarray:
    """Getis-Ord Gi* z-scores for values x with a binary weights matrix W that
    includes self (Ord & Getis 1995, eq. for G*_i as a standard normal z)."""
    x = np.asarray(x, float)
    n = x.size
    xbar = x.mean()
    s = np.sqrt((x ** 2).mean() - xbar ** 2)
    wsum = np.asarray(W.sum(axis=1)).ravel()
    w2sum = np.asarray(W.multiply(W).sum(axis=1)).ravel()
    num = W @ x - xbar * wsum
    den = s * np.sqrt((n * w2sum - wsum ** 2) / (n - 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(den > 0, num / den, 0.0)
    return z


def spacetime_weights(W: sparse.csr_matrix, n_time: int, time_window: int = 1) -> sparse.csr_matrix:
    """Space-time neighbour matrix for a cube flattened time-major
    (index = t * n_space + i). Neighbours of (i, t): spatial neighbours of i
    (incl. i) at times t, t-1, ..., t-time_window."""
    n = W.shape[0]
    blocks = []
    for lag in range(time_window + 1):
        # place W on the sub-diagonal `lag` (bin at t looks back to t-lag)
        shift = sparse.eye(n_time, n_time, k=-lag, format="csr")
        blocks.append(sparse.kron(shift, W, format="csr"))
    return sum(blocks).tocsr()


def spacetime_gi_star(cube: np.ndarray, W: sparse.csr_matrix, time_window: int = 1) -> np.ndarray:
    """Gi* z-score for every bin of a (time, space) cube."""
    T, n = cube.shape
    Wst = spacetime_weights(W, T, time_window)
    return gi_star(cube.ravel(), Wst).reshape(T, n)


def p_from_z(z):
    return 2 * stats.norm.sf(np.abs(z))


def bh_fdr(p: np.ndarray, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg (1995). Returns (reject, q_values)."""
    p = np.asarray(p, float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    q = np.empty(m)
    q[order] = q_sorted
    return q <= alpha, q


def gi_classes(z: np.ndarray) -> np.ndarray:
    """Esri-style Gi* confidence bins after FDR, coded -3..3
    (+/-3: q<0.01, +/-2: q<0.05, +/-1: q<0.10, 0: not significant)."""
    _, q = bh_fdr(p_from_z(z))
    level = np.select([q < 0.01, q < 0.05, q < 0.10], [3, 2, 1], 0)
    return (np.sign(z) * level).astype(int)


# ---------------------------------------------------------------- Moran's I
def moran_series(table: pd.DataFrame, w, permutations: int = 999, seed: int = 0) -> pd.DataFrame:
    """Global Moran's I for each row (time step) of a time x county table."""
    w.transform = "r"
    rows = []
    for t, row in table.iterrows():
        np.random.seed(seed)
        m = Moran(row.to_numpy(float), w, permutations=permutations)
        rows.append({"period": t, "I": m.I, "z_sim": m.z_sim, "p_sim": m.p_sim})
    w.transform = "o"
    return pd.DataFrame(rows).set_index("period")


# ---------------------------------------------------------------- Mann-Kendall
def mann_kendall(x) -> tuple[float, float, float]:
    """Mann (1945) / Kendall (1975) trend test with tie correction.
    Returns (S, z, two-sided p)."""
    x = np.asarray(x, float)
    n = x.size
    diffs = np.sign(x[None, :] - x[:, None])
    s = np.triu(diffs, k=1).sum()
    _, counts = np.unique(x, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - (counts * (counts - 1) * (2 * counts + 5)).sum()) / 18.0
    if var <= 0:
        return float(s), 0.0, 1.0
    z = (s - np.sign(s)) / np.sqrt(var)  # continuity correction
    return float(s), float(z), float(2 * stats.norm.sf(abs(z)))


# ---------------------------------------------------------------- categories
CATEGORIES = [
    "New Hot Spot", "Consecutive Hot Spot", "Intensifying Hot Spot", "Persistent Hot Spot",
    "Diminishing Hot Spot", "Sporadic Hot Spot", "Oscillating Hot Spot", "Historical Hot Spot",
    "New Cold Spot", "Consecutive Cold Spot", "Intensifying Cold Spot", "Persistent Cold Spot",
    "Diminishing Cold Spot", "Sporadic Cold Spot", "Oscillating Cold Spot", "Historical Cold Spot",
    "No Pattern Detected",
]


def _definitions_one_side(sig, other, z, kind, alpha, share):
    """Every category definition met on one side (hot or cold), read literally
    from Esri's "How Emerging Hot Spot Analysis works". `sig` = significant
    bins of this kind, `other` = significant bins of the opposite kind. For cold
    spots z is negated before calling, so "intensifying" always means the
    clustering is getting stronger."""
    T = sig.size
    frac = sig.sum() / T
    met = []
    _, mz, mp = mann_kendall(z)
    trend_up, trend_down = mp < alpha and mz > 0, mp < alpha and mz < 0
    if sig[-1]:
        run = T if sig.all() else T - np.flatnonzero(~sig)[-1] - 1  # final uninterrupted run
        hot_before_run = sig[: T - run].any()
        if not sig[:-1].any():
            met.append("New")
        if run >= 2 and not hot_before_run and frac < share:
            met.append("Consecutive")
        if frac >= share and trend_up:
            met.append("Intensifying")
        if frac >= share and trend_down:
            met.append("Diminishing")
        if frac < share and hot_before_run and not other.any():
            met.append("Sporadic")
        if frac < share and other[:-1].any():
            met.append("Oscillating")
    elif frac >= share:
        met.append("Historical")
    if frac >= share and not (trend_up or trend_down):
        met.append("Persistent")  # Esri's wording does not mention the final step
    return [f"{m} {kind} Spot" for m in met]


# Esri's page gives the definitions but no rule for a location that meets more
# than one. Precedence used here: categories that describe the final time step
# come first (the purpose of the tool is the current state), the more specific
# "no earlier history" categories (New, Consecutive) before the history-based
# ones (Sporadic, Oscillating), hot before cold, and Historical before a
# Persistent reading that lacks a significant final step.
PRECEDENCE = [f"{c} {k} Spot" for k in ("Hot", "Cold") for c in
              ("New", "Consecutive", "Intensifying", "Diminishing", "Sporadic", "Oscillating")]
PRECEDENCE += ["Historical Hot Spot", "Historical Cold Spot", "Persistent Hot Spot", "Persistent Cold Spot"]


def ehsa_definitions(z, hot, cold, alpha: float = 0.05, share: float = 0.9) -> list[str]:
    """All categories whose definition the location meets, in precedence order."""
    z, hot, cold = np.asarray(z, float), np.asarray(hot, bool), np.asarray(cold, bool)
    met = (_definitions_one_side(hot, cold, z, "Hot", alpha, share)
           + _definitions_one_side(cold, hot, -z, "Cold", alpha, share))
    return sorted(met, key=PRECEDENCE.index)


def classify_ehsa(z: np.ndarray, hot: np.ndarray, cold: np.ndarray,
                  alpha: float = 0.05, share: float = 0.9) -> str:
    """Emerging hot spot category for one location (first definition met in
    PRECEDENCE order, else "No Pattern Detected").

    z    : Gi* z-scores over time (oldest first)
    hot  : bool, bin is a significant hot spot (after FDR)
    cold : bool, bin is a significant cold spot
    """
    met = ehsa_definitions(z, hot, cold, alpha, share)
    return met[0] if met else "No Pattern Detected"


def emerging_hot_spots(cube: np.ndarray, W: sparse.csr_matrix, time_window: int = 1,
                       alpha: float = 0.05, fdr: bool = True) -> dict:
    """Run space-time Gi*, FDR per time step (optional), Mann-Kendall per
    location and classify every location. cube is (time, space)."""
    Z = spacetime_gi_star(cube, W, time_window)
    T, n = Z.shape
    hot = np.zeros_like(Z, bool)
    cold = np.zeros_like(Z, bool)
    for t in range(T):
        p = p_from_z(Z[t])
        reject = bh_fdr(p, alpha)[0] if fdr else p < alpha
        hot[t] = reject & (Z[t] > 0)
        cold[t] = reject & (Z[t] < 0)
    cats, also, trend_z, trend_p = [], [], [], []
    for i in range(n):
        met = ehsa_definitions(Z[:, i], hot[:, i], cold[:, i], alpha)
        cats.append(met[0] if met else "No Pattern Detected")
        also.append(met[1:])
        _, mz, mp = mann_kendall(Z[:, i])
        trend_z.append(mz)
        trend_p.append(mp)
    return {"z": Z, "hot": hot, "cold": cold, "category": np.array(cats), "also_meets": also,
            "trend_z": np.array(trend_z), "trend_p": np.array(trend_p)}
