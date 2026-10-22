"""Fast tests on synthetic data (no downloads)."""
import os

os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from covid import data, rates, spacetime as st


# ------------------------------------------------------------ cumulative -> new
def test_cumulative_to_new_simple():
    cum = pd.Series([0, 2, 5, 5, 9], dtype=float)
    new, n_neg = data.cumulative_to_new(cum)
    assert new.tolist() == [0, 2, 3, 0, 4]
    assert n_neg == 0


def test_cumulative_to_new_downward_revision():
    # 10 -> 8 is a correction: the 2 removed cases were never real
    cum = pd.Series([3, 10, 8, 12], dtype=float)
    new, n_neg = data.cumulative_to_new(cum)
    assert n_neg == 1
    assert (new >= 0).all()
    assert new.sum() == cum.iloc[-1]          # totals are preserved
    assert new.tolist() == [3, 5, 0, 4]


def test_cumulative_to_new_dataframe_and_weekly():
    idx = pd.date_range("2021-01-04", periods=14, freq="D")  # Mon .. Sun (2 weeks)
    cum = pd.DataFrame({"a": np.arange(14) * 2.0, "b": [0] * 7 + [5] * 6 + [4]}, index=idx)
    wk, n_neg = data.to_period_new(cum, "W-SUN")
    assert n_neg == 0  # sampling happens after correction; count comes from the weekly diff
    assert wk["a"].tolist() == [12.0, 14.0]
    assert wk["b"].tolist() == [0.0, 4.0]
    assert wk.sum().sum() == cum.iloc[-1].sum()


# ------------------------------------------------------------ rates
def test_per_100k():
    assert rates.per_100k(50, 200_000) == pytest.approx(25.0)
    s = rates.per_100k(pd.Series([1, 10]), pd.Series([1000, 1_000_000]))
    assert s.tolist() == pytest.approx([100.0, 1.0])


def test_eb_matches_marshall_and_shrinks_small_counties():
    events = pd.Series([1, 50, 400, 3, 0], index=list("abcde"), dtype=float)
    pop = pd.Series([60, 20_000, 150_000, 5_000, 800], index=list("abcde"), dtype=float)
    eb = rates.eb_rate(events, pop)
    assert eb.to_numpy() == pytest.approx(rates.marshall_eb(events, pop), rel=1e-6)
    crude = rates.per_100k(events, pop)
    state = events.sum() / pop.sum() * 1e5
    # tiny county "a" moves most of the way toward the state rate
    assert abs(eb["a"] - state) < abs(crude["a"] - state) * 0.2
    # big county "c" barely moves
    assert eb["c"] == pytest.approx(crude["c"], rel=0.05)


# ------------------------------------------------------------ Gi*, FDR
def test_gi_star_detects_cluster_on_a_line():
    n = 30
    W = sparse.diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(n, n), format="csr")
    x = np.ones(n)
    x[10:14] = 10
    z = st.gi_star(x, W)
    assert z[11] > 3 and z[12] > 3
    assert abs(z[25]) < 1


def test_spacetime_weights_look_back_only():
    W = sparse.identity(2, format="csr")
    Wst = st.spacetime_weights(W, n_time=3, time_window=1).toarray()
    # bin (t=1, i=0) -> index 2 should see itself and (t=0, i=0) -> index 0, not t=2
    assert Wst[2, 2] == 1 and Wst[2, 0] == 1 and Wst[2, 4] == 0


def test_bh_fdr():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    reject, q = st.bh_fdr(p, 0.05)
    assert reject.tolist() == [True, True] + [False] * 8
    assert q[0] == pytest.approx(0.01)


# ------------------------------------------------------------ Mann-Kendall
def test_mann_kendall_known_values():
    s, z, p = st.mann_kendall([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    assert s == 45 and z > 0 and p < 0.001
    s, z, p = st.mann_kendall([5, 5, 5, 5])
    assert p == 1.0
    _, z, _ = st.mann_kendall(list(range(10, 0, -1)))
    assert z < 0


# ------------------------------------------------------------ categories
T = 20
F = np.zeros(T, bool)


def mk(hot_idx=(), cold_idx=()):
    h, c = F.copy(), F.copy()
    h[list(hot_idx)] = True
    c[list(cold_idx)] = True
    return h, c


@pytest.mark.parametrize("hot_idx,cold_idx,z,expected", [
    ([], [], np.zeros(T), "No Pattern Detected"),
    ([19], [], np.zeros(T), "New Hot Spot"),
    ([17, 18, 19], [], np.zeros(T), "Consecutive Hot Spot"),
    (range(T), [], np.linspace(2, 6, T), "Intensifying Hot Spot"),
    (range(T), [], np.linspace(6, 2, T), "Diminishing Hot Spot"),
    (range(T), [], np.tile([3.0, 3.5], T // 2), "Persistent Hot Spot"),
    ([3, 8, 19], [], np.zeros(T), "Sporadic Hot Spot"),
    ([3, 19], [10], np.zeros(T), "Oscillating Hot Spot"),
    (range(18), [], np.zeros(T), "Historical Hot Spot"),
    ([], [19], np.zeros(T), "New Cold Spot"),
    ([], range(T), -np.linspace(2, 6, T), "Intensifying Cold Spot"),
    ([], [16, 17, 18, 19], np.zeros(T), "Consecutive Cold Spot"),
    ([5], [2, 19], np.zeros(T), "Oscillating Cold Spot"),
])
def test_classify_ehsa(hot_idx, cold_idx, z, expected):
    h, c = mk(hot_idx, cold_idx)
    assert st.classify_ehsa(z, h, c) == expected


def test_consecutive_requires_no_earlier_hot():
    h, c = mk([2, 17, 18, 19])
    assert st.classify_ehsa(np.zeros(T), h, c) == "Sporadic Hot Spot"


def test_overlapping_definitions_are_reported():
    # hot only in the final step, cold earlier: meets New Hot and Oscillating Hot
    h, c = mk([19], [4])
    met = st.ehsa_definitions(np.zeros(T), h, c)
    assert met == ["New Hot Spot", "Oscillating Hot Spot"]
    assert st.classify_ehsa(np.zeros(T), h, c) == "New Hot Spot"


def test_historical_beats_persistent_without_final_step():
    h, c = mk(range(19))
    met = st.ehsa_definitions(np.tile([3.0, 3.5], T // 2), h, c)
    assert met[0] == "Historical Hot Spot" and "Persistent Hot Spot" in met


def test_precedence_lists_every_category():
    assert set(st.PRECEDENCE) | {"No Pattern Detected"} == set(st.CATEGORIES)


# ------------------------------------------------------------ exclusion rule
def test_implausible_resident_rate():
    cases = pd.Series({"loving": 399.0, "big": 300_000.0, "edge": 100.0})
    pop = pd.Series({"loving": 54.0, "big": 1_000_000.0, "edge": 100.0})
    ex = rates.implausible_resident_rate(cases, pop)
    assert ex.to_dict() == {"loving": True, "big": False, "edge": False}  # equal is not "exceeds"
    assert rates.implausible_resident_rate(cases, pop, threshold=0.25)["big"]
