"""Rates per 100,000 and Empirical Bayes smoothing.

Why smooth: Texas has counties of about 50 people (Loving) and of 4.7 million
(Harris). One case in Loving is about 1,900 per 100,000; one case in Harris
barely registers. Raw rates for tiny counties swing wildly and would dominate
any map or hot spot test. Global Empirical Bayes (Marshall 1991) shrinks each
county's rate toward the statewide rate, strongly when the county is small and
hardly at all when it is large, so the map shows signal rather than noise.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from esda.smoothing import Empirical_Bayes

PER = 100_000


def per_100k(count, pop):
    """Count per 100,000 population (works on scalars, arrays and pandas)."""
    return count / pop * PER


def eb_rate(events: pd.Series, pop: pd.Series) -> pd.Series:
    """Global Empirical Bayes smoothed rate per 100,000 (Marshall 1991)."""
    e = events.to_numpy(float).reshape(-1, 1)
    b = pop.reindex(events.index).to_numpy(float).reshape(-1, 1)
    if e.sum() == 0:
        return pd.Series(0.0, index=events.index)
    r = Empirical_Bayes(e, b).r.ravel()
    return pd.Series(r * PER, index=events.index)


def eb_table(counts: pd.DataFrame, pop: pd.Series) -> pd.DataFrame:
    """EB rate per 100,000 for every row (time step) of a time x county table."""
    return counts.apply(lambda row: eb_rate(row, pop), axis=1)


def marshall_eb(events, pop):
    """Plain-numpy version of Marshall's global EB estimator, used in tests to
    check the esda result."""
    e = np.asarray(events, float)
    n = np.asarray(pop, float)
    r = e / n
    m = e.sum() / n.sum()
    s2 = (n * (r - m) ** 2).sum() / n.sum()
    a = max(s2 - m / n.mean(), 0.0)
    w = a / (a + m / n) if a > 0 else np.zeros_like(n)
    return (w * r + (1 - w) * m) * PER


# Named periods. Variant names describe which lineage dominated US sequencing
# for most of each window; the boundaries are round calendar months chosen to
# bracket each Texas wave, not exact variant-takeover dates.
WAVES = {
    "Summer 2020": ("2020-06-01", "2020-09-30"),
    "Winter 2020-21": ("2020-10-01", "2021-03-31"),
    "Delta": ("2021-07-01", "2021-11-30"),
    "Omicron BA.1": ("2021-12-01", "2022-03-31"),
    "Omicron BA.5 and later": ("2022-06-01", "2023-02-28"),
}


def window_sum(table: pd.DataFrame, start: str, end: str) -> pd.Series:
    """Sum a time x county table over [start, end] (period end dates)."""
    return table.loc[(table.index >= start) & (table.index <= end)].sum()


def implausible_resident_rate(cases_total: pd.Series, pop: pd.Series, threshold: float = 1.0) -> pd.Series:
    """True where cumulative reported cases exceed `threshold` x the resident
    population. A resident-based rate cannot be interpreted there: more cases
    than residents usually means the cases include people who do not live in
    the county (for example a temporary workforce), so the denominator is wrong."""
    ratio = cases_total / pop.reindex(cases_total.index)
    return ratio > threshold
