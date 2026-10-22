"""COVID-19 in Texas counties, 2020-2023: where and when did hot spots emerge?

Downloads public data, computes rates, runs the spatial and space-time
statistics and writes figures/, results/ and docs/index.html.

Primary analysis: counties whose cumulative reported cases exceed their
resident population are excluded from the rate-based statistics (their
resident rates cannot be interpreted) but kept in statewide totals. The same
statistics on all counties are reported as a sensitivity analysis.

    python run.py                      # Texas, March 2020 - February 2023
    python run.py --state 48 --start 2020-03 --end 2023-02
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")  # see covid/__init__.py

import numpy as np
import pandas as pd
from esda.moran import Moran
from libpysal.weights import KNN
from libpysal.weights.util import attach_islands

from covid import data, maps, rates, spacetime as st

ROOT = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--state", default="48", help="two-digit state FIPS (default 48, Texas)")
    p.add_argument("--start", default="2020-03", help="first month, YYYY-MM (default 2020-03)")
    p.add_argument("--end", default="2023-02",
                   help="last complete month, YYYY-MM (default 2023-02; NYT data stop 2023-03-23)")
    p.add_argument("--time-window", type=int, default=1,
                   help="previous months in each space-time neighbourhood (default 1)")
    p.add_argument("--permutations", type=int, default=999, help="Moran's I permutations")
    p.add_argument("--alpha", type=float, default=0.05, help="FDR and Mann-Kendall significance level")
    p.add_argument("--max-cases-per-resident", type=float, default=1.0,
                   help="exclude counties whose cumulative cases / population exceed this (default 1)")
    p.add_argument("--raw", type=Path, default=ROOT / "data" / "raw")
    return p.parse_args()


def r(x, nd=1):
    return float(round(float(x), nd))


def cat_counts(c):
    vc = pd.Series(c).value_counts()
    return {k: int(vc.get(k, 0)) for k in st.CATEGORIES if vc.get(k, 0) > 0}


def analyze(gdf, monthly, pop, keep, waves, a):
    """All rate-based statistics on the counties where `keep` is True.
    Weights are rebuilt on that subset; any county left without a contiguous
    neighbour is joined to its nearest neighbour (k = 1)."""
    sub = gdf[keep].reset_index(drop=True)
    f = sub["fips"].tolist()
    p = pop.reindex(f)
    mon = {k: v[f] for k, v in monthly.items()}

    w = st.queen_weights(sub)
    islands = [sub["name"].iloc[i] for i in w.islands]
    if w.islands:
        w = attach_islands(w, KNN.from_dataframe(sub, k=1))

    eb_month = {k: rates.eb_table(v, p) for k, v in mon.items()}
    total = {k: v.sum() for k, v in mon.items()}
    eb_total = {k: rates.eb_rate(total[k], p) for k in total}

    moran = {k: st.moran_series(eb_month[k], w, a.permutations) for k in eb_month}
    m_full = {}
    for k in ("cases", "deaths"):
        w.transform = "r"
        np.random.seed(0)
        mm = Moran(eb_total[k].to_numpy(), w, permutations=a.permutations)
        w.transform = "o"
        m_full[k] = {"I": r(mm.I, 3), "p_sim": r(mm.p_sim, 4)}

    W = st.binary_matrix_with_self(w)
    rr_wave, gi_wave, gi_counts, gi_counts_d, ref_wave = {}, {}, {}, {}, {}
    for name, (ws, we) in waves.items():
        ev = rates.window_sum(mon["cases"], ws, we)
        eb = rates.eb_rate(ev, p)
        ref_wave[name] = rates.per_100k(ev.sum(), p.sum())
        rr_wave[name] = eb / ref_wave[name]
        cls = st.gi_classes(st.gi_star(eb.to_numpy(), W))
        gi_wave[name] = pd.Series(cls, index=f)
        gi_counts[name] = {"hot": int((cls >= 2).sum()), "cold": int((cls <= -2).sum())}
        evd = rates.window_sum(mon["deaths"], ws, we)
        clsd = st.gi_classes(st.gi_star(rates.eb_rate(evd, p).to_numpy(), W))
        gi_counts_d[name] = {"hot": int((clsd >= 2).sum()), "cold": int((clsd <= -2).sum())}

    # space-time cube: smoothed monthly rate / pooled rate of the analysed counties that month
    ref_month = rates.per_100k(mon["cases"].sum(axis=1), p.sum())
    rel = eb_month["cases"].div(ref_month, axis=0)
    ehsa = st.emerging_hot_spots(rel.to_numpy(), W, a.time_window, a.alpha)
    ref_month_d = rates.per_100k(mon["deaths"].sum(axis=1), p.sum())
    sens = {
        "raw_rate_input_counts": cat_counts(
            st.emerging_hot_spots(eb_month["cases"].to_numpy(), W, a.time_window, a.alpha)["category"]),
        "no_time_window_counts": cat_counts(st.emerging_hot_spots(rel.to_numpy(), W, 0, a.alpha)["category"]),
        "no_fdr_counts": cat_counts(
            st.emerging_hot_spots(rel.to_numpy(), W, a.time_window, a.alpha, fdr=False)["category"]),
    }
    ehsa_d = st.emerging_hot_spots(eb_month["deaths"].div(ref_month_d, axis=0).to_numpy(), W,
                                   a.time_window, a.alpha)
    sens["deaths_counts"] = cat_counts(ehsa_d["category"])

    month_ix = rel.index
    hot, cold = ehsa["hot"], ehsa["cold"]
    first_hot = pd.Series([month_ix[hot[:, i]].min() if hot[:, i].any() else pd.NaT
                           for i in range(len(f))], index=f)

    def wave_of(t):
        if pd.isna(t):
            return "Never a hot spot"
        for name, (ws, we) in waves.items():
            if pd.Timestamp(ws) <= t <= pd.Timestamp(we):
                return name
        return "Between waves"

    return {
        "fips": f, "islands": islands, "eb_month": eb_month, "eb_total": eb_total, "total": total,
        "moran": moran, "moran_full": m_full, "rr_wave": pd.DataFrame(rr_wave), "gi_wave": pd.DataFrame(gi_wave),
        "gi_counts": gi_counts, "gi_counts_deaths": gi_counts_d, "ref_wave": ref_wave, "rel": rel,
        "ehsa": ehsa, "ehsa_d": ehsa_d, "sens": sens,
        "category": pd.Series(ehsa["category"], index=f),
        "months_hot": pd.Series(hot.sum(axis=0), index=f), "months_cold": pd.Series(cold.sum(axis=0), index=f),
        "first_hot": first_hot, "first_wave": first_hot.map(wave_of),
        "hot_pm": pd.Series(hot.sum(axis=1), index=month_ix), "cold_pm": pd.Series(cold.sum(axis=1), index=month_ix),
    }


def moran_summary(m):
    i = m["I"]
    return {"min": r(i.min(), 3), "min_month": str(i.idxmin())[:7],
            "max": r(i.max(), 3), "max_month": str(i.idxmax())[:7],
            "median": r(i.median(), 3),
            "months_significant_p05": int((m["p_sim"] < 0.05).sum()), "months": int(len(m))}


def headline(res, gdf_names):
    """Comparable summary used for the primary vs all-county table."""
    mh = res["months_hot"]
    return {
        "n_counties_analysed": len(res["fips"]),
        "gi_hot_counties_by_wave": {k: v["hot"] for k, v in res["gi_counts"].items()},
        "gi_cold_counties_by_wave": {k: v["cold"] for k, v in res["gi_counts"].items()},
        "counties_hot_at_least_once": int((mh > 0).sum()),
        "counties_never_hot": int((mh == 0).sum()),
        "counties_hot_9plus_months": int((mh >= 9).sum()),
        "hot_counties_peak_month": f"{res['hot_pm'].idxmax():%Y-%m}",
        "hot_counties_peak_count": int(res["hot_pm"].max()),
        "ehsa_counts": cat_counts(res["category"]),
        "morans_i_cases_range": [r(res["moran"]["cases"]["I"].min(), 3), r(res["moran"]["cases"]["I"].max(), 3)],
        "morans_i_cases_median": r(res["moran"]["cases"]["I"].median(), 3),
        "top_hot_counties": [{"county": gdf_names[fp], "months_hot": int(v)}
                             for fp, v in mh.sort_values(ascending=False).head(5).items()],
    }


def main():
    a = parse_args()
    t0 = time.time()
    fig_dir, res_dir, doc_dir = ROOT / "figures", ROOT / "results", ROOT / "docs"
    for d in (fig_dir, res_dir, doc_dir):
        d.mkdir(exist_ok=True)

    print("1/6 data")
    data.download_all(a.raw)
    nyt, nyt_info = data.load_nyt(a.raw, a.state)
    pop = data.load_population(a.raw, a.state)
    gdf = data.load_counties(a.raw, a.state)
    rucc = data.load_rucc(a.raw, a.state).reindex(gdf["fips"])
    metro = data.metro_label(rucc)
    fips = gdf["fips"].tolist()
    names = dict(zip(gdf["fips"], gdf["name"]))
    pop = pop.reindex(fips)
    assert pop.notna().all(), "population missing for some counties"
    state_pop = pop.sum()

    cum = {k: data.daily_cumulative(nyt, k).reindex(columns=fips, fill_value=0.0) for k in ("cases", "deaths")}
    _, n_rev_cases = data.cumulative_to_new(cum["cases"])
    _, n_rev_deaths = data.cumulative_to_new(cum["deaths"])

    start = pd.Period(a.start, "M").start_time
    end = pd.Period(a.end, "M").end_time.normalize()
    period_label = f"{start:%B %Y} - {end:%B %Y}"

    weekly = {k: data.to_period_new(cum[k], "W-SUN")[0] for k in cum}
    monthly = {k: data.to_period_new(cum[k], "ME")[0] for k in cum}
    monthly = {k: v.loc[(v.index >= start) & (v.index <= end)] for k, v in monthly.items()}
    wk_mask = (weekly["cases"].index >= start) & (weekly["cases"].index <= end + pd.Timedelta(days=6))
    state_wk = {k: rates.per_100k(v.loc[wk_mask].sum(axis=1), state_pop) for k, v in weekly.items()}

    # ---- exclusion rule: reported cases in the analysis period vs resident population.
    # (The whole NYT record is not used: Delta County, TX gained 7,842 cases in one day on
    # 2023-03-15, after the analysis period, which looks like a reporting error.)
    final_cases = monthly["cases"].sum()
    excluded = rates.implausible_resident_rate(final_cases, pop, a.max_cases_per_resident).reindex(fips)
    ex_mask = excluded.to_numpy()
    keep = ~ex_mask
    excluded_info = [{"county": names[fp], "fips": fp, "population_2021": int(pop[fp]),
                      "cases_in_period": int(final_cases[fp]),
                      "cases_per_resident": r(final_cases[fp] / pop[fp], 2)}
                     for fp in excluded[excluded].index]
    ratio = (final_cases / pop).sort_values(ascending=False)
    next_highest = [{"county": names[fp], "cases_per_resident": r(v, 3)} for fp, v in ratio.iloc[
        len(excluded_info):len(excluded_info) + 3].items()]
    print(f"   excluded: {[e['county'] for e in excluded_info]}")

    print("2/6 statewide wave timeline (all counties)")
    waves = rates.WAVES
    peaks = {}
    for name, (ws, we) in waves.items():
        pk = {}
        for k in ("cases", "deaths"):
            s = state_wk[k].loc[(state_wk[k].index >= ws) & (state_wk[k].index <= we)]
            pk[f"{k}_peak_week"] = str(s.idxmax().date())
            pk[f"{k}_peak_per_100k"] = r(s.max(), 1 if k == "cases" else 2)
            pk[f"{k}_total_per_100k"] = r(rates.per_100k(rates.window_sum(monthly[k], ws, we).sum(), state_pop))
        peaks[name] = pk
    total_all = {k: v.sum() for k, v in monthly.items()}
    crude_total = {k: rates.per_100k(total_all[k], pop) for k in total_all}

    print("3/6 primary analysis (excluded counties removed)")
    prim = analyze(gdf, monthly, pop, keep, waves, a)
    print("4/6 sensitivity analysis (all counties)")
    allc = analyze(gdf, monthly, pop, np.ones(len(fips), bool), waves, a)
    for k, m in prim["moran"].items():
        m.to_csv(res_dir / f"morans_i_monthly_{k}.csv", float_format="%.4f")

    cats = prim["category"].reindex(fips).fillna("Excluded")
    months_hot = prim["months_hot"].reindex(fips)
    first_wave = prim["first_wave"].reindex(fips).fillna("Excluded")
    ehsa = prim["ehsa"]
    also = pd.Series([";".join(m) for m in ehsa["also_meets"]], index=prim["fips"]).reindex(fips)

    print("5/6 metro / nonmetro (all counties, pooled)")
    ur_rows = {}
    for name, (ws, we) in list(waves.items()) + [("Full period", (str(start.date()), str(end.date())))]:
        row = {}
        for k in ("cases", "deaths"):
            ev = rates.window_sum(monthly[k], ws, we)
            for grp in ("Metro", "Nonmetro"):
                msk = (metro == grp).to_numpy()
                row[(k, grp)] = rates.per_100k(ev[msk].sum(), pop[msk].sum())
        ur_rows[name] = row
    ur = pd.DataFrame(ur_rows).T
    ur.columns = pd.MultiIndex.from_tuples(ur.columns)
    m_keep = metro[keep]
    cats_by_metro = pd.crosstab(prim["category"].to_numpy(), m_keep.to_numpy())

    print("6/6 outputs")
    summary = pd.DataFrame({
        "fips": fips, "county": gdf["name"].to_numpy(), "population_2021": pop.to_numpy().astype(int),
        "rucc_2023": rucc.to_numpy(), "metro": metro.to_numpy(),
        "excluded_from_rate_statistics": ex_mask,
        "cases": total_all["cases"].to_numpy().astype(int), "deaths": total_all["deaths"].to_numpy().astype(int),
        "case_rate_crude": crude_total["cases"].round(1).to_numpy(),
        "case_rate_eb": prim["eb_total"]["cases"].reindex(fips).round(1).to_numpy(),
        "death_rate_crude": crude_total["deaths"].round(1).to_numpy(),
        "death_rate_eb": prim["eb_total"]["deaths"].reindex(fips).round(1).to_numpy(),
        "ehsa_category": cats.to_numpy(),
        "ehsa_trend_z": pd.Series(ehsa["trend_z"], index=prim["fips"]).reindex(fips).round(2).to_numpy(),
        "ehsa_trend_p": pd.Series(ehsa["trend_p"], index=prim["fips"]).reindex(fips).round(4).to_numpy(),
        "months_hot": months_hot.to_numpy(), "months_cold": prim["months_cold"].reindex(fips).to_numpy(),
        "first_hot_month": [f"{t:%Y-%m}" if pd.notna(t) else "" for t in prim["first_hot"].reindex(fips)],
        "first_hot_wave": first_wave.to_numpy(),
        "ehsa_also_meets": also.fillna("").to_numpy(),
        "ehsa_multiple_definitions": (also.fillna("") != "").to_numpy(),
        "ehsa_category_deaths": pd.Series(prim["ehsa_d"]["category"], index=prim["fips"]).reindex(fips)
        .fillna("Excluded").to_numpy(),
        "ehsa_category_all_counties": allc["category"].reindex(fips).to_numpy(),
        "months_hot_all_counties": allc["months_hot"].reindex(fips).to_numpy(),
    })
    for name in waves:
        summary[f"gi_class_{name}"] = prim["gi_wave"][name].reindex(fips).to_numpy()
    summary.to_csv(res_dir / "county_summary.csv", index=False)
    pd.DataFrame({"cases_per_100k": state_wk["cases"].round(2),
                  "deaths_per_100k": state_wk["deaths"].round(3)}).to_csv(
        res_dir / "statewide_weekly.csv", index_label="week_ending")

    # ---------------------------------------------------------- metrics
    an = summary[~summary["excluded_from_rate_statistics"]]
    smallest = an.loc[an["population_2021"].idxmin()]
    mh = prim["months_hot"]
    top_hot = summary[summary["ehsa_category"].str.contains("Hot")] \
        .sort_values("months_hot", ascending=False)[["county", "ehsa_category", "months_hot"]].head(10)
    metrics = {
        "question": "COVID-19 in Texas counties, 2020-2023: where and when did hot spots emerge?",
        "state_fips": a.state, "n_counties": len(fips), "period": period_label,
        "n_months": int(len(prim["rel"])), "population_2021": int(state_pop),
        "exclusion_rule": {
            "rule": f"reported cases in the analysis period > {a.max_cases_per_resident:g} x 2021 "
                    "resident population",
            "reason": ("a resident-based rate is not interpretable; a likely (unverified) explanation is that "
                       "cases among non-resident workers were counted in the county"),
            "excluded_counties": excluded_info,
            "next_highest_cases_per_resident": next_highest,
            "n_counties_analysed": int(keep.sum()),
            "islands_after_exclusion": prim["islands"],
            "kept_in": "statewide totals, wave timeline, metro/nonmetro pooled rates",
            "excluded_from": "Empirical Bayes rate maps, Moran's I, Gi*, space-time cube, emerging hot spots",
            "reference_rate": "pooled rate of the analysed counties (same period)",
        },
        "data_notes": {
            **nyt_info,
            "downward_revisions_daily_cases": n_rev_cases,
            "downward_revisions_daily_deaths": n_rev_deaths,
            "cases_in_period": int(total_all["cases"].sum()), "deaths_in_period": int(total_all["deaths"].sum()),
        },
        "statewide": {
            "cases_per_100k_period": r(rates.per_100k(total_all["cases"].sum(), state_pop)),
            "deaths_per_100k_period": r(rates.per_100k(total_all["deaths"].sum(), state_pop)),
        },
        "eb_smoothing_example": {
            "smallest_analysed_county": smallest["county"],
            "population_2021": int(smallest["population_2021"]), "deaths": int(smallest["deaths"]),
            "death_rate_crude": float(smallest["death_rate_crude"]), "death_rate_eb": float(smallest["death_rate_eb"]),
            "county_death_rate_crude_range": [float(an["death_rate_crude"].min()), float(an["death_rate_crude"].max())],
            "county_death_rate_eb_range": [float(an["death_rate_eb"].min()), float(an["death_rate_eb"].max())],
        },
        "waves": {n: {"start": s, "end": e, **peaks[n]} for n, (s, e) in waves.items()},
        "morans_i_monthly": {k: moran_summary(m) for k, m in prim["moran"].items()},
        "morans_i_full_period": prim["moran_full"],
        "gi_star_by_wave_cases": prim["gi_counts"], "gi_star_by_wave_deaths": prim["gi_counts_deaths"],
        "gi_star_hot_counties_by_wave": {n: sorted(names[fp] for fp, v in prim["gi_wave"][n].items() if v >= 2)
                                         for n in waves},
        "gi_note": "hot/cold counts use FDR-adjusted q < 0.05",
        "emerging_hot_spots": {
            "input": "monthly EB case rate divided by the pooled rate of the analysed counties that month",
            "time_window_months": a.time_window, "alpha": a.alpha,
            "counts": cat_counts(prim["category"]),
            "counties_by_category": {c: sorted(names[fp] for fp, v in prim["category"].items() if v == c)
                                     for c in cat_counts(prim["category"]) if c != "No Pattern Detected"},
            "precedence": st.PRECEDENCE,
            "precedence_note": ("Esri's documentation defines the categories but gives no rule for a "
                                "location meeting more than one; the first definition met in "
                                "'precedence' is assigned and the others are listed in "
                                "county_summary.csv (ehsa_also_meets)."),
            "counties_meeting_multiple_definitions": [
                {"county": names[fp], "category": prim["category"][fp], "also_meets": m}
                for fp, m in zip(prim["fips"], ehsa["also_meets"]) if m],
            "counts_by_metro": {c: {k: int(v) for k, v in row.items()} for c, row in cats_by_metro.iterrows()},
            "max_months_hot": int(mh.max()),
            "counties_hot_at_least_once": int((mh > 0).sum()),
            "counties_cold_at_least_once": int((prim["months_cold"] > 0).sum()),
            "share_bins_hot": r(ehsa["hot"].mean(), 3), "share_bins_cold": r(ehsa["cold"].mean(), 3),
            "longest_hot_counties": [{"county": names[fp], "months_hot": int(v)}
                                     for fp, v in mh.sort_values(ascending=False).head(15).items()],
            "ehsa_hot_categories_by_months_hot": top_hot.to_dict(orient="records"),
            "sensitivity": prim["sens"],
            "hot_months_per_county": {
                "0": int((mh == 0).sum()), "1-2": int(((mh >= 1) & (mh <= 2)).sum()),
                "3-5": int(((mh >= 3) & (mh <= 5)).sum()),
                "6-8": int(((mh >= 6) & (mh <= 8)).sum()), "9+": int((mh >= 9).sum())},
            "first_hot_by_wave": {k: int(v) for k, v in prim["first_wave"].value_counts().items()},
            "first_hot_by_month": {k: int(v) for k, v in prim["first_hot"].dropna().dt.strftime("%Y-%m")
                                   .value_counts().sort_index().items()},
            "hot_counties_per_month_max": int(prim["hot_pm"].max()),
            "hot_counties_per_month_max_month": f"{prim['hot_pm'].idxmax():%Y-%m}",
            "cold_counties_per_month_max": int(prim["cold_pm"].max()),
            "cold_counties_per_month_max_month": f"{prim['cold_pm'].idxmax():%Y-%m}",
            "final_month_hot": int(ehsa["hot"][-1].sum()), "final_month_cold": int(ehsa["cold"][-1].sum()),
            "hot_months_by_metro_mean": {g: r(mh[(m_keep == g).to_numpy()].mean(), 2)
                                         for g in ("Metro", "Nonmetro")},
        },
        "primary_vs_all_counties": {"primary": headline(prim, names), "all_counties": headline(allc, names)},
        "metro_vs_nonmetro_per_100k": {
            n: {f"{k}_{g.lower()}": r(ur.loc[n, (k, g)]) for k in ("cases", "deaths") for g in ("Metro", "Nonmetro")}
            for n in ur.index},
        "n_metro_counties": int((metro == "Metro").sum()),
    }

    # ---------------------------------------------------------- figures
    maps.set_style()
    g3083 = gdf.to_crs(data.STATIC_CRS)
    maps.wave_timeline(state_wk["cases"], state_wk["deaths"], waves, peaks, fig_dir / "wave_timeline.png")
    maps.rr_small_multiples(g3083, prim["rr_wave"].reindex(fips), prim["ref_wave"],
                            fig_dir / "case_rate_by_wave.png", excluded=ex_mask)
    maps.gi_small_multiples(g3083, prim["gi_wave"].reindex(fips), prim["gi_counts"],
                            fig_dir / "gi_star_by_wave.png", excluded=ex_mask)
    maps.ehsa_map(g3083, cats, fig_dir / "emerging_hot_spots.png", period_label,
                  metro_outline=g3083[(metro == "Metro").to_numpy()].dissolve(), excluded=ex_mask)
    maps.moran_chart(prim["moran"]["cases"], prim["moran"]["deaths"], waves, fig_dir / "morans_i_over_time.png")
    maps.urban_rural_chart(ur.loc[list(waves)], fig_dir / "metro_vs_nonmetro.png")
    maps.death_rate_map(g3083, prim["eb_total"]["deaths"].reindex(fips), fig_dir / "death_rate_cumulative.png",
                        period_label, excluded=ex_mask)
    maps.hot_spot_timing(g3083, months_hot.to_numpy(), first_wave, prim["hot_pm"], prim["cold_pm"], waves,
                         fig_dir / "hot_spot_timing.png", excluded=ex_mask)
    maps.hero(g3083, months_hot.to_numpy(), state_wk["cases"],
              {"omicron_peak": peaks["Omicron BA.1"]["cases_peak_per_100k"]}, fig_dir / "hero.png",
              excluded=ex_mask)

    # ---------------------------------------------------------- web map
    simp = gdf.to_crs(data.STATIC_CRS)
    try:
        simp["geometry"] = simp.geometry.simplify_coverage(800)  # shared edges stay shared
    except AttributeError:
        simp["geometry"] = simp.geometry.simplify(800)
    simp = simp.to_crs("EPSG:4326")
    months = [f"{t:%Y-%m}" for t in prim["rel"].index]
    eb_m = prim["eb_month"]["cases"].reindex(columns=fips)
    hot = pd.DataFrame(ehsa["hot"].astype(int) - ehsa["cold"].astype(int), columns=prim["fips"]).reindex(
        columns=fips).fillna(0).astype(int)

    def ints(s, fill=-1):
        return s.reindex(fips).round(0).fillna(fill).astype(int).to_numpy()

    props = pd.DataFrame({"n": gdf["name"].to_numpy(), "m": metro.to_numpy(),
                          "p": pop.to_numpy().astype(int), "x": ex_mask.astype(int),
                          "tc": final_cases.reindex(fips).astype(int).to_numpy(),
                          "c": [st.CATEGORIES.index(c) if c in st.CATEGORIES else -1 for c in cats],
                          "mh": ints(prim["months_hot"]), "mc": ints(prim["months_cold"]),
                          "cr": ints(prim["eb_total"]["cases"]), "dr": ints(prim["eb_total"]["deaths"])})
    for i in range(len(months)):
        props[f"r{i}"] = eb_m.iloc[i].round(0).fillna(-1).astype(int).to_numpy()
        props[f"h{i}"] = hot.iloc[i].to_numpy()
    props["fw"] = first_wave.to_numpy()
    for j, name in enumerate(waves):
        props[f"g{j}"] = prim["gi_wave"][name].reindex(fips).fillna(0).astype(int).to_numpy()
    ebd = prim["eb_total"]["deaths"]
    rate_vals = prim["eb_month"]["cases"].to_numpy()
    meta = {"months": months, "waves": list(waves), "categories": st.CATEGORIES,
            "colors": maps.EHSA_COLORS, "gi_colors": {str(k): v for k, v in maps.GI_COLORS.items()},
            "gi_labels": {str(k): v for k, v in maps.GI_LABELS.items()},
            "rate_breaks": maps.RATE_BREAKS, "rate_colors": maps.RATE_COLORS,
            "death_breaks": np.quantile(ebd, [.2, .4, .6, .8]).round(0).tolist(),
            "death_min": int(round(ebd.min())), "death_max": int(round(ebd.max())),
            "rate_min": int(np.floor(rate_vals.min())), "rate_max": int(np.ceil(rate_vals.max())),
            "period": period_label, "start_month": int(prim["hot_pm"].to_numpy().argmax()),
            "counts": metrics["emerging_hot_spots"]["counts"], "n_excluded": int(ex_mask.sum()),
            "excluded_label": maps.EXCL_LABEL}
    size = maps.web_map(simp, props, meta, doc_dir / "index.html", ROOT / "covid" / "web_template.html")
    metrics["web_map_bytes"] = int(size)
    metrics["runtime_seconds"] = r(time.time() - t0)
    (res_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(f"done in {metrics['runtime_seconds']} s; web map {size/1e3:.0f} kB")
    print(json.dumps(metrics["primary_vs_all_counties"], indent=1))


if __name__ == "__main__":
    main()
