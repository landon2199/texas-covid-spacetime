"""Static figures (matplotlib) and the interactive web map (MapLibre GL JS)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

from .spacetime import CATEGORIES

INK = "#1f2933"
MUTED = "#6b7280"
EDGE = "#ffffff"
NOPAT = "#e9e9e6"
CASE_C = "#b4432f"
DEATH_C = "#3d4a5c"
WAVE_BG = "#f3efe6"
SOURCE = ("Data: The New York Times COVID-19 data (github.com/nytimes/covid-19-data); "
          "U.S. Census Bureau population estimates (Vintage 2023) and 2020 cartographic boundaries.")

# Esri-like emerging hot spot palette: warm = hot, cool = cold.
EHSA_COLORS = {
    "New Hot Spot": "#fcbba1", "Consecutive Hot Spot": "#fb6a4a",
    "Intensifying Hot Spot": "#99000d", "Persistent Hot Spot": "#cb181d",
    "Diminishing Hot Spot": "#fd8d3c", "Sporadic Hot Spot": "#fdd0a2",
    "Oscillating Hot Spot": "#c51b7d", "Historical Hot Spot": "#9e6b5c",
    "New Cold Spot": "#c6dbef", "Consecutive Cold Spot": "#6baed6",
    "Intensifying Cold Spot": "#08306b", "Persistent Cold Spot": "#2171b5",
    "Diminishing Cold Spot": "#41b6c4", "Sporadic Cold Spot": "#deebf7",
    "Oscillating Cold Spot": "#807dba", "Historical Cold Spot": "#6c8299",
    "No Pattern Detected": NOPAT,
}
GI_COLORS = {-3: "#2166ac", -2: "#67a9cf", -1: "#d1e5f0", 0: NOPAT,
             1: "#fddbc7", 2: "#ef8a62", 3: "#b2182b"}
GI_LABELS = {3: "Hot spot, q < 0.01", 2: "Hot spot, q < 0.05", 1: "Hot spot, q < 0.10",
             0: "Not significant", -1: "Cold spot, q < 0.10", -2: "Cold spot, q < 0.05",
             -3: "Cold spot, q < 0.01"}
RR_BREAKS = [0, 0.5, 0.67, 0.8, 0.91, 1.1, 1.25, 1.5, 2.0, 99]
RR_COLORS = ["#2c7bb6", "#6baed6", "#abd9e9", "#e0f0f3", "#f7f7f2",
             "#fee8c8", "#fdbb84", "#e34a33", "#a50f15"]
RATE_BREAKS = [0, 50, 100, 250, 500, 1000, 2500, 5000]
RATE_COLORS = ["#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c",
               "#e6550d", "#a63603", "#6b2106"]


def set_style():
    mpl.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans"],
        "font.size": 10, "text.color": INK, "axes.labelcolor": INK,
        "axes.edgecolor": "#c9c9c4", "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "axes.titlesize": 11, "axes.titlelocation": "left",
        "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def _title(fig, title, subtitle=None, y=0.98):
    fig.text(0.02, y, title, fontsize=15, fontweight="bold", ha="left", va="top")
    if subtitle:
        fig.text(0.02, y - 0.045, subtitle, fontsize=10.5, color=MUTED, ha="left", va="top")


def _source(fig, extra=""):
    fig.text(0.02, 0.005, (extra + " " if extra else "") + SOURCE, fontsize=7, color=MUTED,
             ha="left", va="bottom", wrap=True)


def _map_ax(ax):
    ax.set_axis_off()
    ax.set_aspect("equal")


EXCL_FC, EXCL_EC = "#d9d9d5", "#7d7d78"
EXCL_LABEL = "Excluded: reported cases exceed resident population"


def _hatch(ax, gdf, excluded, fc=EXCL_FC, ec=EXCL_EC, label=True):
    """Draw excluded counties hatched grey; return a legend handle list."""
    if excluded is None or not np.any(excluded):
        return []
    ex = np.asarray(excluded, bool)
    gdf[ex].plot(ax=ax, facecolor=fc, edgecolor=ec, hatch="////", linewidth=0.4)
    return [Patch(fc=fc, ec=ec, hatch="////", lw=0.4, label=f"{EXCL_LABEL} ({int(ex.sum())})")] if label else []


def _shade_waves(ax, waves, label=True, ymax_frac=0.97):
    for name, (a, b) in waves.items():
        ax.axvspan(pd.Timestamp(a), pd.Timestamp(b), color=WAVE_BG, zorder=0, lw=0)
        if label:
            ax.text(pd.Timestamp(a) + (pd.Timestamp(b) - pd.Timestamp(a)) / 2, ymax_frac, name,
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=8, color=MUTED)


# ------------------------------------------------------------------ figures
def wave_timeline(cases_wk: pd.Series, deaths_wk: pd.Series, waves, peaks, out: Path):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.12})
    for ax in (a1, a2):
        _shade_waves(ax, waves, label=(ax is a1))
    a1.fill_between(cases_wk.index, cases_wk.values, color=CASE_C, alpha=0.18, lw=0, step="pre")
    a1.plot(cases_wk.index, cases_wk.values, color=CASE_C, lw=1.4, drawstyle="steps-pre")
    a2.fill_between(deaths_wk.index, deaths_wk.values, color=DEATH_C, alpha=0.18, lw=0, step="pre")
    a2.plot(deaths_wk.index, deaths_wk.values, color=DEATH_C, lw=1.4, drawstyle="steps-pre")
    for name, p in peaks.items():
        t = pd.Timestamp(p["cases_peak_week"])
        a1.annotate(f"{p['cases_peak_per_100k']:,.0f}", (t, p["cases_peak_per_100k"]),
                    xytext=(0, 4), textcoords="offset points", ha="center", fontsize=8.5,
                    color=CASE_C, fontweight="bold")
        t = pd.Timestamp(p["deaths_peak_week"])
        a2.annotate(f"{p['deaths_peak_per_100k']:.1f}", (t, p["deaths_peak_per_100k"]),
                    xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8.5,
                    color=DEATH_C, fontweight="bold")
    a1.set_ylabel("New cases per 100,000\n(per week)")
    a2.set_ylabel("Deaths per 100,000\n(per week)")
    a1.set_ylim(0, cases_wk.max() * 1.12)
    a2.set_ylim(0, deaths_wk.max() * 1.25)
    a2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    a1.grid(axis="y", color="#ececec", lw=0.6)
    a2.grid(axis="y", color="#ececec", lw=0.6)
    _title(fig, "Five waves of reported COVID-19 in Texas",
           "Statewide reported cases and deaths per 100,000 residents, by week ending Sunday; "
           "peak week labelled", y=1.03)
    _source(fig, "Weeks end on Sunday. One-week death spikes (late July 2020, mid-November 2022) are reporting "
                 "artefacts: Texas changed how it counted deaths in July 2020, and the 2022 spike looks like a backlog.")
    fig.subplots_adjust(top=0.9, bottom=0.12)
    fig.savefig(out)
    plt.close(fig)


def rr_small_multiples(gdf, rr: pd.DataFrame, state_rates: dict, out: Path, excluded=None):
    names = list(rr.columns)
    fig, axes = plt.subplots(1, len(names), figsize=(16, 4.6))
    cmap = ListedColormap(RR_COLORS)
    norm = BoundaryNorm(RR_BREAKS, cmap.N)
    for ax, name in zip(axes, names):
        g = gdf.assign(v=rr[name].to_numpy())
        g.plot(column="v", cmap=cmap, norm=norm, ax=ax, edgecolor=EDGE, linewidth=0.2)
        ex_handles = _hatch(ax, gdf, excluded)
        _map_ax(ax)
        ax.set_title(f"{name}\nReference: {state_rates[name]:,.0f} per 100k", fontsize=10.5)
    labels = ["< 0.5", "0.5-0.67", "0.67-0.8", "0.8-0.91", "0.91-1.1",
              "1.1-1.25", "1.25-1.5", "1.5-2", "> 2"]
    handles = [Patch(fc=c, ec="#bbbbbb", lw=0.3, label=l) for c, l in zip(RR_COLORS, labels)]
    handles += ex_handles
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9,
               title="County case rate relative to the pooled rate of the analysed counties in the same wave "
                     "(Empirical Bayes smoothed)", title_fontsize=9.5, bbox_to_anchor=(0.5, 0.04))
    _title(fig, "Where cases ran above and below the state rate, wave by wave",
           "Ratio of each county's smoothed case rate to the pooled rate of the analysed counties; blue = lower, red = higher",
           y=1.05)
    _source(fig)
    fig.subplots_adjust(top=0.82, bottom=0.2, wspace=0.02)
    fig.savefig(out)
    plt.close(fig)


def gi_small_multiples(gdf, classes: pd.DataFrame, counts: dict, out: Path, measure="case", excluded=None):
    names = list(classes.columns)
    fig, axes = plt.subplots(1, len(names), figsize=(16, 4.6))
    for ax, name in zip(axes, names):
        cols = classes[name].map(GI_COLORS).fillna(NOPAT).to_numpy()
        gdf.plot(color=cols, ax=ax, edgecolor=EDGE, linewidth=0.2)
        ex_handles = _hatch(ax, gdf, excluded)
        _map_ax(ax)
        c = counts[name]
        ax.set_title(f"{name}\n{c['hot']} hot, {c['cold']} cold counties", fontsize=10.5)
    handles = [Patch(fc=GI_COLORS[k], ec="#bbbbbb", lw=0.3, label=GI_LABELS[k])
               for k in (3, 2, 1, 0, -1, -2, -3)] + ex_handles
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9,
               title="Getis-Ord Gi* with Benjamini-Hochberg false discovery rate (q = FDR-adjusted p)",
               title_fontsize=9.5, bbox_to_anchor=(0.5, 0.04))
    _title(fig, f"Hot and cold spots of the {measure} rate, wave by wave",
           "Gi* on each wave's Empirical Bayes smoothed rate; neighbourhood = county + queen-contiguous neighbours",
           y=1.05)
    _source(fig)
    fig.subplots_adjust(top=0.82, bottom=0.2, wspace=0.02)
    fig.savefig(out)
    plt.close(fig)


def _ehsa_legend(cats: pd.Series):
    counts = cats.value_counts()
    return [Patch(fc=EHSA_COLORS[c], ec="#bbbbbb", lw=0.3, label=f"{c}  ({counts[c]})")
            for c in CATEGORIES if c in counts.index]


def ehsa_map(gdf, cats: pd.Series, out: Path, period: str, metro_outline=None, excluded=None):
    fig, ax = plt.subplots(figsize=(10, 8.6))
    gdf.plot(color=cats.map(EHSA_COLORS).fillna(NOPAT).to_numpy(), ax=ax, edgecolor=EDGE, linewidth=0.3)
    ex_handles = _hatch(ax, gdf, excluded)
    if metro_outline is not None:
        metro_outline.boundary.plot(ax=ax, color=INK, linewidth=0.6, alpha=0.55)
    _map_ax(ax)
    leg = ax.legend(handles=_ehsa_legend(cats) + ex_handles, loc="lower left", frameon=False, fontsize=9.5,
                    title="Category (number of counties)", title_fontsize=10,
                    bbox_to_anchor=(-0.02, 0.0))
    leg._legend_box.align = "left"
    if metro_outline is not None:
        ax.plot([], [], color=INK, lw=0.8, alpha=0.55, label="Metro county outline")
    _title(fig, "Emerging hot spots of COVID-19 cases in Texas counties",
           f"Space-time Gi* on monthly case rate relative to the analysed counties, {period}. Categories describe the "
           "end of the series,\n"
           "when few counties stood out; see the timing figure for where hot spots appeared earlier.")
    _source(fig, "Dark outline: metro counties (USDA ERS Rural-Urban Continuum Codes 2023, codes 1-3).")
    fig.subplots_adjust(top=0.9, bottom=0.04)
    fig.savefig(out)
    plt.close(fig)


def moran_chart(mc: pd.DataFrame, md: pd.DataFrame, waves, out: Path):
    fig, ax = plt.subplots(figsize=(11, 4.6))
    _shade_waves(ax, waves)
    for df, color, label in ((mc, CASE_C, "Cases per 100,000"), (md, DEATH_C, "Deaths per 100,000")):
        x = pd.to_datetime(df.index) - pd.offsets.MonthBegin(1) + pd.Timedelta(days=14)
        ax.plot(x, df["I"], color=color, lw=1.8, label=label, marker="o", ms=0)
        sig = df["p_sim"] < 0.05
        ax.scatter(x[sig], df["I"][sig], color=color, s=22, zorder=3)
        ax.scatter(x[~sig], df["I"][~sig], facecolor="white", edgecolor=color, s=22, zorder=3)
    ax.axhline(0, color="#999999", lw=0.8)
    ax.set_ylabel("Global Moran's I")
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.grid(axis="y", color="#ececec", lw=0.6)
    ax.scatter([], [], color=MUTED, s=22, label="p < 0.05 (999 permutations)")
    ax.scatter([], [], facecolor="white", edgecolor=MUTED, s=22, label="not significant")
    ax.legend(frameon=False, loc="upper right", ncol=2, fontsize=9, bbox_to_anchor=(1, 0.9))
    ax.set_ylim(min(-0.05, mc["I"].min() - 0.05, md["I"].min() - 0.05), max(mc["I"].max(), md["I"].max()) * 1.35)
    _title(fig, "How spatially clustered was COVID-19 in Texas, month by month?",
           "Global Moran's I of monthly Empirical Bayes smoothed rates, queen contiguity "
           "(row-standardised). Higher = similar rates sit next to each other.", y=1.06)
    _source(fig)
    fig.subplots_adjust(top=0.84, bottom=0.17)
    fig.savefig(out)
    plt.close(fig)


def urban_rural_chart(tab: pd.DataFrame, out: Path):
    """tab: index = wave names, columns MultiIndex (measure, Metro/Nonmetro)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    colors = {"Metro": "#5b6b7f", "Nonmetro": "#c98b3a"}
    for ax, (measure, title) in zip(axes, (("cases", "Reported cases per 100,000"),
                                           ("deaths", "Deaths per 100,000"))):
        sub = tab[measure]
        x = np.arange(len(sub))
        for k, grp in enumerate(("Metro", "Nonmetro")):
            v = sub[grp].to_numpy()
            ax.bar(x + (k - 0.5) * 0.38, v, width=0.36, color=colors[grp], label=grp)
            for xi, vi in zip(x, v):
                ax.text(xi + (k - 0.5) * 0.38, vi, f"{vi:,.0f}", ha="center", va="bottom", fontsize=7.5,
                        color=MUTED)
        ax.set_xticks(x, [n.replace(" and later", "+").replace(" ", "\n", 1) for n in sub.index], fontsize=8.5)
        ax.set_title(title)
        ax.grid(axis="y", color="#ececec", lw=0.6)
    axes[0].legend(frameon=False)
    _title(fig, "Metro and nonmetro Texas, wave by wave",
           "Pooled rates (all counties in each group); metro = USDA ERS RUCC 2023 codes 1-3. "
           "Rates are not age-adjusted.", y=1.08)
    _source(fig)
    fig.subplots_adjust(top=0.82, bottom=0.2, wspace=0.18)
    fig.savefig(out)
    plt.close(fig)


def death_rate_map(gdf, rate: pd.Series, out: Path, period: str, excluded=None):
    fig, ax = plt.subplots(figsize=(10, 8.6))
    q = np.nanquantile(rate, [0, .2, .4, .6, .8, 1])
    q[-1] += 1e-6
    cols = ["#f1eef6", "#bdc9e1", "#74a9cf", "#2b8cbe", "#045a8d"]
    cmap = ListedColormap(cols)
    gdf.assign(v=rate.to_numpy()).plot(column="v", cmap=cmap, norm=BoundaryNorm(q, cmap.N), ax=ax,
                                        edgecolor=EDGE, linewidth=0.3)
    ex_handles = _hatch(ax, gdf, excluded)
    _map_ax(ax)
    handles = [Patch(fc=c, ec="#bbbbbb", lw=0.3, label=f"{q[i]:,.0f} - {q[i+1]:,.0f}")
               for i, c in enumerate(cols)] + ex_handles
    ax.legend(handles=handles, loc="lower left", frameon=False, title="Deaths per 100,000\n(quintiles)",
              fontsize=9.5, title_fontsize=10)
    _title(fig, "COVID-19 deaths per 100,000 residents, by county",
           f"Cumulative, {period}; Empirical Bayes smoothed; not age-adjusted")
    _source(fig)
    fig.subplots_adjust(top=0.9, bottom=0.04)
    fig.savefig(out)
    plt.close(fig)


MONTHS_HOT_BINS = [(0, 0, "Never"), (1, 2, "1-2"), (3, 5, "3-5"), (6, 8, "6-8"), (9, 99, "9 or more")]
MONTHS_HOT_COLORS = [NOPAT, "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"]
MONTHS_HOT_DARK = ["#2b313b", "#fdd49e", "#fc8d59", "#e34a33", "#b30000"]
WAVE_COLORS = {"Summer 2020": "#7f0000", "Winter 2020-21": "#d7301f", "Between waves": "#9e9ac8",
               "Delta": "#fc8d59", "Omicron BA.1": "#fdbb84", "Omicron BA.5 and later": "#fee8c8",
               "Never a hot spot": NOPAT}


def months_hot_class(m):
    m = np.asarray(m)
    out = np.zeros(m.shape, int)
    for k, (lo, hi, _) in enumerate(MONTHS_HOT_BINS):
        out[(m >= lo) & (m <= hi)] = k
    return out


def hot_spot_timing(gdf, months_hot, first_wave: pd.Series, hot_pm: pd.Series, cold_pm: pd.Series,
                    waves, out: Path, excluded=None):
    fig = plt.figure(figsize=(15, 9.4))
    a1 = fig.add_axes([0.0, 0.33, 0.5, 0.56])
    a2 = fig.add_axes([0.5, 0.33, 0.5, 0.56])
    a3 = fig.add_axes([0.06, 0.08, 0.9, 0.17])
    keep = np.ones(len(gdf), bool) if excluded is None else ~np.asarray(excluded, bool)
    k = months_hot_class(months_hot)
    gdf.plot(color=np.array(MONTHS_HOT_COLORS)[k], ax=a1, edgecolor=EDGE, linewidth=0.25)
    ex_handles = _hatch(a1, gdf, excluded)
    _map_ax(a1)
    a1.set_title("Months as a significant hot spot (of 36)", fontsize=12)
    cnt = np.bincount(k[keep], minlength=len(MONTHS_HOT_BINS))
    a1.legend(handles=[Patch(fc=c, ec="#bbbbbb", lw=0.3, label=f"{lab}  ({n})")
                       for c, (_, _, lab), n in zip(MONTHS_HOT_COLORS, MONTHS_HOT_BINS, cnt)] + ex_handles,
              loc="lower left", frameon=False, fontsize=9.5, title="Months (counties)", title_fontsize=10)
    gdf.plot(color=first_wave.map(WAVE_COLORS).fillna(NOPAT).to_numpy(), ax=a2, edgecolor=EDGE, linewidth=0.25)
    _hatch(a2, gdf, excluded, label=False)
    _map_ax(a2)
    a2.set_title("When each county first became a hot spot", fontsize=12)
    vc = first_wave.value_counts()
    order = ["Summer 2020", "Winter 2020-21", "Between waves", "Delta", "Omicron BA.1",
             "Omicron BA.5 and later", "Never a hot spot"]
    a2.legend(handles=[Patch(fc=WAVE_COLORS[w], ec="#bbbbbb", lw=0.3, label=f"{w}  ({vc[w]})")
                       for w in order if w in vc.index],
              loc="lower left", frameon=False, fontsize=9.5, title="First hot month falls in (counties)",
              title_fontsize=10)
    _shade_waves(a3, waves, ymax_frac=0.97)
    x = hot_pm.index - pd.offsets.MonthBegin(1) + pd.Timedelta(days=14)
    a3.bar(x, hot_pm.values, width=24, color=CASE_C, label="hot spot counties")
    a3.bar(x, -cold_pm.values, width=24, color="#2171b5", label="cold spot counties")
    a3.axhline(0, color="#999999", lw=0.6)
    a3.set_ylabel("Counties")
    a3.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
    a3.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    a3.legend(frameon=False, fontsize=9, loc="lower right", ncol=2, bbox_to_anchor=(1, 1.0))
    a3.set_ylim(-cold_pm.max() * 1.3 - 2, hot_pm.max() * 1.35)
    a3.set_title("Significant hot (up) and cold (down) spot counties each month", fontsize=11, pad=6)
    _title(fig, "Where and when hot spots emerged",
           "Space-time Getis-Ord Gi* of each county's monthly case rate relative to the analysed counties (county + neighbours, "
           "this month and last); Benjamini-Hochberg FDR, q < 0.05", y=0.995)
    _source(fig)
    fig.savefig(out)
    plt.close(fig)


def hero(gdf, months_hot, cases_wk: pd.Series, headline: dict, out: Path, excluded=None):
    bg, fg, sub = "#14181f", "#f5f3ee", "#a7adb7"
    fig = plt.figure(figsize=(16, 9), dpi=100, facecolor=bg)
    axm = fig.add_axes([0.40, 0.03, 0.58, 0.94], facecolor=bg)
    k = months_hot_class(months_hot)
    gdf.plot(color=np.array(MONTHS_HOT_DARK)[k], ax=axm, edgecolor=bg, linewidth=0.35)
    n_ex = _hatch(axm, gdf, excluded, fc="#2b313b", ec="#8b93a1")
    keep = np.ones(len(gdf), bool) if excluded is None else ~np.asarray(excluded, bool)
    _map_ax(axm)
    fig.text(0.04, 0.88, "COVID-19 in Texas", color=fg, fontsize=40, fontweight="bold", va="top")
    fig.text(0.04, 0.79, "Where and when did hot spots emerge?", color=fg, fontsize=22, va="top")
    fig.text(0.04, 0.73, "254 counties  ·  March 2020 to February 2023  ·  space-time hot spot analysis",
             color=sub, fontsize=13, va="top")
    axt = fig.add_axes([0.045, 0.36, 0.33, 0.25], facecolor=bg)
    axt.fill_between(cases_wk.index, cases_wk.values, color="#e4572e", alpha=0.35, lw=0)
    axt.plot(cases_wk.index, cases_wk.values, color="#ff7b54", lw=1.4)
    for s in axt.spines.values():
        s.set_visible(False)
    axt.tick_params(colors=sub, labelsize=10, length=0)
    axt.set_yticks([])
    axt.xaxis.set_major_locator(mdates.YearLocator())
    axt.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axt.text(0, 1.04, "Weekly reported cases per 100,000, statewide", transform=axt.transAxes, color=sub,
             fontsize=11)
    axt.text(0.66, 0.8, f"Omicron peak, Jan 2022:\n{headline['omicron_peak']:,.0f} per 100,000",
             transform=axt.transAxes, color=fg, fontsize=11, va="top")
    fig.text(0.045, 0.285, "Months each county was a significant hot spot", color=sub, fontsize=11.5)
    cnt = np.bincount(k[keep], minlength=len(MONTHS_HOT_BINS))
    x = 0.045
    for c, (_, _, lab), n in zip(MONTHS_HOT_DARK, MONTHS_HOT_BINS, cnt):
        fig.patches.append(mpl.patches.Rectangle((x, 0.225), 0.016, 0.03, transform=fig.transFigure, color=c))
        fig.text(x + 0.021, 0.24, lab, color=fg, fontsize=11.5, va="center")
        x += 0.067
    if n_ex:
        fig.patches.append(mpl.patches.Rectangle((0.045, 0.165), 0.016, 0.03, transform=fig.transFigure,
                                                  facecolor="#2b313b", edgecolor="#8b93a1", hatch="////"))
        fig.text(0.066, 0.18, "Excluded: reported cases exceed resident population", color=sub,
                 fontsize=10.5, va="center")
    fig.text(0.04, 0.03, "Getis-Ord Gi* space-time cube, FDR-controlled  ·  Empirical Bayes rates  ·  "
             "data: The New York Times, U.S. Census Bureau", color=sub, fontsize=10)
    with mpl.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(out, dpi=100, facecolor=bg)
    plt.close(fig)


# ------------------------------------------------------------------ web map
def web_map(gdf_wgs, props: pd.DataFrame, meta: dict, out: Path, template: Path):
    """Write a self-contained MapLibre page. gdf_wgs: simplified EPSG:4326
    geometries; props: one row per county, columns become feature properties."""
    feats = []
    for geom, (_, row) in zip(gdf_wgs.geometry, props.iterrows()):
        gj = json.loads(json.dumps(geom.__geo_interface__))
        gj["coordinates"] = _round_coords(gj["coordinates"])
        p = {k: (v.item() if hasattr(v, "item") else v) for k, v in row.items()}
        feats.append({"type": "Feature", "geometry": gj, "properties": p})
    fc = {"type": "FeatureCollection", "features": feats}
    html = template.read_text(encoding="utf-8")
    html = html.replace("/*__DATA__*/null", json.dumps(fc, separators=(",", ":")))
    html = html.replace("/*__META__*/null", json.dumps(meta, separators=(",", ":")))
    out.write_text(html, encoding="utf-8")
    return out.stat().st_size


def _round_coords(c, nd=4):
    if isinstance(c[0], (int, float)):
        return [round(c[0], nd), round(c[1], nd)]
    return [_round_coords(x, nd) for x in c]
