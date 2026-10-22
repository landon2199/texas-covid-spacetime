"""Download and load the public inputs: NYT county case/death counts, Census
population estimates, Census county boundaries and USDA Rural-Urban Continuum
Codes. Also turns cumulative counts into new counts per period."""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

NYT_URL = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties-{year}.csv"
POP_URL = ("https://www2.census.gov/programs-surveys/popest/datasets/2020-2023/"
           "counties/totals/co-est2023-alldata.csv")
SHP_URL = "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_500k.zip"
RUCC_URL = ("https://ers.usda.gov/sites/default/files/_laserfiche/DataFiles/53251/"
            "Ruralurbancontinuumcodes2023.csv")
NYT_YEARS = (2020, 2021, 2022, 2023)

# Texas Centric Albers Equal Area; used for every static map.
STATIC_CRS = "EPSG:3083"


def _fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
        tmp.replace(dest)
    return dest


def download_all(raw: Path) -> None:
    for y in NYT_YEARS:
        _fetch(NYT_URL.format(year=y), raw / f"us-counties-{y}.csv")
    _fetch(POP_URL, raw / "co-est2023-alldata.csv")
    _fetch(SHP_URL, raw / "cb_2020_us_county_500k.zip")
    _fetch(RUCC_URL, raw / "ruralurbancodes2023.csv")


def load_nyt(raw: Path, state_fips: str) -> tuple[pd.DataFrame, dict]:
    """Daily cumulative cases/deaths for one state's counties.

    Rows with no FIPS (NYT's "Unknown" county) cannot be mapped; their final
    cumulative totals are returned in the info dict so they can be reported.
    """
    frames = [pd.read_csv(raw / f"us-counties-{y}.csv", dtype={"fips": str},
                          parse_dates=["date"]) for y in NYT_YEARS]
    df = pd.concat(frames, ignore_index=True)
    df = df[df["fips"].str.startswith(state_fips, na=False) |
            (df["fips"].isna() & df["state"].eq(_state_name(df, state_fips)))]
    # NYT "Unknown" rows have no FIPS; "Pending County Assignment" uses xx999
    df.loc[df["fips"].str.endswith("999", na=False), "fips"] = np.nan
    unknown = df[df["fips"].isna()]
    last_unknown = unknown.sort_values("date").groupby("county").tail(1)
    info = {
        "unknown_county_cases_final": int(last_unknown["cases"].sum()),
        "unknown_county_deaths_final": int(last_unknown["deaths"].fillna(0).sum()),
        "first_date": str(df["date"].min().date()),
        "last_date": str(df["date"].max().date()),
    }
    df = df[df["fips"].notna()].copy()
    df["deaths"] = df["deaths"].fillna(0)
    return df[["date", "fips", "cases", "deaths"]], info


def _state_name(df: pd.DataFrame, state_fips: str) -> str:
    m = df.loc[df["fips"].str.startswith(state_fips, na=False), "state"]
    return m.mode().iloc[0]


def daily_cumulative(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Wide date x county table of cumulative counts on a complete daily grid.
    Before a county's first report the count is 0; gaps are forward-filled."""
    wide = df.pivot_table(index="date", columns="fips", values=col, aggfunc="sum")
    full = pd.date_range(wide.index.min(), wide.index.max(), freq="D")
    return wide.reindex(full).ffill().fillna(0.0)


def cumulative_to_new(cum: pd.DataFrame | pd.Series) -> tuple:
    """New counts per step from cumulative counts.

    NYT cumulative series sometimes go *down* when a health department removes
    duplicates or moves cases to another county. A drop means the earlier
    cumulative values were too high, so each value is replaced by the lowest
    value reported on or after that date (a backward running minimum) before
    differencing. New counts are then never negative and they add up exactly
    to the final cumulative total (simply clipping negative differences at 0
    would count the removed cases twice).
    Returns (new_counts, number_of_downward_revisions).
    """
    raw_diff = cum.diff()
    neg = raw_diff < 0
    n_neg = int(neg.to_numpy().sum())
    adj = cum.iloc[::-1].cummin().iloc[::-1]
    new = adj.diff()
    new.iloc[0] = adj.iloc[0]
    return new, n_neg


def to_period_new(cum_daily: pd.DataFrame, freq: str) -> tuple[pd.DataFrame, int]:
    """Sample cumulative counts at the end of each period (e.g. 'W-SUN' weeks,
    'ME' months) and difference them. Only complete periods are kept.
    Downward revisions are corrected on the daily series first."""
    last = cum_daily.index.max()
    cum_daily = cum_daily.iloc[::-1].cummin().iloc[::-1]
    sampled = cum_daily.resample(freq).last()
    if sampled.index[-1] > last:  # incomplete final period
        sampled = sampled.iloc[:-1]
    return cumulative_to_new(sampled)


def load_population(raw: Path, state_fips: str, year: int = 2021) -> pd.Series:
    pop = pd.read_csv(raw / "co-est2023-alldata.csv", encoding="latin-1",
                      dtype={"STATE": str, "COUNTY": str})
    pop = pop[(pop["SUMLEV"] == 50) & (pop["STATE"] == state_fips)]
    s = pop.set_index(pop["STATE"] + pop["COUNTY"])[f"POPESTIMATE{year}"]
    s.index.name = "fips"
    return s.astype(float)


def load_counties(raw: Path, state_fips: str) -> gpd.GeoDataFrame:
    zpath = raw / "cb_2020_us_county_500k.zip"
    gdf = gpd.read_file(f"zip://{zpath.as_posix()}")
    gdf = gdf[gdf["STATEFP"] == state_fips].copy()
    gdf["fips"] = gdf["GEOID"]
    return gdf[["fips", "NAME", "geometry"]].rename(columns={"NAME": "name"}) \
        .sort_values("fips").reset_index(drop=True)


def load_rucc(raw: Path, state_fips: str) -> pd.Series:
    """USDA ERS 2023 Rural-Urban Continuum Code per county (1-3 metro, 4-9 nonmetro)."""
    r = pd.read_csv(raw / "ruralurbancodes2023.csv", dtype={"FIPS": str}, encoding="latin-1")
    r = r[(r["Attribute"] == "RUCC_2023") & r["FIPS"].str.startswith(state_fips)]
    return r.set_index("FIPS")["Value"].astype(int).rename("rucc")


def metro_label(rucc: pd.Series) -> pd.Series:
    return pd.Series(np.where(rucc <= 3, "Metro", "Nonmetro"), index=rucc.index)
