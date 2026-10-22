"""COVID-19 space-time analysis of US state counties (default: Texas)."""
import os

# Some conda-forge MKL builds on Windows crash on the first matrix multiply when
# the Intel OpenMP runtime is missing. The sequential layer avoids that and is
# plenty fast for 254 counties.
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
