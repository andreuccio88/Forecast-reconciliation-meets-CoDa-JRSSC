######################################################################
## config.R
## Minimal configuration for the replication package.
######################################################################

PROJECT_ROOT <- getwd()
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
LEGACY_DIR   <- file.path(PROJECT_ROOT, "legacy")
OUT_DIR      <- file.path(PROJECT_ROOT, "outputs")
RDS_DIR      <- file.path(OUT_DIR, "rds")
TABLE_DIR    <- file.path(OUT_DIR, "tables")
FIGURE_DIR   <- file.path(OUT_DIR, "figures")

## Raw mortality data and life tables
DATA_RAW_FILE <- file.path(DATA_DIR, "Ita_Dx 80+.RData")
LT_DIR        <- file.path(PROJECT_ROOT, "LT")

## Output used by 01_point_forecasting.R
DATA_FILE <- file.path(DATA_DIR, "Output - baseline 1970-2000.RData")

## Forecast design
BASELINE_YEARS <- 1970:2000
FORECAST_H <- 15
FORECAST_END_YEAR <- 2015

## If TRUE, 00_generate_output_legacy.R writes the generated Output also to DATA_FILE.
## This lets 01_point_forecasting.R run immediately on the generated Output.
WRITE_GENERATED_OUTPUT_AS_DATA_FILE <- TRUE

REQUIRED_PACKAGES <- c(
  "dplyr", "tidyr", "stringr", "ggplot2", "readr", "rio",
  "FoReco", "patchwork", "tidyverse", "stats", "strucchange",
  "tables", "calibrate", "compositions", "abind", "MLmetrics",
  "FactoMineR", "tseries", "vars", "rTensor", "forecast", "r.jive"
)

for (dd in c(DATA_DIR, LT_DIR, OUT_DIR, RDS_DIR, TABLE_DIR, FIGURE_DIR)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}
