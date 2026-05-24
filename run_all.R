######################################################################
## run_all.R
## Replication script for JRSSC revision package
##
## Minimal-output version:
## after the pipeline finishes, outputs/ contains only the CSV and TEX
## files needed to reproduce the manuscript tables/results.
######################################################################

rm(list = ls())

## -------------------------------------------------------------------
## Project paths
## -------------------------------------------------------------------

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
CODE_DIR <- file.path(PROJECT_ROOT, "R")
DATA_DIR <- file.path(PROJECT_ROOT, "data")
LT_DIR   <- file.path(PROJECT_ROOT, "LT")
OUT_DIR  <- file.path(PROJECT_ROOT, "outputs")

## Minimal output controls

CREATE_FIGURES <- TRUE
PRUNE_TO_MANUSCRIPT_TABLES <- FALSE

MANUSCRIPT_TABLE_FILES <- c(
  ## Point forecast accuracy: Table 2 and Appendix Table 5
  "RMSE_all_structures.csv",
  "MAE_all_structures.csv",
  "point_RMSE_table.tex",
  "point_MAE_table.tex",

  ## Model-based CoDA bootstrap CRPS: Table 3
  "CRPS_all_structures.csv",
  "interval_CRPS_table.tex",

  ## Residual block-bootstrap robustness CRPS: Table 4
  "BB_CRPS_all_structures.csv",
  "BB_interval_CRPS_table.tex"
)

cat("\n============================================================\n")
cat("JRSSC revision replication pipeline\n")
cat("Project root:", PROJECT_ROOT, "\n")
cat("Code folder: ", CODE_DIR, "\n")
cat("Data folder: ", DATA_DIR, "\n")
cat("Minimal outputs:", PRUNE_TO_MANUSCRIPT_TABLES, "\n")
cat("============================================================\n")

## -------------------------------------------------------------------
## Check folder structure
## -------------------------------------------------------------------

if (!dir.exists(CODE_DIR)) {
  stop(
    "Cannot find the R/ folder.\n",
    "Please run this script from the project root, not from inside R/.\n\n",
    "Current working directory is:\n",
    PROJECT_ROOT
  )
}

if (!dir.exists(DATA_DIR)) {
  stop(
    "Cannot find the data/ folder.\n",
    "Please create data/ in the project root and place the input .RData files there.\n\n",
    "Current working directory is:\n",
    PROJECT_ROOT
  )
}

if (!dir.exists(LT_DIR)) {
  stop(
    "Cannot find the LT/ folder.\n",
    "Please create LT/ in the project root and place the HMD life-table files there.\n\n",
    "Current working directory is:\n",
    PROJECT_ROOT
  )
}

## -------------------------------------------------------------------
## Clean and recreate output folders
## -------------------------------------------------------------------

if (dir.exists(OUT_DIR)) {
  cat("\nRemoving previous outputs folder...\n")
  unlink(OUT_DIR, recursive = TRUE)
}

dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "rds"), recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------
## Check required scripts
## -------------------------------------------------------------------

required_scripts <- c(
  "base_fcast.R",
  "reco_fcast.R",
  "gen_boot.R",
  "reco_boot.R",
  "Base function.R",
  "CoDa functions.R",
  "coda_boot.R",
  "config.R",
  "helpers.R",
  "residual_block_bootstrap_crps_legacy_fixed.R"
)

script_paths <- file.path(CODE_DIR, required_scripts)
missing_scripts <- required_scripts[!file.exists(script_paths)]

if (length(missing_scripts) > 0) {
  stop(
    "The following required scripts are missing in R/:\n",
    paste(missing_scripts, collapse = "\n"),
    "\n\nR/ folder is:\n",
    CODE_DIR
  )
}

cat("\nAll required scripts found.\n")

## -------------------------------------------------------------------
## Check required data and life-table files
## -------------------------------------------------------------------

required_data <- c("Ita_Dx 80+.RData")
required_lt <- c("mltper_5x1.txt", "fltper_5x1.txt", "bltper_5x1.txt")

lt_paths <- file.path(LT_DIR, required_lt)
missing_lt <- required_lt[!file.exists(lt_paths)]

if (length(missing_lt) > 0) {
  stop(
    "The following required life-table files are missing in LT/:\n",
    paste(missing_lt, collapse = "\n"),
    "\n\nLT folder is:\n",
    LT_DIR
  )
}
cat("All required life-table files found.\n")

data_paths <- file.path(DATA_DIR, required_data)
missing_data <- required_data[!file.exists(data_paths)]

if (length(missing_data) > 0) {
  stop(
    "The following required data files are missing in data/:\n",
    paste(missing_data, collapse = "\n"),
    "\n\nData folder is:\n",
    DATA_DIR
  )
}
cat("All required data files found.\n")

## -------------------------------------------------------------------
## Check required packages
## -------------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "dplyr",
  "tidyr",
  "ggplot2",
  "forecast",
  "FoReco",
  "rio",
  "compositions"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running the pipeline:\n",
    paste(missing_packages, collapse = ", "),
    "\n\nExample:\ninstall.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))"
  )
}

cat("All required packages are installed.\n")

## -------------------------------------------------------------------
## Helper to run each script safely
## -------------------------------------------------------------------

run_script <- function(script_name, working_dir) {
  script_path <- file.path(CODE_DIR, script_name)

  if (!file.exists(script_path)) {
    stop("Cannot find script: ", script_path)
  }

  cat("\n------------------------------------------------------------\n")
  cat("Running:", script_path, "\n")
  cat("Working directory:", working_dir, "\n")
  cat("------------------------------------------------------------\n")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(working_dir)

  script_env <- new.env(parent = globalenv())
  script_env$PROJECT_ROOT <- PROJECT_ROOT
  script_env$CODE_DIR <- CODE_DIR
  script_env$DATA_DIR <- DATA_DIR
  script_env$LT_DIR <- LT_DIR
  script_env$OUT_DIR <- OUT_DIR
  script_env$CREATE_FIGURES <- CREATE_FIGURES
  script_env$PRUNE_TO_MANUSCRIPT_TABLES <- PRUNE_TO_MANUSCRIPT_TABLES

  source(script_path, local = script_env)

  invisible(TRUE)
}

prune_outputs <- function() {
  table_dir <- file.path(OUT_DIR, "tables")

  if (dir.exists(table_dir)) {
    existing <- list.files(table_dir, full.names = FALSE, recursive = FALSE)
    remove <- setdiff(existing, MANUSCRIPT_TABLE_FILES)

    if (length(remove) > 0) {
      unlink(file.path(table_dir, remove), recursive = TRUE, force = TRUE)
    }
  }

  ## Remove non-table artifacts from final outputs. They are caches or plots,
  ## not required to reproduce the manuscript tables.
  top_level <- list.files(OUT_DIR, full.names = TRUE, recursive = FALSE, all.files = FALSE)
  remove_top_level <- top_level[basename(top_level) != "tables"]
  if (length(remove_top_level) > 0) {
    unlink(remove_top_level, recursive = TRUE, force = TRUE)
  }

  invisible(TRUE)
}

## -------------------------------------------------------------------
## Pipeline
## -------------------------------------------------------------------

cat("\n[1/5] Building base CoDA point forecasts\n")
run_script("base_fcast.R", working_dir = CODE_DIR)

cat("\n[2/5] Running point forecast reconciliation and creating point tables\n")
run_script("reco_fcast.R", working_dir = CODE_DIR)

cat("\n[3/5] Generating CoDA bootstrap predictive paths\n")
run_script("gen_boot.R", working_dir = PROJECT_ROOT)

cat("\n[4/5] Reconciling CoDA bootstrap paths and computing CRPS tables\n")
run_script("reco_boot.R", working_dir = PROJECT_ROOT)

cat("\n[5/5] Running residual block-bootstrap robustness check\n")
BB_NSIM <- 10000
BB_CRPS_EVAL_DRAWS <- NULL
BB_BLOCK_SIZE <- 15
BB_CENTER_RESIDUALS <- TRUE
BB_PREFIX <- "BB_"
run_script("residual_block_bootstrap_crps_legacy_fixed.R", working_dir = PROJECT_ROOT)

if (isTRUE(PRUNE_TO_MANUSCRIPT_TABLES)) {
  cat("\nPruning outputs to manuscript CSV/TEX tables only...\n")
  prune_outputs()
}

## -------------------------------------------------------------------
## Final report
## -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Pipeline completed successfully.\n")
cat("Final table outputs: ", file.path(OUT_DIR, "tables"), "\n", sep = "")
cat("============================================================\n\n")

cat("Kept CSV/TEX files:\n")
for (ff in MANUSCRIPT_TABLE_FILES) {
  fpath <- file.path(OUT_DIR, "tables", ff)
  if (file.exists(fpath)) {
    cat(" - outputs/tables/", ff, "\n", sep = "")
  } else {
    cat(" - MISSING after run: outputs/tables/", ff, "\n", sep = "")
  }
}
cat("\n")

