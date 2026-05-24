######################################################################
## Generating unreconcilied base forecasts
######################################################################

rm(list = ls())
options(scipen = 999)

PROJECT_ROOT <- getwd()

## If you accidentally run this while already inside R/, go back one level.
if (basename(PROJECT_ROOT) == "R") {
  setwd("..")
  PROJECT_ROOT <- getwd()
}

PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

R_DIR <- file.path(PROJECT_ROOT, "R")
DATA_DIR <- file.path(PROJECT_ROOT, "data")
LT_SOURCE_DIR <- file.path(PROJECT_ROOT, "LT")
OUT_DIR <- file.path(PROJECT_ROOT, "outputs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DATA_FILE <- file.path(DATA_DIR, "Ita_Dx 80+.RData")


required_packages <- c(
  "tidyverse",
  "stats",
  "strucchange",
  "tables",
  "calibrate",
  "compositions",
  "abind",
  "MLmetrics",
  "FactoMineR",
  "tseries",
  "vars",
  "rTensor",
  "forecast",
  "r.jive"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

######################################################################
## Basic checks
######################################################################

if (!file.exists(DATA_FILE)) {
  stop("Cannot find raw data file: ", DATA_FILE, call. = FALSE)
}

if (!dir.exists(LT_SOURCE_DIR)) {
  stop("Cannot find LT folder: ", LT_SOURCE_DIR, call. = FALSE)
}

needed_lt <- c("mltper_5x1.txt", "fltper_5x1.txt", "bltper_5x1.txt")

missing_lt <- needed_lt[
  !file.exists(file.path(LT_SOURCE_DIR, needed_lt))
]

if (length(missing_lt) > 0) {
  stop(
    "Missing LT files in LT/: ",
    paste(missing_lt, collapse = ", "),
    "\nLT folder is: ",
    LT_SOURCE_DIR,
    call. = FALSE
  )
}

if (!file.exists(file.path(R_DIR, "CoDa functions.R"))) {
  stop("Cannot find R/CoDa functions.R", call. = FALSE)
}

if (!file.exists(file.path(R_DIR, "Base function.R"))) {
  stop("Cannot find R/Base function.R", call. = FALSE)
}

file.copy(
  from = file.path(R_DIR, "CoDa functions.R"),
  to = file.path(PROJECT_ROOT, "CoDa functions.r"),
  overwrite = TRUE
)

## Root-level LT folder expected by the original CoDa function

LT_TARGET_DIR <- file.path(PROJECT_ROOT, "LT")
dir.create(LT_TARGET_DIR, recursive = TRUE, showWarnings = FALSE)

for (ff in needed_lt) {
  
  src <- file.path(LT_SOURCE_DIR, ff)
  dst <- file.path(LT_TARGET_DIR, ff)
  
  if (!file.exists(src)) {
    stop("Cannot find LT source file: ", src, call. = FALSE)
  }
  
  src_norm <- normalizePath(src, winslash = "/", mustWork = TRUE)
  
  if (file.exists(dst)) {
    dst_norm <- normalizePath(dst, winslash = "/", mustWork = TRUE)
  } else {
    dst_norm <- normalizePath(dirname(dst), winslash = "/", mustWork = TRUE)
    dst_norm <- file.path(dst_norm, basename(dst))
  }
  
  # On Windows paths are case-insensitive, so compare lower-case paths
  if (tolower(src_norm) == tolower(dst_norm)) {
    message("LT file already in place, skipping copy: ", ff)
    next
  }
  
  ok <- file.copy(
    from = src,
    to = dst,
    overwrite = TRUE
  )
  
  if (!isTRUE(ok)) {
    stop(
      "Could not copy LT file:\n",
      "from: ", src, "\n",
      "to:   ", dst,
      call. = FALSE
    )
  }
}


######################################################################
## Source legacy functions
######################################################################

source(file.path(R_DIR, "CoDa functions.R"))
source(file.path(R_DIR, "Base function.R"))

######################################################################
## Load raw data
######################################################################

load(DATA_FILE)

if (!exists("ITA_Dx")) {
  stop("The file data/Ita_Dx 80+.RData does not contain object ITA_Dx.", call. = FALSE)
}

######################################################################
## Build input dt exactly as in the original code
######################################################################

dt <- rbind(
  ITA_Dx,
  ITA_Dx %>%
    group_by(Year, Age, Cause_Rev) %>%
    summarise(Dx = sum(Dx), .groups = "drop") %>%
    mutate(Sex = "T")
) %>%
  arrange(Cause_Rev, Sex, Age, Year) %>%
  pivot_wider(names_from = Age, values_from = Dx) %>%
  filter(Year < 2016)

######################################################################
## Parameters
######################################################################

baseline <- 1970:2000
ih <- 15

######################################################################
## Run CoDa models
######################################################################

message("Generating bottom-level CoDA forecasts: M and F by cause...")

Output_bl_M <- CoDa(dt, sex = "M", hstr = 1, ih = ih, years.fit = baseline)
Output_bl_F <- CoDa(dt, sex = "F", hstr = 1, ih = ih, years.fit = baseline)

message("Generating middle-level CoDA forecasts: M and F all-cause...")

Output_ml_M <- CoDa(dt, sex = "M", hstr = 2, ih = ih, years.fit = baseline)
Output_ml_F <- CoDa(dt, sex = "F", hstr = 2, ih = ih, years.fit = baseline)

message("Generating middle-level CoDA forecasts: T by cause...")

Output_ml_T <- CoDa(dt, sex = "T", hstr = 1, ih = ih, years.fit = baseline)

message("Generating top-level CoDA forecasts: T all-cause...")

Output_tl <- CoDa(dt, sex = "T", hstr = 2, ih = ih, years.fit = baseline)

######################################################################
## Merge results exactly as in the original code
######################################################################

Output_bl <- rbind(
  Output_bl_M %>%
    mutate(Sex = "M", .after = Age),
  Output_bl_F %>%
    mutate(Sex = "F", .after = Age)
) %>%
  arrange(Type, Sex, Year, Cause_cod, Age) %>%
  mutate(Year = as.numeric(Year))

Output_ml <- rbind(
  Output_ml_M %>%
    mutate(Sex = "M", .after = Age),
  Output_ml_F %>%
    mutate(Sex = "F", .after = Age),
  Output_ml_T %>%
    mutate(Sex = "T", .after = Age)
)

Output_tl <- Output_tl %>%
  mutate(Sex = "T", .after = Age)

## Original legacy correction
Output_ml[, 9] <- Output_ml[, 6]
Output_ml <- Output_ml[, -6]

Output <- rbind(
  Output_bl,
  Output_ml,
  Output_tl
) %>%
  arrange(Type, Year, Sex, Cause_cod, Age)

######################################################################
## Save Output
######################################################################

## Main generated Output
saveRDS(
  Output,
  file = file.path(OUT_DIR, "Output_generated_legacy.rds")
)

save(
  Output,
  file = file.path(OUT_DIR, "Output_generated_legacy.RData")
)

## Standard file used by the next scripts: point forecasting, tables, figures
forecast_results <- Output

saveRDS(
  forecast_results,
  file = file.path(OUT_DIR, "forecast_results.rds")
)

save(
  forecast_results,
  file = file.path(OUT_DIR, "forecast_results.RData")
)

message("\nDone.")
message("Output saved as:")
message(file.path(OUT_DIR, "Output_generated_legacy.rds"))
message(file.path(OUT_DIR, "Output_generated_legacy.RData"))

message("\nForecast results saved for the next scripts as:")
message(file.path(OUT_DIR, "forecast_results.rds"))
message(file.path(OUT_DIR, "forecast_results.RData"))