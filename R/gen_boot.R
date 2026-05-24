######################################################################
## Generate bootstrap paths
######################################################################

rm(list = ls())
options(scipen = 999)

PROJECT_ROOT <- getwd()

DATA_DIR <- file.path(PROJECT_ROOT, "data")
LT_DIR <- file.path(PROJECT_ROOT, "LT")
OUT_DIR <- file.path(PROJECT_ROOT, "outputs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "CoDa functions.R"))
source(file.path("R", "coda_boot.R"))

library(tidyverse)
library(compositions)
library(forecast)

load(file.path(DATA_DIR, "ITA_Dx 80+.RData"))

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

boot_paths <- generate_coda_bootstrap_paths(
  dt = dt,
  ih = 15,
  years.fit = 1970:2000,
  lt_dir = LT_DIR,
  nr_rank = 3,
  n_epsilon = 100,
  n_kappa = 100,
  seed = 1234,
  save_file = file.path(OUT_DIR, "bootstrap_coda_paths.rds")
)