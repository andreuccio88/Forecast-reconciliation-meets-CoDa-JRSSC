######################################################################
## Point forecasts: reconciliation, RMSE/MAE tables, and official figures.
######################################################################

rm(list = ls())
options(scipen = 999)

######################################################################
## Project paths
######################################################################

PROJECT_ROOT <- getwd()

## If accidentally run from inside R/, go back to project root.
if (basename(PROJECT_ROOT) == "R") {
  setwd("..")
  PROJECT_ROOT <- getwd()
}

R_DIR <- file.path(PROJECT_ROOT, "R")
OUT_DIR <- file.path(PROJECT_ROOT, "outputs")
RDS_DIR <- file.path(OUT_DIR, "rds")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

######################################################################
## Packages
######################################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "readr",
  "FoReco",
  "patchwork",
  "xtable"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

######################################################################
## Helpers
######################################################################

source(file.path(R_DIR, "helpers.R"))

######################################################################
## Load generated Output
######################################################################

OUTPUT_FILE <- file.path(OUT_DIR, "forecast_results.rds")

if (!file.exists(OUTPUT_FILE)) {
  stop(
    "Cannot find generated Output file: ",
    OUTPUT_FILE,
    "\nRun R/00_generate_output.R first.",
    call. = FALSE
  )
}

Output <- readRDS(OUTPUT_FILE) |>
  dplyr::filter(Year >= 1970)

ages <- unique(Output$Age)
causes <- unique(Output$Cause)
cause_specific <- causes[1:8]
cause_total <- causes[9]

stopifnot(length(cause_specific) == 8)
stopifnot(!is.na(cause_total))

names <- build_series_names()
S_g1 <- build_S_g1(names$rows_g1, names$bottom_by_sex)
S_g2 <- build_S_g2(names$rows_g2, names$bottom_interlaced)
S_gT <- build_S_gT(names$rows_gT, names$bottom_interlaced)

res_g1 <- run_point_structure(Output, ages, cause_specific, cause_total, names$rows_g1, S_g1, 1:3,  "g1_sex")
res_g2 <- run_point_structure(Output, ages, cause_specific, cause_total, names$rows_g2, S_g2, 1:9,  "g2_cause")
res_gT <- run_point_structure(Output, ages, cause_specific, cause_total, names$rows_gT, S_gT, 1:11, "gT_full")

saveRDS(
  list(res_g1 = res_g1, res_g2 = res_g2, res_gT = res_gT, names = names, ages = ages, causes = causes),
  file.path(RDS_DIR, "point_forecasts.rds")
)

rmse_tabs <- build_point_tables(res_gT, res_g1, res_g2, names, "RMSE")
mae_tabs  <- build_point_tables(res_gT, res_g1, res_g2, names, "MAE")

readr::write_csv(rmse_tabs$gT,   file.path(TABLE_DIR, "RMSE_full_grouped.csv"))
readr::write_csv(rmse_tabs$g1,   file.path(TABLE_DIR, "RMSE_sex_hierarchy.csv"))
readr::write_csv(rmse_tabs$g2,   file.path(TABLE_DIR, "RMSE_cause_hierarchy.csv"))
readr::write_csv(rmse_tabs$full, file.path(TABLE_DIR, "RMSE_all_structures.csv"))
readr::write_csv(mae_tabs$gT,    file.path(TABLE_DIR, "MAE_full_grouped.csv"))
readr::write_csv(mae_tabs$g1,    file.path(TABLE_DIR, "MAE_sex_hierarchy.csv"))
readr::write_csv(mae_tabs$g2,    file.path(TABLE_DIR, "MAE_cause_hierarchy.csv"))
readr::write_csv(mae_tabs$full,  file.path(TABLE_DIR, "MAE_all_structures.csv"))

write_manuscript_table(rmse_tabs$gT, rmse_tabs$g1, rmse_tabs$g2, "RMSE", file.path(TABLE_DIR, "point_RMSE_table.tex"))
write_manuscript_table(mae_tabs$gT,  mae_tabs$g1,  mae_tabs$g2,  "MAE",  file.path(TABLE_DIR, "point_MAE_table.tex"))

######################################################################
## Figures
######################################################################

## 1. Top-level RMSE distribution

df_top <- dplyr::bind_rows(
  data.frame(Age = ages, Method = "CoDA", Type = "Total", RMSE = as.numeric(res_gT$RMSE$CoDA["Total", ])),
  data.frame(Age = ages, Method = "OLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$OLS["Total", ])),
  data.frame(Age = ages, Method = "WLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$WLS["Total", ])),
  data.frame(Age = ages, Method = "MinT", Type = "Total", RMSE = as.numeric(res_gT$RMSE$MinT["Total", ])),
  data.frame(Age = ages, Method = "CoDA", Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$CoDA["Total", ])),
  data.frame(Age = ages, Method = "OLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$OLS["Total", ])),
  data.frame(Age = ages, Method = "WLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$WLS["Total", ])),
  data.frame(Age = ages, Method = "MinT", Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$MinT["Total", ])),
  data.frame(Age = ages, Method = "CoDA", Type = "Cause", RMSE = as.numeric(res_g2$RMSE$CoDA["Total", ])),
  data.frame(Age = ages, Method = "OLS",  Type = "Cause", RMSE = as.numeric(res_g2$RMSE$OLS["Total", ])),
  data.frame(Age = ages, Method = "WLS",  Type = "Cause", RMSE = as.numeric(res_g2$RMSE$WLS["Total", ])),
  data.frame(Age = ages, Method = "MinT", Type = "Cause", RMSE = as.numeric(res_g2$RMSE$MinT["Total", ]))
) |>
  dplyr::mutate(
    Type = factor(Type, levels = c("Total", "Sex", "Cause")),
    Method = factor(Method, levels = c("CoDA", "MinT", "OLS", "WLS"))
  )

p_top <- ggplot2::ggplot(df_top, ggplot2::aes(x = Type, y = RMSE, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::coord_cartesian(ylim = c(0, 17000)) +
  ggplot2::labs(x = "Grouped structure", y = "RMSE") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "right")
save_plot(p_top, "total-rmse", FIGURE_DIR, 8.5, 5.5)

## 2. Heatmap top-level RSE
p_heat <- build_heatmap_top(res_gT, res_g1, res_g2, ages, TABLE_DIR, FIGURE_DIR)
print(p_heat)

## 3. Middle-level sex totals

df_middle_sex <- dplyr::bind_rows(
  data.frame(Age = ages, Sex = "Male",   Method = "CoDA", Type = "CoDA",  RMSE = as.numeric(res_g1$RMSE$CoDA["Total - Males", ])),
  data.frame(Age = ages, Sex = "Female", Method = "CoDA", Type = "CoDA",  RMSE = as.numeric(res_g1$RMSE$CoDA["Total - Females", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "OLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$OLS["Total - Males", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "WLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$WLS["Total - Males", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "MinT", Type = "Total", RMSE = as.numeric(res_gT$RMSE$MinT["Total - Males", ])),
  data.frame(Age = ages, Sex = "Female", Method = "OLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$OLS["Total - Females", ])),
  data.frame(Age = ages, Sex = "Female", Method = "WLS",  Type = "Total", RMSE = as.numeric(res_gT$RMSE$WLS["Total - Females", ])),
  data.frame(Age = ages, Sex = "Female", Method = "MinT", Type = "Total", RMSE = as.numeric(res_gT$RMSE$MinT["Total - Females", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "OLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$OLS["Total - Males", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "WLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$WLS["Total - Males", ])),
  data.frame(Age = ages, Sex = "Male",   Method = "MinT", Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$MinT["Total - Males", ])),
  data.frame(Age = ages, Sex = "Female", Method = "OLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$OLS["Total - Females", ])),
  data.frame(Age = ages, Sex = "Female", Method = "WLS",  Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$WLS["Total - Females", ])),
  data.frame(Age = ages, Sex = "Female", Method = "MinT", Type = "Sex",   RMSE = as.numeric(res_g1$RMSE$MinT["Total - Females", ]))
) |>
  dplyr::mutate(
    Sex = factor(Sex, levels = c("Female", "Male")),
    Method = factor(Method, levels = c("CoDA", "MinT", "OLS", "WLS")),
    Type = factor(Type, levels = c("CoDA", "Sex", "Total"))
  )

readr::write_csv(
  df_middle_sex |>
    dplyr::group_by(Sex, Method, Type) |>
    dplyr::summarise(median_RMSE = median(RMSE, na.rm = TRUE), .groups = "drop"),
  file.path(TABLE_DIR, "diagnostic_middle_sex_plot_medians.csv")
)

p_middle_sex <- ggplot2::ggplot(df_middle_sex, ggplot2::aes(x = Method, y = RMSE, fill = Type)) +
  ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(0.75)) +
  ggplot2::coord_cartesian(ylim = c(0, 7500)) +
  ggplot2::facet_wrap(~Sex) +
  ggplot2::labs(x = "Method", y = "RMSE") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
save_plot(p_middle_sex, "middle-rmse", FIGURE_DIR, 8.5, 5.5)

## 4. Middle-level cause totals, median by cause

df_middle_cause_median <- dplyr::bind_rows(
  data.frame(Cause = names$cause_totals, Type = "CoDA",        Method = "CoDA", RMSE = apply(res_g2$RMSE$CoDA[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Total",       Method = "OLS",  RMSE = apply(res_gT$RMSE$OLS[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Total",       Method = "WLS",  RMSE = apply(res_gT$RMSE$WLS[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Total",       Method = "MinT", RMSE = apply(res_gT$RMSE$MinT[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Cause-based", Method = "OLS",  RMSE = apply(res_g2$RMSE$OLS[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Cause-based", Method = "WLS",  RMSE = apply(res_g2$RMSE$WLS[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Cause = names$cause_totals, Type = "Cause-based", Method = "MinT", RMSE = apply(res_g2$RMSE$MinT[names$cause_totals, , drop = FALSE], 1, median, na.rm = TRUE))
) |>
  dplyr::mutate(
    Type = factor(Type, levels = c("CoDA", "Total", "Cause-based")),
    Method = factor(Method, levels = c("CoDA", "MinT", "OLS", "WLS"))
  )

p_middle_cause_median <- ggplot2::ggplot(df_middle_cause_median, ggplot2::aes(x = Type, y = RMSE, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge2(preserve = "single")) +
  ggplot2::coord_cartesian(ylim = c(0, 750)) +
  ggplot2::labs(x = "Grouped structure", y = "RMSE") +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
save_plot(p_middle_cause_median, "middle2-rmse-tot", FIGURE_DIR, 8.5, 5.5)

## 5. Middle-level cause totals by cause, across ages

df_middle_cause_age <- dplyr::bind_rows(
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Grouped structure (Total)", Method = "OLS",  RMSE = row_vec(res_gT$RMSE$OLS[names$cause_totals, , drop = FALSE])),
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Grouped structure (Total)", Method = "WLS",  RMSE = row_vec(res_gT$RMSE$WLS[names$cause_totals, , drop = FALSE])),
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Grouped structure (Total)", Method = "MinT", RMSE = row_vec(res_gT$RMSE$MinT[names$cause_totals, , drop = FALSE])),
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Cause-based hierarchy",  Method = "OLS",  RMSE = row_vec(res_g2$RMSE$OLS[names$cause_totals, , drop = FALSE])),
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Cause-based hierarchy",  Method = "WLS",  RMSE = row_vec(res_g2$RMSE$WLS[names$cause_totals, , drop = FALSE])),
  data.frame(Age = rep(ages, times = 8), Cause = rep(paste0("Cause ", 1:8), each = length(ages)), Type = "Cause-based hierarchy",  Method = "MinT", RMSE = row_vec(res_g2$RMSE$MinT[names$cause_totals, , drop = FALSE]))
) |>
  dplyr::mutate(
    Cause = factor(Cause, levels = paste0("Cause ", 1:8)),
    Method = factor(Method, levels = c("MinT", "OLS", "WLS"))
  )

p_middle_cause_age <- ggplot2::ggplot(df_middle_cause_age, ggplot2::aes(x = Cause, y = RMSE, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::coord_cartesian(ylim = c(0, 5000)) +
  ggplot2::facet_wrap(~Type, nrow = 1) +
  ggplot2::labs(x = "Cause of death", y = "RMSE") +
  ggplot2::theme_minimal()
save_plot(p_middle_cause_age, "middle2-rmse", FIGURE_DIR, 11, 5.5)

## 6. Bottom-level series

df_bottom <- dplyr::bind_rows(
  data.frame(Series = names$bottom_interlaced, Type = "CoDA",        Method = "CoDA", RMSE = apply(res_gT$RMSE$CoDA[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Total",       Method = "OLS",  RMSE = apply(res_gT$RMSE$OLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Total",       Method = "WLS",  RMSE = apply(res_gT$RMSE$WLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Total",       Method = "MinT", RMSE = apply(res_gT$RMSE$MinT[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Sex",         Method = "OLS",  RMSE = apply(res_g1$RMSE$OLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Sex",         Method = "WLS",  RMSE = apply(res_g1$RMSE$WLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Sex",         Method = "MinT", RMSE = apply(res_g1$RMSE$MinT[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Cause-based", Method = "OLS",  RMSE = apply(res_g2$RMSE$OLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Cause-based", Method = "WLS",  RMSE = apply(res_g2$RMSE$WLS[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE)),
  data.frame(Series = names$bottom_interlaced, Type = "Cause-based", Method = "MinT", RMSE = apply(res_g2$RMSE$MinT[names$bottom_interlaced, , drop = FALSE], 1, median, na.rm = TRUE))
) |>
  dplyr::mutate(
    Type = factor(Type, levels = c("CoDA", "Total", "Sex", "Cause-based")),
    Method = factor(Method, levels = c("CoDA", "MinT", "OLS", "WLS"))
  )

p_bottom <- ggplot2::ggplot(df_bottom, ggplot2::aes(x = Type, y = RMSE, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge2(preserve = "single")) +
  ggplot2::coord_cartesian(ylim = c(0, 500)) +
  ggplot2::labs(x = "Grouped structure", y = "RMSE") +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(legend.title = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
save_plot(p_bottom, "bottom", FIGURE_DIR, 8.5, 5.5)

cat("\nPoint forecasts completed. Tables and figures saved in outputs/.\n")
