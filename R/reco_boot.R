######################################################################
##
## Reconcile CoDA bootstrap paths and compute median CRPS tables.
##
######################################################################

rm(list = ls())
options(scipen = 999)

######################################################################
## Project paths
######################################################################

PROJECT_ROOT <- getwd()

if (basename(PROJECT_ROOT) == "R") {
  setwd("..")
  PROJECT_ROOT <- getwd()
}

R_DIR <- file.path(PROJECT_ROOT, "R")
OUT_DIR <- file.path(PROJECT_ROOT, "outputs")
RDS_DIR <- file.path(OUT_DIR, "rds")
TABLE_DIR <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

BOOT_FILE <- file.path(OUT_DIR, "bootstrap_coda_paths.rds")

if (!file.exists(BOOT_FILE)) {
  stop(
    "Cannot find bootstrap file: ", BOOT_FILE,
    "\nRun R/02_generate_bootstrap_paths.R first.",
    call. = FALSE
  )
}

######################################################################
## Packages
######################################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "FoReco"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_packages, library, character.only = TRUE))

######################################################################
## Load helpers
######################################################################

if (!file.exists(file.path(R_DIR, "helpers.R"))) {
  stop("Cannot find R/helpers.R", call. = FALSE)
}

source(file.path(R_DIR, "helpers.R"))

######################################################################
## Local fallbacks, in case helpers.R does not contain them
######################################################################

if (!exists("make_metric_matrix")) {
  
  make_metric_matrix <- function(base, ols, wls, mint, series_order) {
    data.frame(
      Series = series_order,
      CoDA = as.numeric(apply(base[series_order, , drop = FALSE], 1, median, na.rm = TRUE)),
      OLS  = as.numeric(apply(ols[series_order, , drop = FALSE], 1, median, na.rm = TRUE)),
      WLS  = as.numeric(apply(wls[series_order, , drop = FALSE], 1, median, na.rm = TRUE)),
      MinT = as.numeric(apply(mint[series_order, , drop = FALSE], 1, median, na.rm = TRUE)),
      check.names = FALSE
    )
  }
}

if (!exists("latex_num")) {
  
  latex_num <- function(x, is_best = FALSE) {
    out <- sprintf("%.2f", x)
    ifelse(is_best, paste0("\\textbf{", out, "}"), out)
  }
}

if (!exists("latex_metric_row")) {
  
  latex_metric_row <- function(series, vals) {
    vals <- as.numeric(vals)
    best <- seq_along(vals) == which.min(vals)
    paste0(
      "& ", series, " & ",
      paste(latex_num(vals, best), collapse = " & "),
      " \\\\"
    )
  }
}

if (!exists("write_manuscript_table")) {
  
  write_manuscript_table <- function(gT, g1, g2, metric = c("RMSE", "MAE", "CRPS"), file) {
    
    metric <- match.arg(metric)
    numcols <- c("CoDA", "OLS", "WLS", "MinT")
    
    caption <- paste0(
      "\\small Median forecast accuracy results: ",
      metric,
      " computed for the three levels of the grouped and hierarchical structures. ",
      "The best method is highlighted in bold."
    )
    
    label <- paste0("tab:", tolower(metric), "_median")
    
    lines <- c(
      "\\begin{center}",
      "\\tabcolsep 0.315in",
      "\\renewcommand\\arraystretch{0.98}",
      "\\begin{small}",
      "\\begin{longtable}{@{}llrrrr@{}}",
      paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
      "\\toprule",
      "Hierarchy & Series & CoDA & OLS & WLS & MinT \\\\",
      "\\midrule",
      "\\endfirsthead",
      "\\toprule",
      "Hierarchy & Series & CoDA & OLS & WLS & MinT \\\\",
      "\\midrule",
      "\\endhead",
      "\\hline \\multicolumn{6}{r}{{Continued on next page}} \\\\",
      "\\endfoot",
      "\\endlastfoot",
      
      "Grouped & \\multicolumn{5}{l}{\\underline{Top-level:}}\\\\",
      latex_metric_row("Total", gT[1, numcols]),
      "\\cmidrule{2-6}",
      "& \\multicolumn{5}{l}{\\underline{Middle-level:}}\\\\"
    )
    
    for (ii in 2:9) {
      lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))
    }
    
    for (ii in 10:11) {
      lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))
    }
    
    lines <- c(
      lines,
      "\\cmidrule{2-6}",
      "& \\multicolumn{5}{l}{\\underline{Bottom-level:}}\\\\"
    )
    
    for (ii in 12:27) {
      lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))
    }
    
    lines <- c(
      lines,
      "\\midrule",
      "Sex & \\multicolumn{5}{l}{\\underline{Top-level:}}\\\\",
      paste0("$\\downarrow$ ", substring(latex_metric_row("Total", g1[1, numcols]), 2)),
      "\\cmidrule{2-6}",
      "Cause & \\multicolumn{5}{l}{\\underline{Middle-level:}}\\\\",
      latex_metric_row(g1$Series[2], g1[2, numcols]),
      latex_metric_row(g1$Series[3], g1[3, numcols]),
      "\\cmidrule{2-6}",
      "& \\multicolumn{5}{l}{\\underline{Bottom-level:}}\\\\"
    )
    
    for (ii in 4:19) {
      lines <- c(lines, latex_metric_row(g1$Series[ii], g1[ii, numcols]))
    }
    
    lines <- c(
      lines,
      "\\midrule",
      "Cause & \\multicolumn{5}{l}{\\underline{Top-level:}}\\\\",
      paste0("$\\downarrow$ ", substring(latex_metric_row("Total", g2[1, numcols]), 2)),
      "\\cmidrule{2-6}",
      "Sex & \\multicolumn{5}{l}{\\underline{Middle-level:}}\\\\"
    )
    
    for (ii in 2:9) {
      lines <- c(lines, latex_metric_row(g2$Series[ii], g2[ii, numcols]))
    }
    
    lines <- c(
      lines,
      "\\cmidrule{2-6}",
      "& \\multicolumn{5}{l}{\\underline{Bottom-level:}}\\\\"
    )
    
    for (ii in 10:25) {
      lines <- c(lines, latex_metric_row(g2$Series[ii], g2[ii, numcols]))
    }
    
    lines <- c(
      lines,
      "\\bottomrule",
      "\\end{longtable}",
      "\\end{small}",
      "\\end{center}"
    )
    
    writeLines(lines, file)
  }
}

######################################################################
## S matrices, same order as point forecasts
######################################################################

build_boot_series_names <- function() {
  
  bottom_by_sex <- c(
    paste0("Cause ", 1:8, " - Males"),
    paste0("Cause ", 1:8, " - Females")
  )
  
  bottom_interlaced <- as.vector(rbind(
    paste0("Cause ", 1:8, " - Males"),
    paste0("Cause ", 1:8, " - Females")
  ))
  
  cause_totals <- paste0("Cause ", 1:8, " - Total")
  
  rows_g1 <- c(
    "Total",
    "Total - Males",
    "Total - Females",
    bottom_by_sex
  )
  
  rows_g2 <- c(
    "Total",
    cause_totals,
    bottom_interlaced
  )
  
  rows_gT <- c(
    "Total",
    cause_totals,
    "Total - Males",
    "Total - Females",
    bottom_interlaced
  )
  
  list(
    bottom_by_sex = bottom_by_sex,
    bottom_interlaced = bottom_interlaced,
    cause_totals = cause_totals,
    rows_g1 = rows_g1,
    rows_g2 = rows_g2,
    rows_gT = rows_gT
  )
}

if (!exists("build_S_g1")) {
  
  build_S_g1 <- function(rows_g1, bottom_by_sex) {
    S <- matrix(0, nrow = length(rows_g1), ncol = length(bottom_by_sex))
    rownames(S) <- rows_g1
    colnames(S) <- bottom_by_sex
    
    S["Total", ] <- 1
    S["Total - Males", paste0("Cause ", 1:8, " - Males")] <- 1
    S["Total - Females", paste0("Cause ", 1:8, " - Females")] <- 1
    S[bottom_by_sex, bottom_by_sex] <- diag(length(bottom_by_sex))
    
    S
  }
}

if (!exists("build_S_g2")) {
  
  build_S_g2 <- function(rows_g2, bottom_interlaced) {
    S <- matrix(0, nrow = length(rows_g2), ncol = length(bottom_interlaced))
    rownames(S) <- rows_g2
    colnames(S) <- bottom_interlaced
    
    S["Total", ] <- 1
    
    for (i in 1:8) {
      S[paste0("Cause ", i, " - Total"),
        c(paste0("Cause ", i, " - Males"),
          paste0("Cause ", i, " - Females"))] <- 1
    }
    
    S[bottom_interlaced, bottom_interlaced] <- diag(length(bottom_interlaced))
    S
  }
}

if (!exists("build_S_gT")) {
  
  build_S_gT <- function(rows_gT, bottom_interlaced) {
    S <- matrix(0, nrow = length(rows_gT), ncol = length(bottom_interlaced))
    rownames(S) <- rows_gT
    colnames(S) <- bottom_interlaced
    
    S["Total", ] <- 1
    
    for (i in 1:8) {
      S[paste0("Cause ", i, " - Total"),
        c(paste0("Cause ", i, " - Males"),
          paste0("Cause ", i, " - Females"))] <- 1
    }
    
    S["Total - Males", paste0("Cause ", 1:8, " - Males")] <- 1
    S["Total - Females", paste0("Cause ", 1:8, " - Females")] <- 1
    
    S[bottom_interlaced, bottom_interlaced] <- diag(length(bottom_interlaced))
    S
  }
}

######################################################################
## CRPS helpers
######################################################################

crps_empirical_one <- function(y, x) {
  ## y: scalar observed value
  ## x: vector of predictive draws
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (!is.finite(y) || length(x) < 2) {
    return(NA_real_)
  }
  
  B <- length(x)
  xs <- sort(x)
  
  term1 <- mean(abs(xs - y))
  
  ## Efficient pairwise absolute distance:
  ## mean_{i,j}|X_i - X_j|
  pair_sum <- 2 * sum((2 * seq_len(B) - B - 1) * xs)
  term2 <- pair_sum / (B^2)
  
  term1 - 0.5 * term2
}


crps_matrix_from_draws <- function(actual, draws) {
  ## actual: H x n_series matrix
  ## draws:  B x H x n_series array
  ## returns: n_series vector, mean CRPS over horizons
  
  if (!is.matrix(actual)) {
    stop("actual must be a matrix.", call. = FALSE)
  }
  
  if (length(dim(draws)) != 3) {
    stop("draws must be a B x H x n_series array.", call. = FALSE)
  }
  
  B <- dim(draws)[1]
  H <- dim(draws)[2]
  n_series <- dim(draws)[3]
  
  if (nrow(actual) != H || ncol(actual) != n_series) {
    stop("Dimensions of actual and draws are inconsistent.", call. = FALSE)
  }
  
  out <- numeric(n_series)
  names(out) <- colnames(actual)
  
  for (j in seq_len(n_series)) {
    crps_h <- numeric(H)
    
    for (h in seq_len(H)) {
      crps_h[h] <- crps_empirical_one(
        y = actual[h, j],
        x = draws[, h, j]
      )
    }
    
    out[j] <- mean(crps_h, na.rm = TRUE)
  }
  
  out
}


compute_crps_age_matrix <- function(actual_list, draws_list, series_names) {
  ## actual_list: list by age, each H x n_series
  ## draws_list: list by age, each B x H x n_series
  ## output: n_series x n_ages matrix
  
  n_ages <- length(draws_list)
  out <- matrix(
    NA_real_,
    nrow = length(series_names),
    ncol = n_ages,
    dimnames = list(series_names, names(draws_list))
  )
  
  for (a in seq_along(draws_list)) {
    out[, a] <- crps_matrix_from_draws(
      actual = actual_list[[a]],
      draws = draws_list[[a]]
    )[series_names]
    
    message("CRPS completed for age ", names(draws_list)[a])
  }
  
  out
}

######################################################################
## Bootstrap reconciliation
######################################################################

reconcile_one_draw_matrix <- function(base_mat, S, agg_rows, comb, residuals) {
  ## base_mat: H x n_series
  ## S: full summing matrix
  ## agg_rows: aggregation rows used by FoReco::csrec
  ## residuals: in-sample residual matrix, T x n_series
  
  out <- FoReco::csrec(
    base = base_mat,
    agg_mat = S[agg_rows, , drop = FALSE],
    comb = comb,
    nn = "strc_osqp",
    res = residuals
  )
  
  out <- as.matrix(out)
  colnames(out) <- colnames(base_mat)
  out
}


reconcile_bootstrap_array <- function(draws, S, agg_rows, comb, residuals) {
  ## draws: B x H x n_series
  ## returns B x H x n_series
  
  B <- dim(draws)[1]
  H <- dim(draws)[2]
  n_series <- dim(draws)[3]
  
  out <- array(
    NA_real_,
    dim = dim(draws),
    dimnames = dimnames(draws)
  )
  
  for (b in seq_len(B)) {
    
    base_mat <- draws[b, , , drop = FALSE]
    base_mat <- matrix(
      base_mat,
      nrow = H,
      ncol = n_series,
      dimnames = list(dimnames(draws)[[2]], dimnames(draws)[[3]])
    )
    
    rec_mat <- reconcile_one_draw_matrix(
      base_mat = base_mat,
      S = S,
      agg_rows = agg_rows,
      comb = comb,
      residuals = residuals
    )
    
    out[b, , ] <- rec_mat
    
    if (b %% 100 == 0) {
      message("  reconciled bootstrap path ", b, " / ", B)
    }
  }
  
  out
}


reconcile_structure <- function(draws_list, residuals_list, S, agg_rows, label) {
  
  out <- list(
    CoDA = draws_list,
    OLS = vector("list", length(draws_list)),
    WLS = vector("list", length(draws_list)),
    MinT = vector("list", length(draws_list))
  )
  
  names(out$OLS) <- names(out$WLS) <- names(out$MinT) <- names(draws_list)
  
  for (a in seq_along(draws_list)) {
    
    message("\nReconciling ", label, " | age ", names(draws_list)[a], " | OLS")
    out$OLS[[a]] <- reconcile_bootstrap_array(
      draws = draws_list[[a]],
      S = S,
      agg_rows = agg_rows,
      comb = "ols",
      residuals = residuals_list[[a]]
    )
    
    message("\nReconciling ", label, " | age ", names(draws_list)[a], " | WLS")
    out$WLS[[a]] <- reconcile_bootstrap_array(
      draws = draws_list[[a]],
      S = S,
      agg_rows = agg_rows,
      comb = "wls",
      residuals = residuals_list[[a]]
    )
    
    message("\nReconciling ", label, " | age ", names(draws_list)[a], " | MinT")
    out$MinT[[a]] <- reconcile_bootstrap_array(
      draws = draws_list[[a]],
      S = S,
      agg_rows = agg_rows,
      comb = "shr",
      residuals = residuals_list[[a]]
    )
  }
  
  out
}

######################################################################
## Load bootstrap paths
######################################################################

boot <- readRDS(BOOT_FILE)

nm <- if (!is.null(boot$names)) boot$names else build_boot_series_names()

S_g1 <- build_S_g1(nm$rows_g1, nm$bottom_by_sex)
S_g2 <- build_S_g2(nm$rows_g2, nm$bottom_interlaced)
S_gT <- build_S_gT(nm$rows_gT, nm$bottom_interlaced)

######################################################################
## Reconcile all bootstrap paths
######################################################################

message("\nStarting bootstrap reconciliation. This may take a while, because apparently we chose civilisation over closed-form shortcuts.\n")

rec_g1 <- reconcile_structure(
  draws_list = boot$g1$draws,
  residuals_list = boot$g1$residuals,
  S = S_g1,
  agg_rows = 1:3,
  label = "g1 sex hierarchy"
)

rec_g2 <- reconcile_structure(
  draws_list = boot$g2$draws,
  residuals_list = boot$g2$residuals,
  S = S_g2,
  agg_rows = 1:9,
  label = "g2 cause hierarchy"
)

rec_gT <- reconcile_structure(
  draws_list = boot$gT$draws,
  residuals_list = boot$gT$residuals,
  S = S_gT,
  agg_rows = 1:11,
  label = "gT full grouped"
)

######################################################################
## Compute CRPS matrices: series x age
######################################################################

message("\nComputing CRPS matrices.\n")

CRPS_gT <- list(
  CoDA = compute_crps_age_matrix(boot$gT$actual, rec_gT$CoDA, nm$rows_gT),
  OLS  = compute_crps_age_matrix(boot$gT$actual, rec_gT$OLS,  nm$rows_gT),
  WLS  = compute_crps_age_matrix(boot$gT$actual, rec_gT$WLS,  nm$rows_gT),
  MinT = compute_crps_age_matrix(boot$gT$actual, rec_gT$MinT, nm$rows_gT)
)

CRPS_g1 <- list(
  CoDA = compute_crps_age_matrix(boot$g1$actual, rec_g1$CoDA, nm$rows_g1),
  OLS  = compute_crps_age_matrix(boot$g1$actual, rec_g1$OLS,  nm$rows_g1),
  WLS  = compute_crps_age_matrix(boot$g1$actual, rec_g1$WLS,  nm$rows_g1),
  MinT = compute_crps_age_matrix(boot$g1$actual, rec_g1$MinT, nm$rows_g1)
)

CRPS_g2 <- list(
  CoDA = compute_crps_age_matrix(boot$g2$actual, rec_g2$CoDA, nm$rows_g2),
  OLS  = compute_crps_age_matrix(boot$g2$actual, rec_g2$OLS,  nm$rows_g2),
  WLS  = compute_crps_age_matrix(boot$g2$actual, rec_g2$WLS,  nm$rows_g2),
  MinT = compute_crps_age_matrix(boot$g2$actual, rec_g2$MinT, nm$rows_g2)
)

######################################################################
## Build CRPS tables, same structure as RMSE
######################################################################

g1_display <- c(
  "Total",
  "Total - Males",
  "Total - Females",
  nm$bottom_interlaced
)

crps_gT <- make_metric_matrix(
  CRPS_gT$CoDA,
  CRPS_gT$OLS,
  CRPS_gT$WLS,
  CRPS_gT$MinT,
  nm$rows_gT
)

crps_g1 <- make_metric_matrix(
  CRPS_g1$CoDA,
  CRPS_g1$OLS,
  CRPS_g1$WLS,
  CRPS_g1$MinT,
  g1_display
)

crps_g2 <- make_metric_matrix(
  CRPS_g2$CoDA,
  CRPS_g2$OLS,
  CRPS_g2$WLS,
  CRPS_g2$MinT,
  nm$rows_g2
)

add_metadata <- function(df, hierarchy, levels) {
  df |>
    dplyr::mutate(
      Hierarchy = hierarchy,
      Level = levels,
      .before = Series
    )
}

crps_full <- dplyr::bind_rows(
  add_metadata(
    crps_gT,
    "Full grouped structure",
    c("Top", rep("Middle (cause)", 8), rep("Middle (sex)", 2), rep("Bottom", 16))
  ),
  add_metadata(
    crps_g1,
    "Sex-based hierarchy",
    c("Top", rep("Middle (sex)", 2), rep("Bottom", 16))
  ),
  add_metadata(
    crps_g2,
    "Cause-based hierarchy",
    c("Top", rep("Middle (cause)", 8), rep("Bottom", 16))
  )
)

######################################################################
## Save outputs
######################################################################

readr::write_csv(crps_gT, file.path(TABLE_DIR, "CRPS_full_grouped.csv"))
readr::write_csv(crps_g1, file.path(TABLE_DIR, "CRPS_sex_hierarchy.csv"))
readr::write_csv(crps_g2, file.path(TABLE_DIR, "CRPS_cause_hierarchy.csv"))
readr::write_csv(crps_full, file.path(TABLE_DIR, "CRPS_all_structures.csv"))




## -------------------------------------------------------------------
## CRPS ratios relative to CoDA
## -------------------------------------------------------------------

add_crps_ratios <- function(df) {
  df |>
    dplyr::mutate(
      OLS_ratio  = OLS / CoDA,
      WLS_ratio  = WLS / CoDA,
      MinT_ratio = MinT / CoDA
    )
}

crps_gT_ratio <- add_crps_ratios(crps_gT)
crps_g1_ratio <- add_crps_ratios(crps_g1)
crps_g2_ratio <- add_crps_ratios(crps_g2)
crps_full_ratio <- add_crps_ratios(crps_full)

readr::write_csv(
  crps_gT_ratio,
  file.path(TABLE_DIR, "CRPS_full_grouped_ratios.csv")
)

readr::write_csv(
  crps_g1_ratio,
  file.path(TABLE_DIR, "CRPS_sex_hierarchy_ratios.csv")
)

readr::write_csv(
  crps_g2_ratio,
  file.path(TABLE_DIR, "CRPS_cause_hierarchy_ratios.csv")
)

readr::write_csv(
  crps_full_ratio,
  file.path(TABLE_DIR, "CRPS_all_structures_ratios.csv")
)

## Compact summary by hierarchy and level
crps_ratio_summary_by_level <- crps_full_ratio |>
  dplyr::group_by(Hierarchy, Level) |>
  dplyr::summarise(
    median_OLS_ratio  = median(OLS_ratio, na.rm = TRUE),
    median_WLS_ratio  = median(WLS_ratio, na.rm = TRUE),
    median_MinT_ratio = median(MinT_ratio, na.rm = TRUE),
    mean_OLS_ratio    = mean(OLS_ratio, na.rm = TRUE),
    mean_WLS_ratio    = mean(WLS_ratio, na.rm = TRUE),
    mean_MinT_ratio   = mean(MinT_ratio, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  crps_ratio_summary_by_level,
  file.path(TABLE_DIR, "CRPS_main_table_ratios_by_level.csv")
)


write_manuscript_table(
  crps_gT,
  crps_g1,
  crps_g2,
  "CRPS",
  file.path(TABLE_DIR, "interval_CRPS_table.tex")
)

saveRDS(
  list(
    rec_g1 = rec_g1,
    rec_g2 = rec_g2,
    rec_gT = rec_gT,
    CRPS_g1 = CRPS_g1,
    CRPS_g2 = CRPS_g2,
    CRPS_gT = CRPS_gT,
    crps_g1 = crps_g1,
    crps_g2 = crps_g2,
    crps_gT = crps_gT,
    crps_full = crps_full,
    names = nm,
    boot_info = list(
      ages = boot$ages,
      years.fit = boot$years.fit,
      years.for = boot$years.for,
      B = boot$B
    )
  ),
  file.path(RDS_DIR, "bootstrap_reconciled_crps.rds")
)

cat("\nDone.\n")
cat("Reconciled bootstrap paths and CRPS tables saved in outputs/.\n")