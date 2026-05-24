######################################################################
## R/helpers.R
## Consolidated helpers for point-forecast replication.
## This file replaces the previous helpers/structures/point_engine/
## tables_figures split. One file. Fewer doors for bugs to hide behind.
######################################################################



######################################################################
## From: helpers.R
######################################################################

######################################################################
## R/helpers.R
######################################################################

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

get_dx <- function(data, age_value, sex_value, cause_value, type_value) {
  out <- data |>
    dplyr::filter(
      Age == age_value,
      Sex == sex_value,
      Cause == cause_value,
      Type == type_value
    ) |>
    dplyr::arrange(as.numeric(Year)) |>
    dplyr::pull(dx)

  if (length(out) == 0) {
    stop(
      "No data found for Age=", age_value,
      ", Sex=", sex_value,
      ", Cause=", cause_value,
      ", Type=", type_value,
      call. = FALSE
    )
  }

  out
}

rmse_vec <- function(err) {
  apply(err, 2, function(x) sqrt(mean(x^2, na.rm = TRUE)))
}

mae_vec <- function(err) {
  apply(err, 2, function(x) mean(abs(x), na.rm = TRUE))
}

row_vec <- function(M) {
  as.vector(t(M))
}

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

latex_num <- function(x, is_best = FALSE) {
  out <- sprintf("%.2f", x)
  ifelse(is_best, paste0("\\textbf{", out, "}"), out)
}

latex_metric_row <- function(series, vals) {
  vals <- as.numeric(vals)
  best <- seq_along(vals) == which.min(vals)
  paste0(
    "& ", series, " & ",
    paste(latex_num(vals, best), collapse = " & "),
    " \\\\"  # LaTeX line break
  )
}



######################################################################
## From: structures.R
######################################################################

build_series_names <- function() {
  bottom_by_sex <- c(
    paste0("Cause ", 1:8, " - Males"),
    paste0("Cause ", 1:8, " - Females")
  )

  bottom_interlaced <- as.vector(rbind(
    paste0("Cause ", 1:8, " - Males"),
    paste0("Cause ", 1:8, " - Females")
  ))

  cause_totals <- paste0("Cause ", 1:8, " - Total")

  rows_g1 <- c("Total", "Total - Males", "Total - Females", bottom_by_sex)
  rows_g2 <- c("Total", cause_totals, bottom_interlaced)
  rows_gT <- c("Total", cause_totals, "Total - Males", "Total - Females", bottom_interlaced)

  list(
    bottom_by_sex = bottom_by_sex,
    bottom_interlaced = bottom_interlaced,
    cause_totals = cause_totals,
    rows_g1 = rows_g1,
    rows_g2 = rows_g2,
    rows_gT = rows_gT
  )
}

build_S_g1 <- function(rows_g1, bottom_by_sex) {
  Smat <- matrix(0, nrow = 19, ncol = 16)
  rownames(Smat) <- rows_g1
  colnames(Smat) <- bottom_by_sex

  Smat["Total", ] <- 1
  Smat["Total - Males", paste0("Cause ", 1:8, " - Males")] <- 1
  Smat["Total - Females", paste0("Cause ", 1:8, " - Females")] <- 1
  Smat[bottom_by_sex, bottom_by_sex] <- diag(16)

  Smat
}

build_S_g2 <- function(rows_g2, bottom_interlaced) {
  Smat <- matrix(0, nrow = 25, ncol = 16)
  rownames(Smat) <- rows_g2
  colnames(Smat) <- bottom_interlaced

  Smat["Total", ] <- 1

  for (i in 1:8) {
    Smat[paste0("Cause ", i, " - Total"),
         c(paste0("Cause ", i, " - Males"), paste0("Cause ", i, " - Females"))] <- 1
  }

  Smat[bottom_interlaced, bottom_interlaced] <- diag(16)
  Smat
}

build_S_gT <- function(rows_gT, bottom_interlaced) {
  Smat <- matrix(0, nrow = 27, ncol = 16)
  rownames(Smat) <- rows_gT
  colnames(Smat) <- bottom_interlaced

  Smat["Total", ] <- 1

  for (i in 1:8) {
    Smat[paste0("Cause ", i, " - Total"),
         c(paste0("Cause ", i, " - Males"), paste0("Cause ", i, " - Females"))] <- 1
  }

  Smat["Total - Males", paste0("Cause ", 1:8, " - Males")] <- 1
  Smat["Total - Females", paste0("Cause ", 1:8, " - Females")] <- 1
  Smat[bottom_interlaced, bottom_interlaced] <- diag(16)

  Smat
}



######################################################################
## From: point_engine.R
######################################################################

######################################################################
## R/point_engine.R
##
## Engine for point forecast reconciliation.
## Returns accuracy matrices plus actual and forecast matrices, because
## heatmaps and interval diagnostics need the raw forecast paths. Shocking.
######################################################################

build_matrices_for_age <- function(Output, age_value, cause_specific, cause_total, rows, structure) {
  oosl <- length(get_dx(Output, age_value, "M", cause_specific[1], "Forecasted"))
  isl  <- length(get_dx(Output, age_value, "M", cause_specific[1], "Fitted"))

  Yhat_unr <- matrix(NA_real_, nrow = oosl, ncol = length(rows), dimnames = list(NULL, rows))
  Yhat_is  <- matrix(NA_real_, nrow = isl,  ncol = length(rows), dimnames = list(NULL, rows))
  X        <- matrix(NA_real_, nrow = isl,  ncol = length(rows), dimnames = list(NULL, rows))
  Y        <- matrix(NA_real_, nrow = oosl, ncol = length(rows), dimnames = list(NULL, rows))

  fill_series <- function(M, type_value) {
    M[, "Total"] <- get_dx(Output, age_value, "T", cause_total, type_value)

    if ("Total - Males" %in% colnames(M)) {
      M[, "Total - Males"] <- get_dx(Output, age_value, "M", cause_total, type_value)
      M[, "Total - Females"] <- get_dx(Output, age_value, "F", cause_total, type_value)
    }

    for (i in 1:8) {
      cname_total <- paste0("Cause ", i, " - Total")
      cname_male  <- paste0("Cause ", i, " - Males")
      cname_fem   <- paste0("Cause ", i, " - Females")

      if (cname_total %in% colnames(M)) {
        M[, cname_total] <- get_dx(Output, age_value, "T", cause_specific[i], type_value)
      }
      if (cname_male %in% colnames(M)) {
        M[, cname_male] <- get_dx(Output, age_value, "M", cause_specific[i], type_value)
      }
      if (cname_fem %in% colnames(M)) {
        M[, cname_fem] <- get_dx(Output, age_value, "F", cause_specific[i], type_value)
      }
    }

    M
  }

  Yhat_unr <- fill_series(Yhat_unr, "Forecasted")
  Yhat_is  <- fill_series(Yhat_is,  "Fitted")
  X        <- fill_series(X,        "Observed_In_sample")
  Y        <- fill_series(Y,        "Observed_Out_of_sample")

  list(Yhat_unr = Yhat_unr, Yhat_is = Yhat_is, X = X, Y = Y)
}

run_point_structure <- function(Output, ages, cause_specific, cause_total, rows, Smat, agg_rows, structure_name) {
  n_series <- length(rows)
  n_ages <- length(ages)

  RMSE_base <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  RMSE_OLS  <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  RMSE_WLS  <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  RMSE_MinT <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))

  MAE_base <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  MAE_OLS  <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  MAE_WLS  <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))
  MAE_MinT <- matrix(NA_real_, nrow = n_series, ncol = n_ages, dimnames = list(rows, ages))

  Yhat_unr_list     <- vector("list", n_ages)
  Yhat_recOLS_list  <- vector("list", n_ages)
  Yhat_recWLS_list  <- vector("list", n_ages)
  Yhat_recMinT_list <- vector("list", n_ages)
  Y_list            <- vector("list", n_ages)
  Resid_list        <- vector("list", n_ages)

  top_actual <- NULL
  top_coda <- NULL
  top_mint <- NULL

  for (eta in seq_along(ages)) {
    age_value <- ages[eta]
    mats <- build_matrices_for_age(Output, age_value, cause_specific, cause_total, rows, structure_name)
    Resid <- mats$Yhat_is - mats$X

    Yhat_OLS <- FoReco::csrec(
      base = mats$Yhat_unr,
      agg_mat = Smat[agg_rows, , drop = FALSE],
      comb = "ols",
      nn = "strc_osqp",
      res = Resid
    )
    Yhat_WLS <- FoReco::csrec(
      base = mats$Yhat_unr,
      agg_mat = Smat[agg_rows, , drop = FALSE],
      comb = "wls",
      nn = "strc_osqp",
      res = Resid
    )
    Yhat_MinT <- FoReco::csrec(
      base = mats$Yhat_unr,
      agg_mat = Smat[agg_rows, , drop = FALSE],
      comb = "shr",
      nn = "strc_osqp",
      res = Resid
    )

    colnames(Yhat_OLS) <- rows
    colnames(Yhat_WLS) <- rows
    colnames(Yhat_MinT) <- rows

    err_base <- mats$Y - mats$Yhat_unr
    err_OLS  <- mats$Y - Yhat_OLS
    err_WLS  <- mats$Y - Yhat_WLS
    err_MinT <- mats$Y - Yhat_MinT

    RMSE_base[, eta] <- rmse_vec(err_base)
    RMSE_OLS[, eta]  <- rmse_vec(err_OLS)
    RMSE_WLS[, eta]  <- rmse_vec(err_WLS)
    RMSE_MinT[, eta] <- rmse_vec(err_MinT)

    MAE_base[, eta] <- mae_vec(err_base)
    MAE_OLS[, eta]  <- mae_vec(err_OLS)
    MAE_WLS[, eta]  <- mae_vec(err_WLS)
    MAE_MinT[, eta] <- mae_vec(err_MinT)

    Yhat_unr_list[[eta]]     <- mats$Yhat_unr
    Yhat_recOLS_list[[eta]]  <- Yhat_OLS
    Yhat_recWLS_list[[eta]]  <- Yhat_WLS
    Yhat_recMinT_list[[eta]] <- Yhat_MinT
    Y_list[[eta]]            <- mats$Y
    Resid_list[[eta]]        <- Resid

    if (eta == 1) {
      top_actual <- matrix(NA_real_, nrow = nrow(mats$Y), ncol = length(ages))
      top_coda <- matrix(NA_real_, nrow = nrow(mats$Y), ncol = length(ages))
      top_mint <- matrix(NA_real_, nrow = nrow(mats$Y), ncol = length(ages))
    }
    top_actual[, eta] <- mats$Y[, "Total"]
    top_coda[, eta] <- mats$Yhat_unr[, "Total"]
    top_mint[, eta] <- Yhat_MinT[, "Total"]

    message(structure_name, " completed for age ", age_value)
  }

  list(
    RMSE = list(CoDA = RMSE_base, OLS = RMSE_OLS, WLS = RMSE_WLS, MinT = RMSE_MinT),
    MAE = list(CoDA = MAE_base, OLS = MAE_OLS, WLS = MAE_WLS, MinT = MAE_MinT),
    Forecasts = list(CoDA = Yhat_unr_list, OLS = Yhat_recOLS_list, WLS = Yhat_recWLS_list, MinT = Yhat_recMinT_list),
    Actual = Y_list,
    Residuals = Resid_list,
    top = list(actual = top_actual, coda = top_coda, mint = top_mint),
    label = structure_name
  )
}



######################################################################
## From: tables_figures.R
######################################################################

######################################################################
## R/tables_figures.R
##
## Helpers for point-forecast tables and official figures.
######################################################################

build_point_tables <- function(res_gT, res_g1, res_g2, names, metric = c("RMSE", "MAE")) {
  metric <- match.arg(metric)

  M_gT <- res_gT[[metric]]
  M_g1 <- res_g1[[metric]]
  M_g2 <- res_g2[[metric]]

  gT <- make_metric_matrix(M_gT$CoDA, M_gT$OLS, M_gT$WLS, M_gT$MinT, names$rows_gT)
  g1_display <- c("Total", "Total - Males", "Total - Females", names$bottom_interlaced)
  g1 <- make_metric_matrix(M_g1$CoDA, M_g1$OLS, M_g1$WLS, M_g1$MinT, g1_display)
  g2 <- make_metric_matrix(M_g2$CoDA, M_g2$OLS, M_g2$WLS, M_g2$MinT, names$rows_g2)

  full <- dplyr::bind_rows(
    gT |>
      dplyr::mutate(
        Hierarchy = "Full grouped structure",
        Level = c("Top", rep("Middle (cause)", 8), rep("Middle (sex)", 2), rep("Bottom", 16)),
        .before = Series
      ),
    g1 |>
      dplyr::mutate(
        Hierarchy = "Sex-based hierarchy",
        Level = c("Top", rep("Middle (sex)", 2), rep("Bottom", 16)),
        .before = Series
      ),
    g2 |>
      dplyr::mutate(
        Hierarchy = "Cause-based hierarchy",
        Level = c("Top", rep("Middle (cause)", 8), rep("Bottom", 16)),
        .before = Series
      )
  )

  list(gT = gT, g1 = g1, g2 = g2, full = full)
}

write_manuscript_table <- function(gT, g1, g2, metric = c("RMSE", "MAE", "CRPS"), file) {
  metric <- match.arg(metric)
  numcols <- c("CoDA", "OLS", "WLS", "MinT")
  br <- "\\\\"

  caption <- paste0(
    "\\small Median forecast accuracy results: ", metric,
    " computed for the three levels of the grouped and hierarchical structures. The best method is highlighted in bold."
  )
  label <- paste0("tab:", tolower(metric), "_median")

  lines <- c(
    "\\begin{center}",
    "\\tabcolsep 0.315in",
    "\\renewcommand\\arraystretch{0.98}",
    "\\begin{small}",
    "\\begin{longtable}{@{}llrrrr@{}}",
    paste0("\\caption{", caption, "}\\label{", label, "} ", br),
    "\\toprule",
    paste0("Hierarchy & Series & CoDA & OLS & WLS & MinT ", br),
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    paste0("Hierarchy & Series & CoDA & OLS & WLS & MinT ", br),
    "\\midrule",
    "\\endhead",
    paste0("\\hline \\multicolumn{6}{r}{{Continued on next page}} ", br),
    "\\endfoot",
    "\\endlastfoot",
    paste0("Grouped & \\multicolumn{5}{l}{\\underline{Top-level:}}", br),
    latex_metric_row("Total", gT[1, numcols]),
    "\\cmidrule{2-6}",
    paste0("& \\multicolumn{5}{l}{\\underline{Middle-level:}}", br)
  )

  for (ii in 2:9) lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))
  for (ii in 10:11) lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))

  lines <- c(lines, "\\cmidrule{2-6}", paste0("& \\multicolumn{5}{l}{\\underline{Bottom-level:}}", br))
  for (ii in 12:27) lines <- c(lines, latex_metric_row(gT$Series[ii], gT[ii, numcols]))

  lines <- c(
    lines,
    "\\midrule",
    paste0("Sex & \\multicolumn{5}{l}{\\underline{Top-level:}}", br),
    paste0("$\\downarrow$ ", substring(latex_metric_row("Total", g1[1, numcols]), 2)),
    "\\cmidrule{2-6}",
    paste0("Cause & \\multicolumn{5}{l}{\\underline{Middle-level:}}", br),
    latex_metric_row(g1$Series[2], g1[2, numcols]),
    latex_metric_row(g1$Series[3], g1[3, numcols]),
    "\\cmidrule{2-6}",
    paste0("& \\multicolumn{5}{l}{\\underline{Bottom-level:}}", br)
  )
  for (ii in 4:19) lines <- c(lines, latex_metric_row(g1$Series[ii], g1[ii, numcols]))

  lines <- c(
    lines,
    "\\midrule",
    paste0("Cause & \\multicolumn{5}{l}{\\underline{Top-level:}}", br),
    paste0("$\\downarrow$ ", substring(latex_metric_row("Total", g2[1, numcols]), 2)),
    "\\cmidrule{2-6}",
    paste0("Sex & \\multicolumn{5}{l}{\\underline{Middle-level:}}", br)
  )
  for (ii in 2:9) lines <- c(lines, latex_metric_row(g2$Series[ii], g2[ii, numcols]))
  lines <- c(lines, "\\cmidrule{2-6}", paste0("& \\multicolumn{5}{l}{\\underline{Bottom-level:}}", br))
  for (ii in 10:25) lines <- c(lines, latex_metric_row(g2$Series[ii], g2[ii, numcols]))

  lines <- c(lines, "\\bottomrule", "\\end{longtable}", "\\end{small}", "\\end{center}")
  writeLines(lines, file)
}

save_plot <- function(plot, filename, figure_dir, width = 8.5, height = 5.5) {

  ## Minimal-output replication packages do not need figure files.
  ## Set CREATE_FIGURES <- TRUE in run_all.R to restore plot exports.
  if (exists("CREATE_FIGURES", inherits = TRUE) &&
      identical(get("CREATE_FIGURES", inherits = TRUE), FALSE)) {
    return(invisible(plot))
  }

  if (!dir.exists(figure_dir)) {
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
}

matrix_to_heat_df <- function(mat, method_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("Year", "Age", "RSE")
  df$Year <- as.numeric(as.character(df$Year))
  df$Age <- as.numeric(as.character(df$Age))
  df$Method <- method_name
  df
}

build_heatmap_top <- function(res_gT, res_g1, res_g2, ages, table_dir, figure_dir) {
  
  top_actual  <- res_gT$top$actual
  top_coda    <- res_gT$top$coda
  top_mint_gT <- res_gT$top$mint
  top_mint_g1 <- res_g1$top$mint
  top_mint_g2 <- res_g2$top$mint
  
  stopifnot(
    is.matrix(top_actual),
    identical(dim(top_actual), dim(top_coda)),
    identical(dim(top_actual), dim(top_mint_gT)),
    identical(dim(top_actual), dim(top_mint_g1)),
    identical(dim(top_actual), dim(top_mint_g2)),
    ncol(top_actual) == length(ages)
  )
  
  H <- nrow(top_actual)
  forecast_years <- 2001:(2000 + H)
  
  rownames(top_actual)  <- forecast_years
  rownames(top_coda)    <- forecast_years
  rownames(top_mint_gT) <- forecast_years
  rownames(top_mint_g1) <- forecast_years
  rownames(top_mint_g2) <- forecast_years
  
  colnames(top_actual)  <- ages
  colnames(top_coda)    <- ages
  colnames(top_mint_gT) <- ages
  colnames(top_mint_g1) <- ages
  colnames(top_mint_g2) <- ages
  
  se_coda    <- (top_actual - top_coda)^2
  se_grouped <- (top_actual - top_mint_gT)^2
  se_g1      <- (top_actual - top_mint_g1)^2
  se_g2      <- (top_actual - top_mint_g2)^2
  
  eps <- 1e-8
  
  gain_grouped <- se_coda / (se_grouped + eps)
  gain_g1      <- se_coda / (se_g1 + eps)
  gain_g2      <- se_coda / (se_g2 + eps)
  
  heat_df <- dplyr::bind_rows(
    matrix_to_heat_df(gain_grouped, "Grouped"),
    matrix_to_heat_df(gain_g1,      "Sex to Cause"),
    matrix_to_heat_df(gain_g2,      "Cause to Sex")
  ) |>
    dplyr::mutate(
      Method = factor(
        Method,
        levels = c("Grouped", "Sex to Cause", "Cause to Sex")
      ),
      RSE_cap = pmin(RSE, 5)
    )
  
  if (!dir.exists(table_dir)) {
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  readr::write_csv(
    heat_df,
    file.path(table_dir, "diagnostic_heatmap_top_CoDA_over_MinT_gain.csv")
  )
  
  p_heat <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = Year, y = Age, fill = RSE_cap)
  ) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::facet_wrap(~Method, nrow = 1) +
    ggplot2::scale_fill_gradientn(
      colours = c("white", "mistyrose", "red", "darkred"),
      limits = c(0, 5),
      breaks = c(0, 1, 2, 3, 4, 5),
      name = "RSE(CoDA) /\nRSE(MinT)"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(min(forecast_years), max(forecast_years), by = 3),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(min(as.numeric(ages)), max(as.numeric(ages)), by = 10),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Year", y = "Age") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey85", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.35, "lines"),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.25, "cm")
    )
  
  save_plot(p_heat, "Heatmap-Top", figure_dir, 9, 3.8)
  
  p_heat
}