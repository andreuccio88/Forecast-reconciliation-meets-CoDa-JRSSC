######################################################################
## residual_block_bootstrap_crps_legacy.R
##
## Distribution-free probabilistic reconciliation via residual block
## bootstrap for the legacy JRSSC folder structure:
##
##   project_root/
##   ├── run_all.R
##   ├── data/ or Data/
##   └── R/
##       ├── base_fcast.R
##       ├── reco_fcast.R
##       ├── gen_boot.R
##       ├── reco_boot.R
##       ├── helpers.R
##       └── residual_block_bootstrap_crps_legacy.R   <-- put this file here
##
## Run after the point-forecast step has created:
##   outputs/rds/point_forecasts.rds
##
## The script creates residual block-bootstrap predictive draws from:
##   base forecast + sampled in-sample residual block
## and then reconciles each draw using OLS/WLS/MinT maps.
##
## Main metric: CRPS and CRPS ratio relative to unreconciled CoDA,
## matching the logic used for the existing probabilistic/CRPS tables:
##   horizon average -> median over ages.
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

R_DIR      <- file.path(PROJECT_ROOT, "R")
OUT_DIR    <- file.path(PROJECT_ROOT, "outputs")
RDS_DIR    <- file.path(OUT_DIR, "rds")
TABLE_DIR  <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")

for (dd in c(OUT_DIR, RDS_DIR, TABLE_DIR, FIGURE_DIR)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

######################################################################
## User options
######################################################################

## Final paper: 10000. For testing: set BB_NSIM <- 500 or 1000 before source().
BB_NSIM <- if (exists("BB_NSIM", inherits = TRUE)) get("BB_NSIM", inherits = TRUE) else 10000L

## Block length. With H = 15, this samples complete 15-year residual paths.
## For robustness: rerun with 3, 5, 10, 15.
BB_BLOCK_SIZE <- if (exists("BB_BLOCK_SIZE", inherits = TRUE)) get("BB_BLOCK_SIZE", inherits = TRUE) else 15L

## Use all draws for CRPS by default. For faster testing: BB_CRPS_EVAL_DRAWS <- 1000.
BB_CRPS_EVAL_DRAWS <- if (exists("BB_CRPS_EVAL_DRAWS", inherits = TRUE)) get("BB_CRPS_EVAL_DRAWS", inherits = TRUE) else NULL

BB_CENTER_RESIDUALS <- if (exists("BB_CENTER_RESIDUALS", inherits = TRUE)) get("BB_CENTER_RESIDUALS", inherits = TRUE) else TRUE
BB_TRUNCATE_BASE_AT_ZERO <- if (exists("BB_TRUNCATE_BASE_AT_ZERO", inherits = TRUE)) get("BB_TRUNCATE_BASE_AT_ZERO", inherits = TRUE) else TRUE

## After linear reconciliation, truncate bottom-level draws at zero and rebuild all upper levels.
## This preserves coherence and non-negativity, but it is a post-processing step.
BB_NONNEGATIVE_RECONCILED_POSTPROCESS <- if (exists("BB_NONNEGATIVE_RECONCILED_POSTPROCESS", inherits = TRUE)) {
  get("BB_NONNEGATIVE_RECONCILED_POSTPROCESS", inherits = TRUE)
} else {
  TRUE
}

## TRUE tries to infer the exact linear reconciliation maps from FoReco::csrec using identity base forecasts.
## If it fails, the script falls back to explicit OLS/WLS/MinT formulas.
BB_USE_FORECO_MAPPING <- if (exists("BB_USE_FORECO_MAPPING", inherits = TRUE)) get("BB_USE_FORECO_MAPPING", inherits = TRUE) else TRUE

## Used only by the formula fallback for MinT if FoReco mapping is unavailable.
BB_MINT_SHRINK_LAMBDA <- if (exists("BB_MINT_SHRINK_LAMBDA", inherits = TRUE)) get("BB_MINT_SHRINK_LAMBDA", inherits = TRUE) else 0.20

BB_SEED <- if (exists("BB_SEED", inherits = TRUE)) get("BB_SEED", inherits = TRUE) else 20260506L
BB_PREFIX <- if (exists("BB_PREFIX", inherits = TRUE)) get("BB_PREFIX", inherits = TRUE) else "BB_"

######################################################################
## Packages and helpers
######################################################################

required_packages <- c("dplyr", "tidyr", "readr", "ggplot2", "FoReco", "MASS")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

if (!file.exists(file.path(R_DIR, "helpers.R"))) {
  stop("Cannot find R/helpers.R. Put this script in the legacy project root/R folder and run from the project root.", call. = FALSE)
}
source(file.path(R_DIR, "helpers.R"))

POINT_FILE <- file.path(RDS_DIR, "point_forecasts.rds")
if (!file.exists(POINT_FILE)) {
  stop(
    "Cannot find outputs/rds/point_forecasts.rds.\n",
    "Run the point-forecast steps first, e.g. source('run_all.R') through step 2, or run R/base_fcast.R and R/reco_fcast.R.",
    call. = FALSE
  )
}

pf <- readRDS(POINT_FILE)

nm <- pf$names
ages <- as.numeric(as.character(pf$ages))
FORECAST_H <- nrow(pf$res_g1$Forecasts$CoDA[[1]])
BASE_END_YEAR <- 2000L

S_g1 <- build_S_g1(nm$rows_g1, nm$bottom_by_sex)
S_g2 <- build_S_g2(nm$rows_g2, nm$bottom_interlaced)
S_gT <- build_S_gT(nm$rows_gT, nm$bottom_interlaced)

structures <- list(
  g1 = list(
    label = "Sex-based hierarchy",
    res = pf$res_g1,
    S = S_g1,
    agg_rows = 1:3,
    row_order = nm$rows_g1,
    display_order = c("Total", "Total - Males", "Total - Females", nm$bottom_interlaced)
  ),
  g2 = list(
    label = "Cause-based hierarchy",
    res = pf$res_g2,
    S = S_g2,
    agg_rows = 1:9,
    row_order = nm$rows_g2,
    display_order = nm$rows_g2
  ),
  gT = list(
    label = "Full grouped structure",
    res = pf$res_gT,
    S = S_gT,
    agg_rows = 1:11,
    row_order = nm$rows_gT,
    display_order = nm$rows_gT
  )
)

method_levels <- c("base", "OLS", "WLS", "MinT")
method_labels <- c(base = "CoDA", OLS = "OLS", WLS = "WLS", MinT = "MinT")

######################################################################
## Utility functions
######################################################################

safe_div <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

write_csv_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
}

save_rds_safe <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
}

safe_solve <- function(A, B = NULL) {
  if (is.null(B)) {
    out <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  } else {
    out <- tryCatch(solve(A, B), error = function(e) MASS::ginv(A) %*% B)
  }
  out
}

sample_crps <- function(draws, y) {
  ## Empirical CRPS from predictive draws, computed in O(L log L):
  ## mean |X-y| - 1/(2L^2) sum_i sum_j |X_i-X_j|.
  x <- as.numeric(draws)
  x <- x[is.finite(x)]
  y <- as.numeric(y)
  if (!is.finite(y) || length(x) == 0) return(NA_real_)

  x <- sort(x)
  n <- length(x)
  term1 <- mean(abs(x - y))
  coeff <- 2 * seq_len(n) - n - 1
  pair_half <- sum(coeff * x) / (n^2)
  term1 - pair_half
}

interval_score <- function(y, l, u, alpha) {
  (u - l) + (2 / alpha) * (l - y) * as.integer(y < l) +
    (2 / alpha) * (y - u) * as.integer(y > u)
}

summarise_interval_draws <- function(draws, y, alpha) {
  l <- stats::quantile(draws, probs = alpha / 2, na.rm = TRUE, names = FALSE, type = 8)
  u <- stats::quantile(draws, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE, type = 8)
  tibble::tibble(
    alpha = alpha,
    nominal_coverage = 1 - alpha,
    lower = l,
    upper = u,
    covered = as.integer(y >= l & y <= u),
    width = u - l,
    interval_score = interval_score(y, l, u, alpha)
  )
}

sample_residual_path <- function(resid_mat, H, block_size) {
  ## resid_mat: T x n matrix. Output: H x n sampled consecutive residual path.
  resid_mat <- as.matrix(resid_mat)
  Tn <- nrow(resid_mat)
  n <- ncol(resid_mat)

  if (block_size < 1) stop("block_size must be >= 1", call. = FALSE)
  if (block_size > Tn) stop("block_size cannot exceed the number of in-sample residual rows.", call. = FALSE)

  max_start <- Tn - block_size + 1
  out <- matrix(NA_real_, nrow = H, ncol = n)
  colnames(out) <- colnames(resid_mat)

  filled <- 0
  while (filled < H) {
    start <- sample.int(max_start, size = 1)
    block <- resid_mat[start:(start + block_size - 1), , drop = FALSE]
    take <- min(nrow(block), H - filled)
    out[(filled + 1):(filled + take), ] <- block[seq_len(take), , drop = FALSE]
    filled <- filled + take
  }
  out
}

classify_level <- function(series, structure) {
  if (series == "Total") return("Top")

  is_bottom <- grepl("^Cause [0-9]+ - (Males|Females)$", series)
  if (is_bottom) return("Bottom")

  if (grepl("^Cause [0-9]+ - Total$", series)) return("Middle (cause)")
  if (series %in% c("Total - Males", "Total - Females")) return("Middle (sex)")

  "Other"
}

nonnegative_postprocess_draw <- function(draw_mat, S) {
  ## draw_mat: H x n. Take bottom-level rows, truncate, rebuild all nodes via S.
  n <- ncol(draw_mat)
  nb <- ncol(S)
  bot_idx <- (n - nb + 1):n
  bottom <- pmax(draw_mat[, bot_idx, drop = FALSE], 0)
  rebuilt <- t(S %*% t(bottom))
  colnames(rebuilt) <- colnames(draw_mat)
  rebuilt
}

make_M_matrices_formula <- function(S, Resid_for_G) {
  ## Explicit linear reconciliation maps M = S G.
  ## OLS/WLS/MinT are formula-based fallbacks if FoReco mapping fails.
  n <- nrow(S)
  nb <- ncol(S)

  G_OLS <- safe_solve(t(S) %*% S, t(S))
  M_OLS <- S %*% G_OLS

  v <- colMeans(Resid_for_G^2, na.rm = TRUE)
  v[!is.finite(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  min_pos <- min(v[v > 0], na.rm = TRUE)
  if (!is.finite(min_pos)) min_pos <- 1
  v <- pmax(v, min_pos * 1e-6)
  w <- 1 / v
  Sw <- sweep(S, 1, w, "*")
  right_wls <- sweep(t(S), 2, w, "*")
  G_WLS <- safe_solve(t(S) %*% Sw, right_wls)
  M_WLS <- S %*% G_WLS

  W <- stats::cov(Resid_for_G, use = "pairwise.complete.obs")
  if (!is.matrix(W) || any(!is.finite(W))) {
    W <- diag(v, n)
  }
  W <- as.matrix(W)
  W_diag <- diag(diag(W), nrow = n, ncol = n)
  lambda <- as.numeric(BB_MINT_SHRINK_LAMBDA)
  lambda <- min(max(lambda, 0), 1)
  W_shr <- (1 - lambda) * W + lambda * W_diag

  ridge <- max(1e-8, mean(diag(W_shr), na.rm = TRUE) * 1e-8)
  W_shr <- W_shr + diag(ridge, n)
  W_inv <- safe_solve(W_shr)
  G_MinT <- safe_solve(t(S) %*% W_inv %*% S, t(S) %*% W_inv)
  M_MinT <- S %*% G_MinT

  list(OLS = M_OLS, WLS = M_WLS, MinT = M_MinT, source = "formula_fallback")
}

make_M_matrices <- function(S, agg_rows, Resid_for_G) {
  n <- nrow(S)
  rows <- rownames(S)

  if (isTRUE(BB_USE_FORECO_MAPPING)) {
    out <- tryCatch({
      Ibase <- diag(n)
      colnames(Ibase) <- rows
      rownames(Ibase) <- paste0("basis_", seq_len(n))

      combs <- c(OLS = "ols", WLS = "wls", MinT = "shr")
      M_list <- list()

      for (mth in names(combs)) {
        recI <- FoReco::csrec(
          base = Ibase,
          agg_mat = S[agg_rows, , drop = FALSE],
          comb = combs[[mth]],
          res = Resid_for_G
        )
        M <- t(as.matrix(recI))
        rownames(M) <- rows
        colnames(M) <- rows
        M_list[[mth]] <- M
      }

      M_list$source <- "FoReco_identity_map"
      M_list
    }, error = function(e) {
      message("FoReco identity-map construction failed; using formula fallback. Reason: ", e$message)
      NULL
    })

    if (!is.null(out)) return(out)
  }

  make_M_matrices_formula(S, Resid_for_G)
}

make_residual_block_bootstrap_draws <- function(base_mean,
                                                resid_mat,
                                                S,
                                                M_list,
                                                nsim,
                                                block_size,
                                                center_residuals,
                                                truncate_base_at_zero,
                                                nonnegative_reconciled,
                                                seed) {
  ## base_mean: H x n matrix of base CoDA forecasts.
  ## resid_mat: T x n matrix of in-sample residuals, observed - fitted.
  ## M_list: OLS/WLS/MinT n x n reconciliation maps.

  H <- nrow(base_mean)
  n <- ncol(base_mean)

  if (ncol(resid_mat) != n) stop("resid_mat and base_mean have incompatible number of series.", call. = FALSE)
  if (nrow(S) != n) stop("S and base_mean have incompatible number of series.", call. = FALSE)
  if (block_size > nrow(resid_mat)) stop("block_size > number of residual rows.", call. = FALSE)

  if (isTRUE(center_residuals)) {
    resid_mat <- sweep(resid_mat, 2, colMeans(resid_mat, na.rm = TRUE), "-")
  }

  set.seed(seed)

  out <- list(
    base = array(NA_real_, dim = c(H, n, nsim), dimnames = list(NULL, colnames(base_mean), NULL)),
    OLS  = array(NA_real_, dim = c(H, n, nsim), dimnames = list(NULL, colnames(base_mean), NULL)),
    WLS  = array(NA_real_, dim = c(H, n, nsim), dimnames = list(NULL, colnames(base_mean), NULL)),
    MinT = array(NA_real_, dim = c(H, n, nsim), dimnames = list(NULL, colnames(base_mean), NULL))
  )

  for (b in seq_len(nsim)) {
    E_path <- sample_residual_path(resid_mat, H = H, block_size = block_size)
    Yb <- base_mean + E_path

    if (isTRUE(truncate_base_at_zero)) {
      Yb <- pmax(Yb, 0)
    }
    out$base[, , b] <- Yb

    for (mth in c("OLS", "WLS", "MinT")) {
      rec <- t(M_list[[mth]] %*% t(Yb))
      colnames(rec) <- colnames(Yb)

      if (isTRUE(nonnegative_reconciled)) {
        rec <- nonnegative_postprocess_draw(rec, S)
      }
      out[[mth]][, , b] <- rec
    }
  }

  out
}

make_crps_matrix <- function(crps_by_series_age, structure, method, row_order, ages) {
  M <- matrix(NA_real_, nrow = length(row_order), ncol = length(ages), dimnames = list(row_order, as.character(ages)))
  tmp <- crps_by_series_age |>
    dplyr::filter(.data$Structure == structure, .data$Method == method)

  for (ii in seq_len(nrow(tmp))) {
    s <- tmp$Series[ii]
    a <- as.character(tmp$Age[ii])
    if (s %in% rownames(M) && a %in% colnames(M)) {
      M[s, a] <- tmp$CRPS[ii]
    }
  }
  M
}

######################################################################
## Select draw indices for CRPS
######################################################################

if (is.null(BB_CRPS_EVAL_DRAWS)) {
  crps_draw_idx <- seq_len(BB_NSIM)
} else {
  set.seed(BB_SEED + 1000L)
  crps_draw_idx <- sort(sample.int(BB_NSIM, min(BB_NSIM, BB_CRPS_EVAL_DRAWS)))
}

cat("\n============================================================\n")
cat("Residual block-bootstrap probabilistic reconciliation\n")
cat("============================================================\n")
cat("Project root: ", PROJECT_ROOT, "\n", sep = "")
cat("nsim: ", BB_NSIM, "\n", sep = "")
cat("block_size: ", BB_BLOCK_SIZE, "\n", sep = "")
cat("center_residuals: ", BB_CENTER_RESIDUALS, "\n", sep = "")
cat("truncate_base_at_zero: ", BB_TRUNCATE_BASE_AT_ZERO, "\n", sep = "")
cat("nonnegative_reconciled_postprocess: ", BB_NONNEGATIVE_RECONCILED_POSTPROCESS, "\n", sep = "")
cat("CRPS draws used: ", length(crps_draw_idx), "\n\n", sep = "")

######################################################################
## Main loop
######################################################################

crps_long_list <- list()
pi_long_list <- list()
diag_list <- list()
crps_idx <- 1L
pi_idx <- 1L
diag_idx <- 1L

for (str_name in names(structures)) {
  info <- structures[[str_name]]
  cat("Structure: ", str_name, " - ", info$label, "\n", sep = "")

  for (a_idx in seq_along(ages)) {
    age <- ages[a_idx]
    cat("  Age ", a_idx, "/", length(ages), " (", age, ")\n", sep = "")

    base_mean <- as.matrix(info$res$Forecasts$CoDA[[a_idx]])
    Y <- as.matrix(info$res$Actual[[a_idx]])

    ## In helpers.R, residuals were stored as fitted - observed.
    ## For bootstrap draws we need observed - fitted, because we add residuals to the forecast.
    Resid_for_G <- as.matrix(info$res$Residuals[[a_idx]])
    Resid_for_boot <- -Resid_for_G

    M_list <- make_M_matrices(
      S = info$S,
      agg_rows = info$agg_rows,
      Resid_for_G = Resid_for_G
    )

    draws <- make_residual_block_bootstrap_draws(
      base_mean = base_mean,
      resid_mat = Resid_for_boot,
      S = info$S,
      M_list = M_list,
      nsim = BB_NSIM,
      block_size = BB_BLOCK_SIZE,
      center_residuals = BB_CENTER_RESIDUALS,
      truncate_base_at_zero = BB_TRUNCATE_BASE_AT_ZERO,
      nonnegative_reconciled = BB_NONNEGATIVE_RECONCILED_POSTPROCESS,
      seed = BB_SEED + 1000L * match(str_name, names(structures)) + a_idx
    )

    for (mth in method_levels) {
      arr <- draws[[mth]]
      series <- colnames(base_mean)

      diag_list[[diag_idx]] <- tibble::tibble(
        Engine = "Residual block bootstrap",
        Structure = str_name,
        Structure_label = info$label,
        Method = mth,
        Method_label = method_labels[[mth]],
        Age = age,
        block_size = BB_BLOCK_SIZE,
        nsim = BB_NSIM,
        mapping_source = ifelse(is.null(M_list$source), NA_character_, M_list$source),
        share_negative_draws = mean(arr < 0, na.rm = TRUE),
        min_draw = suppressWarnings(min(arr, na.rm = TRUE)),
        max_draw = suppressWarnings(max(arr, na.rm = TRUE))
      )
      diag_idx <- diag_idx + 1L

      for (s in seq_along(series)) {
        level_s <- classify_level(series[s], str_name)

        for (h in seq_len(FORECAST_H)) {
          y <- Y[h, s]
          d <- arr[h, s, crps_draw_idx]

          crps_long_list[[crps_idx]] <- tibble::tibble(
            Engine = "Residual block bootstrap",
            Structure = str_name,
            Structure_label = info$label,
            Method = mth,
            Method_label = method_labels[[mth]],
            Level = level_s,
            Series = series[s],
            Age = age,
            Horizon = h,
            Year = BASE_END_YEAR + h,
            Y = y,
            CRPS = sample_crps(d, y),
            block_size = BB_BLOCK_SIZE,
            nsim = BB_NSIM
          )
          crps_idx <- crps_idx + 1L

          pi_tmp <- dplyr::bind_rows(
            summarise_interval_draws(d, y, alpha = 0.20),
            summarise_interval_draws(d, y, alpha = 0.05)
          ) |>
            dplyr::mutate(
              Engine = "Residual block bootstrap",
              Structure = str_name,
              Structure_label = info$label,
              Method = mth,
              Method_label = method_labels[[mth]],
              Level = level_s,
              Series = series[s],
              Age = age,
              Horizon = h,
              Year = BASE_END_YEAR + h,
              Y = y,
              block_size = BB_BLOCK_SIZE,
              nsim = BB_NSIM,
              .before = 1
            )

          pi_long_list[[pi_idx]] <- pi_tmp
          pi_idx <- pi_idx + 1L
        }
      }
    }

    rm(draws)
    gc(verbose = FALSE)
  }
}

crps_long <- dplyr::bind_rows(crps_long_list) |>
  dplyr::mutate(
    Method = factor(.data$Method, levels = method_levels),
    Method_label = factor(.data$Method_label, levels = c("CoDA", "OLS", "WLS", "MinT")),
    Structure_label = factor(
      .data$Structure_label,
      levels = c("Sex-based hierarchy", "Cause-based hierarchy", "Full grouped structure")
    )
  )

pi_long <- dplyr::bind_rows(pi_long_list) |>
  dplyr::mutate(
    Method = factor(.data$Method, levels = method_levels),
    Method_label = factor(.data$Method_label, levels = c("CoDA", "OLS", "WLS", "MinT")),
    Structure_label = factor(
      .data$Structure_label,
      levels = c("Sex-based hierarchy", "Cause-based hierarchy", "Full grouped structure")
    )
  )

negative_diagnostics <- dplyr::bind_rows(diag_list) |>
  dplyr::mutate(
    Method = factor(.data$Method, levels = method_levels),
    Method_label = factor(.data$Method_label, levels = c("CoDA", "OLS", "WLS", "MinT")),
    Structure_label = factor(
      .data$Structure_label,
      levels = c("Sex-based hierarchy", "Cause-based hierarchy", "Full grouped structure")
    )
  )

save_rds_safe(crps_long, file.path(RDS_DIR, paste0(BB_PREFIX, "CRPS_long_horizon.rds")))
write_csv_safe(crps_long, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_long_horizon.csv")))
save_rds_safe(pi_long, file.path(RDS_DIR, paste0(BB_PREFIX, "PI_long_horizon.rds")))
write_csv_safe(pi_long, file.path(TABLE_DIR, paste0(BB_PREFIX, "PI_long_horizon.csv")))
write_csv_safe(negative_diagnostics, file.path(TABLE_DIR, paste0(BB_PREFIX, "negative_diagnostics.csv")))

######################################################################
## CRPS aggregation: mean across horizons, median across ages
######################################################################

crps_by_series_age <- crps_long |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$Level, .data$Series, .data$Age,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    n_horizons = dplyr::n(),
    CRPS = mean(.data$CRPS, na.rm = TRUE),
    .groups = "drop"
  )

save_rds_safe(crps_by_series_age, file.path(RDS_DIR, paste0(BB_PREFIX, "CRPS_by_series_age.rds")))
write_csv_safe(crps_by_series_age, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_by_series_age.csv")))

base_crps_by_age <- crps_by_series_age |>
  dplyr::filter(.data$Method == "base") |>
  dplyr::select(
    .data$Structure, .data$Level, .data$Series, .data$Age,
    base_CRPS = .data$CRPS
  )

crps_ratio_by_series_age <- crps_by_series_age |>
  dplyr::left_join(
    base_crps_by_age,
    by = c("Structure", "Level", "Series", "Age")
  ) |>
  dplyr::mutate(
    CRPS_ratio = safe_div(.data$CRPS, .data$base_CRPS),
    CRPS_ratio = ifelse(.data$Method == "base", 1, .data$CRPS_ratio)
  )

write_csv_safe(crps_ratio_by_series_age, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_ratio_by_series_age.csv")))

crps_median_age_by_series <- crps_ratio_by_series_age |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$Level, .data$Series,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    median_CRPS = median(.data$CRPS, na.rm = TRUE),
    mean_CRPS = mean(.data$CRPS, na.rm = TRUE),
    IQR_CRPS = IQR(.data$CRPS, na.rm = TRUE),
    median_CRPS_ratio = median(.data$CRPS_ratio, na.rm = TRUE),
    mean_CRPS_ratio = mean(.data$CRPS_ratio, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(crps_median_age_by_series, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_median_age_by_series.csv")))

crps_table_like_rmse <- crps_median_age_by_series |>
  dplyr::select(
    .data$Structure_label, .data$Structure,
    .data$Level, .data$Series,
    .data$Method_label, .data$median_CRPS
  ) |>
  tidyr::pivot_wider(
    names_from = .data$Method_label,
    values_from = .data$median_CRPS
  ) |>
  dplyr::mutate(
    best_method = dplyr::case_when(
      .data$CoDA <= .data$OLS & .data$CoDA <= .data$WLS & .data$CoDA <= .data$MinT ~ "CoDA",
      .data$OLS  <= .data$CoDA & .data$OLS  <= .data$WLS & .data$OLS  <= .data$MinT ~ "OLS",
      .data$WLS  <= .data$CoDA & .data$WLS  <= .data$OLS & .data$WLS  <= .data$MinT ~ "WLS",
      TRUE ~ "MinT"
    ),
    OLS_over_CoDA = safe_div(.data$OLS, .data$CoDA),
    WLS_over_CoDA = safe_div(.data$WLS, .data$CoDA),
    MinT_over_CoDA = safe_div(.data$MinT, .data$CoDA)
  ) |>
  dplyr::arrange(.data$Structure, .data$Level, .data$Series)

write_csv_safe(crps_table_like_rmse, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_table_like_RMSE_median_age.csv")))

crps_summary_by_level <- crps_ratio_by_series_age |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$Level,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_CRPS = median(.data$CRPS, na.rm = TRUE),
    mean_CRPS = mean(.data$CRPS, na.rm = TRUE),
    median_CRPS_ratio = median(.data$CRPS_ratio, na.rm = TRUE),
    mean_CRPS_ratio = mean(.data$CRPS_ratio, na.rm = TRUE),
    IQR_CRPS_ratio = IQR(.data$CRPS_ratio, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$Structure, .data$Level, .data$Method)

write_csv_safe(crps_summary_by_level, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_summary_by_level.csv")))

crps_summary_overall <- crps_ratio_by_series_age |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    median_CRPS = median(.data$CRPS, na.rm = TRUE),
    mean_CRPS = mean(.data$CRPS, na.rm = TRUE),
    median_CRPS_ratio = median(.data$CRPS_ratio, na.rm = TRUE),
    mean_CRPS_ratio = mean(.data$CRPS_ratio, na.rm = TRUE),
    IQR_CRPS_ratio = IQR(.data$CRPS_ratio, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$Structure, .data$Method)

write_csv_safe(crps_summary_overall, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_summary_overall.csv")))

crps_main_table <- crps_summary_by_level |>
  dplyr::select(
    Engine,
    Structure = .data$Structure_label,
    Level = .data$Level,
    Method = .data$Method_label,
    n_cells,
    median_CRPS_ratio,
    mean_CRPS_ratio,
    IQR_CRPS_ratio,
    block_size,
    nsim
  )

write_csv_safe(crps_main_table, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_main_table_ratios_by_level.csv")))

######################################################################
## CRPS tables in same wide/LaTeX style as the existing CRPS table
######################################################################

CRPS_g1 <- list(
  CoDA = make_crps_matrix(crps_by_series_age, "g1", "base", structures$g1$row_order, ages),
  OLS  = make_crps_matrix(crps_by_series_age, "g1", "OLS",  structures$g1$row_order, ages),
  WLS  = make_crps_matrix(crps_by_series_age, "g1", "WLS",  structures$g1$row_order, ages),
  MinT = make_crps_matrix(crps_by_series_age, "g1", "MinT", structures$g1$row_order, ages)
)
CRPS_g2 <- list(
  CoDA = make_crps_matrix(crps_by_series_age, "g2", "base", structures$g2$row_order, ages),
  OLS  = make_crps_matrix(crps_by_series_age, "g2", "OLS",  structures$g2$row_order, ages),
  WLS  = make_crps_matrix(crps_by_series_age, "g2", "WLS",  structures$g2$row_order, ages),
  MinT = make_crps_matrix(crps_by_series_age, "g2", "MinT", structures$g2$row_order, ages)
)
CRPS_gT <- list(
  CoDA = make_crps_matrix(crps_by_series_age, "gT", "base", structures$gT$row_order, ages),
  OLS  = make_crps_matrix(crps_by_series_age, "gT", "OLS",  structures$gT$row_order, ages),
  WLS  = make_crps_matrix(crps_by_series_age, "gT", "WLS",  structures$gT$row_order, ages),
  MinT = make_crps_matrix(crps_by_series_age, "gT", "MinT", structures$gT$row_order, ages)
)

bb_crps_gT <- make_metric_matrix(CRPS_gT$CoDA, CRPS_gT$OLS, CRPS_gT$WLS, CRPS_gT$MinT, structures$gT$display_order)
bb_crps_g1 <- make_metric_matrix(CRPS_g1$CoDA, CRPS_g1$OLS, CRPS_g1$WLS, CRPS_g1$MinT, structures$g1$display_order)
bb_crps_g2 <- make_metric_matrix(CRPS_g2$CoDA, CRPS_g2$OLS, CRPS_g2$WLS, CRPS_g2$MinT, structures$g2$display_order)

bb_crps_full <- dplyr::bind_rows(
  bb_crps_gT |>
    dplyr::mutate(
      Hierarchy = "Full grouped structure",
      Level = c("Top", rep("Middle (cause)", 8), rep("Middle (sex)", 2), rep("Bottom", 16)),
      .before = Series
    ),
  bb_crps_g1 |>
    dplyr::mutate(
      Hierarchy = "Sex-based hierarchy",
      Level = c("Top", rep("Middle (sex)", 2), rep("Bottom", 16)),
      .before = Series
    ),
  bb_crps_g2 |>
    dplyr::mutate(
      Hierarchy = "Cause-based hierarchy",
      Level = c("Top", rep("Middle (cause)", 8), rep("Bottom", 16)),
      .before = Series
    )
)

write_csv_safe(bb_crps_gT, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_full_grouped.csv")))
write_csv_safe(bb_crps_g1, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_sex_hierarchy.csv")))
write_csv_safe(bb_crps_g2, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_cause_hierarchy.csv")))
write_csv_safe(bb_crps_full, file.path(TABLE_DIR, paste0(BB_PREFIX, "CRPS_all_structures.csv")))

if (exists("write_manuscript_table")) {
  write_manuscript_table(
    bb_crps_gT,
    bb_crps_g1,
    bb_crps_g2,
    "CRPS",
    file.path(TABLE_DIR, paste0(BB_PREFIX, "interval_CRPS_table.tex"))
  )
}

######################################################################
## Prediction interval summaries
######################################################################

pi_by_series_age <- pi_long |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$Level, .data$Series, .data$Age,
    .data$alpha, .data$nominal_coverage,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    empirical_coverage = mean(.data$covered, na.rm = TRUE),
    mean_width = mean(.data$width, na.rm = TRUE),
    median_width = median(.data$width, na.rm = TRUE),
    mean_interval_score = mean(.data$interval_score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(pi_by_series_age, file.path(TABLE_DIR, paste0(BB_PREFIX, "PI_by_series_age.csv")))

pi_summary_by_level <- pi_by_series_age |>
  dplyr::group_by(
    .data$Engine,
    .data$Structure, .data$Structure_label,
    .data$Method, .data$Method_label,
    .data$Level,
    .data$alpha, .data$nominal_coverage,
    .data$block_size, .data$nsim
  ) |>
  dplyr::summarise(
    median_empirical_coverage = median(.data$empirical_coverage, na.rm = TRUE),
    mean_empirical_coverage = mean(.data$empirical_coverage, na.rm = TRUE),
    median_width = median(.data$median_width, na.rm = TRUE),
    mean_width = mean(.data$mean_width, na.rm = TRUE),
    median_interval_score = median(.data$mean_interval_score, na.rm = TRUE),
    mean_interval_score = mean(.data$mean_interval_score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(pi_summary_by_level, file.path(TABLE_DIR, paste0(BB_PREFIX, "PI_summary_by_level.csv")))

######################################################################
## Optional comparison with existing model-based CoDA bootstrap CRPS table
######################################################################

model_crps_path <- file.path(TABLE_DIR, "CRPS_all_structures.csv")
if (file.exists(model_crps_path)) {
  model_crps <- readr::read_csv(model_crps_path, show_col_types = FALSE) |>
    dplyr::rename_with(~ paste0("model_based_", .x), c("CoDA", "OLS", "WLS", "MinT"))

  compare_model_vs_bb <- bb_crps_full |>
    dplyr::rename_with(~ paste0("residual_BB_", .x), c("CoDA", "OLS", "WLS", "MinT")) |>
    dplyr::left_join(model_crps, by = c("Hierarchy", "Level", "Series")) |>
    dplyr::mutate(
      CoDA_resBB_over_model = safe_div(.data$residual_BB_CoDA, .data$model_based_CoDA),
      OLS_resBB_over_model  = safe_div(.data$residual_BB_OLS,  .data$model_based_OLS),
      WLS_resBB_over_model  = safe_div(.data$residual_BB_WLS,  .data$model_based_WLS),
      MinT_resBB_over_model = safe_div(.data$residual_BB_MinT, .data$model_based_MinT)
    )

  write_csv_safe(compare_model_vs_bb, file.path(TABLE_DIR, paste0(BB_PREFIX, "COMPARE_model_based_vs_residualBB_median_CRPS.csv")))
}

######################################################################
## Save compact RDS bundle
######################################################################

save_rds_safe(
  list(
    CRPS_g1 = CRPS_g1,
    CRPS_g2 = CRPS_g2,
    CRPS_gT = CRPS_gT,
    bb_crps_g1 = bb_crps_g1,
    bb_crps_g2 = bb_crps_g2,
    bb_crps_gT = bb_crps_gT,
    bb_crps_full = bb_crps_full,
    crps_by_series_age = crps_by_series_age,
    crps_ratio_by_series_age = crps_ratio_by_series_age,
    crps_summary_by_level = crps_summary_by_level,
    pi_summary_by_level = pi_summary_by_level,
    negative_diagnostics = negative_diagnostics,
    options = list(
      BB_NSIM = BB_NSIM,
      BB_BLOCK_SIZE = BB_BLOCK_SIZE,
      BB_CENTER_RESIDUALS = BB_CENTER_RESIDUALS,
      BB_TRUNCATE_BASE_AT_ZERO = BB_TRUNCATE_BASE_AT_ZERO,
      BB_NONNEGATIVE_RECONCILED_POSTPROCESS = BB_NONNEGATIVE_RECONCILED_POSTPROCESS,
      BB_USE_FORECO_MAPPING = BB_USE_FORECO_MAPPING,
      BB_MINT_SHRINK_LAMBDA = BB_MINT_SHRINK_LAMBDA,
      BB_SEED = BB_SEED
    )
  ),
  file.path(RDS_DIR, paste0(BB_PREFIX, "residual_block_bootstrap_results.rds"))
)

cat("\nDone. Residual block-bootstrap CRPS tables and figures created.\n")
cat("Main tables:\n")
cat(" - outputs/tables/", BB_PREFIX, "CRPS_all_structures.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "interval_CRPS_table.tex\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "CRPS_table_like_RMSE_median_age.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "CRPS_main_table_ratios_by_level.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "CRPS_summary_by_level.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "CRPS_ratio_by_series_age.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "PI_summary_by_level.csv\n", sep = "")
cat(" - outputs/tables/", BB_PREFIX, "negative_diagnostics.csv\n", sep = "")
