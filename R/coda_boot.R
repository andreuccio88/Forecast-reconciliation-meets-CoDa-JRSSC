######################################################################
## R/bootstrap_coda.R
##
## Model-based bootstrap for legacy CoDA forecasts.
##
## Implements the Appendix B logic:
## 1. Fit CoDA in clr space.
## 2. Resample clr residuals to account for parameter uncertainty.
## 3. Re-estimate SVD parameters on bootstrap fitted clr tables.
## 4. Simulate future kappa paths via ARIMA residual bootstrap.
## 5. Transform back with clrInv and perturbation by alpha.
##
## Output:
## - CoDA bootstrap paths for g1, g2, and gT structures.
## - Arrays have dimensions:
##     B x H x n_series
##   where B = n_epsilon * n_kappa.
######################################################################


######################################################################
## Basic helpers
######################################################################

close_rows <- function(x, radix = 100000) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  x[x <= 0] <- 0.25
  
  rs <- rowSums(x, na.rm = TRUE)
  if (any(!is.finite(rs) | rs <= 0)) {
    stop("Some rows cannot be closed because their sum is non-positive.", call. = FALSE)
  }
  
  sweep(x, 1, rs, FUN = "/") * radix
}


resample_clr_residual_table <- function(E) {
  ## Appendix B says that a new residual table is created by assigning
  ## randomly chosen row and column residuals to each time-age cell.
  ## Since clr rows must sum to zero, we re-centre each simulated row.
  
  E <- as.matrix(E)
  n <- nrow(E)
  q <- ncol(E)
  
  row_id <- matrix(sample.int(n, n * q, replace = TRUE), nrow = n, ncol = q)
  col_id <- matrix(sample.int(q, n * q, replace = TRUE), nrow = n, ncol = q)
  
  E_star <- matrix(E[cbind(as.vector(row_id), as.vector(col_id))],
                   nrow = n, ncol = q)
  
  E_star <- sweep(E_star, 1, rowMeans(E_star, na.rm = TRUE), FUN = "-")
  E_star
}


safe_auto_arima <- function(x) {
  tryCatch(
    forecast::auto.arima(x, max.d = 1),
    error = function(e) {
      forecast::Arima(x, order = c(0, 1, 0), include.drift = TRUE)
    }
  )
}


safe_simulate_arima <- function(model, h) {
  out <- tryCatch(
    as.numeric(stats::simulate(model, nsim = h, future = TRUE, bootstrap = TRUE)),
    error = function(e) {
      as.numeric(forecast::forecast(model, h = h)$mean)
    }
  )
  
  if (length(out) != h || any(!is.finite(out))) {
    out <- as.numeric(forecast::forecast(model, h = h)$mean)
  }
  
  out
}


######################################################################
## Legacy data constructor: reproduces the old CoDa() input matrix
######################################################################

build_legacy_coda_input <- function(dt,
                                    sex = c("M", "F", "T"),
                                    hstr = c(1, 2),
                                    ih = 15,
                                    years.fit = 1970:2000,
                                    lt_dir = "LT") {
  
  sex <- match.arg(sex)
  hstr <- as.integer(match.arg(as.character(hstr), choices = c("1", "2")))
  
  ages <- as.numeric(colnames(dt[, -c(1:3)]))
  m <- length(ages)
  years <- as.numeric(unique(dt$Year))
  n <- length(years)
  
  CoD <- unique(dt$Cause_Rev)
  k <- length(CoD)
  
  nomiCoD <- c(
    "Infectious diseases",
    "Cancer Smoking",
    "Cancer No Smoking",
    "Diseases of the circulatory system",
    "Diseases of the respiratory system",
    "External diseases",
    "Diseases of the digestive system",
    "Other diseases"
  )
  
  if (k != 8) {
    stop("Expected 8 causes in Cause_Rev, found ", k, ".", call. = FALSE)
  }
  
  dati0 <- as.matrix(subset(dt, Sex == sex)[, 1:m + 3])
  
  dati1 <- array(
    c(t(dati0)),
    dim = c(m, n, k),
    dimnames = list(ages, years, nomiCoD)
  )
  
  ## Legacy code keeps ages 25,...,80+
  dati <- dati1[7:18, , ]
  ages2 <- ages[7:18]
  
  lt_file <- switch(
    sex,
    "M" = file.path(lt_dir, "mltper_5x1.txt"),
    "F" = file.path(lt_dir, "fltper_5x1.txt"),
    "T" = file.path(lt_dir, "bltper_5x1.txt")
  )
  
  if (!file.exists(lt_file)) {
    stop("Cannot find LT file: ", lt_file, call. = FALSE)
  }
  
  LT.1 <- read.table(lt_file, skip = 1, header = TRUE, stringsAsFactors = FALSE)
  
  dx.nat <- matrix(
    LT.1$dx,
    length(unique(LT.1$Age)),
    length(unique(LT.1$Year))
  )
  
  colnames(dx.nat) <- unique(LT.1$Year)
  rownames(dx.nat) <- unique(LT.1$Age)
  
  dx.nat1 <- dx.nat[, colnames(dx.nat) %in% years, drop = FALSE]
  rownames(dx.nat1) <- unique(LT.1$Age)
  
  ## Legacy aggregation into 80+
  dx.nat3 <- rbind(
    dx.nat1[7:17, , drop = FALSE],
    "80+" = dx.nat1[18, ] +
      dx.nat1[19, ] +
      dx.nat1[20, ] +
      dx.nat1[21, ] +
      dx.nat1[22, ] +
      dx.nat1[23, ] +
      dx.nat1[24, ]
  )
  
  dati.total <- dati[, , 1]
  for (jj in 2:8) {
    dati.total <- dati.total + dati[, , jj]
  }
  
  dati.rela1 <- array(
    NA_real_,
    dim = c(ncol(dati[, , 1]), nrow(dati[, , 1]), 8),
    dimnames = list(years, ages2, nomiCoD)
  )
  
  dati.temp <- dati * 0
  
  for (i in 1:8) {
    dati.temp[, , i] <- ifelse(dati[, , i] == 0, 0.25, dati[, , i])
    dati.rela1[, , i] <- t(dx.nat3 * (dati.temp[, , i] / dati.total))
  }
  
  ## Legacy closure to radix 100000
  for (j in 1:n) {
    somma <- sum(dati.rela1[j, , ], na.rm = TRUE)
    if (!is.finite(somma) || somma <= 0) {
      stop("Non-positive yearly composition for year ", years[j], ".", call. = FALSE)
    }
    dati.rela1[j, , ] <- dati.rela1[j, , ] * (100000 / somma)
  }
  
  years.fitfor <- c(years.fit, (max(years.fit) + 1):(max(years.fit) + ih))
  years.for <- (max(years.fit) + 1):(max(years.fit) + ih)
  
  idx_fit <- match(years.fit, years)
  idx_for <- match(years.for, years)
  
  if (any(is.na(idx_fit)) || any(is.na(idx_for))) {
    stop("The data do not contain all fitting and forecast years.", call. = FALSE)
  }
  
  if (hstr == 1) {
    dx.com <- do.call(
      cbind,
      lapply(1:8, function(i) dati.rela1[, , i])
    )
    
    colnames(dx.com) <- unlist(
      lapply(nomiCoD, function(cc) paste(ages2, cc, sep = "."))
    )
    
    obs_array <- dati.rela1
    causes <- nomiCoD
  }
  
  if (hstr == 2) {
    allcause <- dati.rela1[, , 1]
    for (i in 2:8) {
      allcause <- allcause + dati.rela1[, , i]
    }
    
    colnames(allcause) <- paste(colnames(allcause), "All-causes", sep = ".")
    
    dx.com <- allcause
    obs_array <- array(
      allcause,
      dim = c(length(years), length(ages2), 1),
      dimnames = list(years, ages2, "All-causes")
    )
    
    causes <- "All-causes"
  }
  
  list(
    dx.com = dx.com,
    obs_array = obs_array,
    years = years,
    years.fit = years.fit,
    years.for = years.for,
    idx_fit = idx_fit,
    idx_for = idx_for,
    ages = ages2,
    causes = causes,
    sex = sex,
    hstr = hstr
  )
}


######################################################################
## Core CoDA bootstrap in clr space
######################################################################

coda_ct_bootstrap <- function(dx,
                              ih,
                              k,
                              years,
                              ages,
                              ses,
                              nr_rank = 3,
                              n_epsilon = 100,
                              n_kappa = 100,
                              seed = 1234) {
  
  set.seed(seed)
  
  dx <- close_rows(dx, radix = 100000)
  
  n <- length(years)
  m <- length(ages)
  q <- m * k
  B <- n_epsilon * n_kappa
  
  if (nrow(dx) != n || ncol(dx) != q) {
    stop("Incorrect dx dimensions in coda_ct_bootstrap().", call. = FALSE)
  }
  
  ## Original CoDA fit
  close.dx <- compositions::acomp(dx)
  ax <- compositions::geometricmeanCol(close.dx)
  dx.cent <- close.dx - ax
  clr.cent <- as.matrix(unclass(compositions::clr(dx.cent)))
  
  sv <- svd(clr.cent, nu = nr_rank, nv = nr_rank)
  
  bx <- sv$v[, seq_len(nr_rank), drop = FALSE]
  kt <- sweep(
    sv$u[, seq_len(nr_rank), drop = FALSE],
    2,
    sv$d[seq_len(nr_rank)],
    FUN = "*"
  )
  
  fit_clr <- kt %*% t(bx)
  res_clr <- clr.cent - fit_clr
  
  ## Fitted dx from original model, useful for diagnostics
  fit_acomp <- compositions::clrInv(fit_clr) + ax
  dx_fit <- array(
    unclass(fit_acomp) * 100000,
    dim = c(n, m, k),
    dimnames = list(as.character(years), as.character(ages), ses)
  )
  
  dx_obs <- array(
    dx,
    dim = c(n, m, k),
    dimnames = list(as.character(years), as.character(ages), ses)
  )
  
  ## Keep the original ARIMA orders as reference
  arima_original <- vector("list", nr_rank)
  for (r in seq_len(nr_rank)) {
    arima_original[[r]] <- safe_auto_arima(kt[, r])
  }
  
  dx_boot <- array(
    NA_real_,
    dim = c(ih, m, k, B),
    dimnames = list(
      as.character((max(years) + 1):(max(years) + ih)),
      as.character(ages),
      ses,
      paste0("b", seq_len(B))
    )
  )
  
  boot_id <- 0
  
  for (b_eps in seq_len(n_epsilon)) {
    
    ## Parameter uncertainty:
    ## simulate a new fitted clr table and re-estimate CoDA parameters.
    E_star <- resample_clr_residual_table(res_clr)
    clr_star <- fit_clr + E_star
    
    sv_star <- svd(clr_star, nu = nr_rank, nv = nr_rank)
    
    bx_star <- sv_star$v[, seq_len(nr_rank), drop = FALSE]
    kt_star <- sweep(
      sv_star$u[, seq_len(nr_rank), drop = FALSE],
      2,
      sv_star$d[seq_len(nr_rank)],
      FUN = "*"
    )
    
    ## Sign alignment, otherwise SVD flips make ARIMA paths look possessed.
    for (r in seq_len(nr_rank)) {
      if (sum(bx_star[, r] * bx[, r], na.rm = TRUE) < 0) {
        bx_star[, r] <- -bx_star[, r]
        kt_star[, r] <- -kt_star[, r]
      }
    }
    
    arima_star <- vector("list", nr_rank)
    
    for (r in seq_len(nr_rank)) {
      arima_star[[r]] <- safe_auto_arima(kt_star[, r])
    }
    
    for (b_kap in seq_len(n_kappa)) {
      
      boot_id <- boot_id + 1
      
      future_clr <- matrix(0, nrow = ih, ncol = q)
      
      ## Time-index extrapolation uncertainty
      for (r in seq_len(nr_rank)) {
        kt_future <- safe_simulate_arima(arima_star[[r]], ih)
        future_clr <- future_clr + kt_future %o% bx_star[, r]
      }
      
      future_acomp <- compositions::clrInv(future_clr) + ax
      future_dx <- unclass(future_acomp) * 100000
      
      dx_boot[, , , boot_id] <- array(
        future_dx,
        dim = c(ih, m, k)
      )
    }
  }
  
  list(
    dx_boot = dx_boot,
    dx_fit = dx_fit,
    dx_obs = dx_obs,
    ax = ax,
    bx = bx,
    kt = kt,
    fit_clr = fit_clr,
    res_clr = res_clr,
    sv = sv,
    arima_original = arima_original,
    n_epsilon = n_epsilon,
    n_kappa = n_kappa,
    B = B,
    nr_rank = nr_rank
  )
}


######################################################################
## Wrapper: legacy CoDA bootstrap for one sex and one hierarchy level
######################################################################

CoDa_boot_legacy <- function(dt,
                             sex = c("M", "F", "T"),
                             hstr = c(1, 2),
                             ih = 15,
                             years.fit = 1970:2000,
                             lt_dir = "LT",
                             nr_rank = 3,
                             n_epsilon = 100,
                             n_kappa = 100,
                             seed = 1234) {
  
  sex <- match.arg(sex)
  hstr <- as.integer(match.arg(as.character(hstr), choices = c("1", "2")))
  
  comp <- build_legacy_coda_input(
    dt = dt,
    sex = sex,
    hstr = hstr,
    ih = ih,
    years.fit = years.fit,
    lt_dir = lt_dir
  )
  
  fit_years <- comp$years.fit
  
  boot <- coda_ct_bootstrap(
    dx = comp$dx.com[comp$idx_fit, , drop = FALSE],
    ih = ih,
    k = length(comp$causes),
    years = fit_years,
    ages = comp$ages,
    ses = comp$causes,
    nr_rank = nr_rank,
    n_epsilon = n_epsilon,
    n_kappa = n_kappa,
    seed = seed
  )
  
  actual_oos <- comp$obs_array[as.character(comp$years.for), , , drop = FALSE]
  actual_fit <- comp$obs_array[as.character(comp$years.fit), , , drop = FALSE]
  
  list(
    boot = boot,
    actual_oos = actual_oos,
    actual_fit = actual_fit,
    ages = comp$ages,
    causes = comp$causes,
    years.fit = comp$years.fit,
    years.for = comp$years.for,
    sex = sex,
    hstr = hstr
  )
}


######################################################################
## Series names for reconciliation structures
######################################################################

build_bootstrap_series_names <- function() {
  
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


extract_boot_slice <- function(obj, age_id, cause_id) {
  ## obj$boot$dx_boot has dimension H x age x cause x B.
  ## Return B x H.
  t(obj$boot$dx_boot[, age_id, cause_id, ])
}


extract_actual_slice <- function(obj, age_id, cause_id) {
  ## obj$actual_oos has dimension H x age x cause.
  as.numeric(obj$actual_oos[, age_id, cause_id])
}


extract_fit_slice <- function(obj, age_id, cause_id) {
  ## obj$boot$dx_fit has dimension n_fit x age x cause.
  as.numeric(obj$boot$dx_fit[, age_id, cause_id])
}


extract_obs_fit_slice <- function(obj, age_id, cause_id) {
  ## obj$actual_fit has dimension n_fit x age x cause.
  as.numeric(obj$actual_fit[, age_id, cause_id])
}


######################################################################
## Assemble bootstrap paths for g1, g2, gT
######################################################################

assemble_coda_bootstrap_structures <- function(bl_M,
                                               bl_F,
                                               ml_M,
                                               ml_F,
                                               ml_T,
                                               tl_T) {
  
  nm <- build_bootstrap_series_names()
  
  ages <- tl_T$ages
  H <- length(tl_T$years.for)
  B <- tl_T$boot$B
  n_ages <- length(ages)
  
  out_g1 <- vector("list", n_ages)
  out_g2 <- vector("list", n_ages)
  out_gT <- vector("list", n_ages)
  
  actual_g1 <- vector("list", n_ages)
  actual_g2 <- vector("list", n_ages)
  actual_gT <- vector("list", n_ages)
  
  residuals_g1 <- vector("list", n_ages)
  residuals_g2 <- vector("list", n_ages)
  residuals_gT <- vector("list", n_ages)
  
  for (a in seq_along(ages)) {
    
    ##################################################################
    ## g1: sex-based hierarchy
    ##################################################################
    
    draws1 <- array(
      NA_real_,
      dim = c(B, H, length(nm$rows_g1)),
      dimnames = list(NULL, as.character(tl_T$years.for), nm$rows_g1)
    )
    
    actual1 <- matrix(
      NA_real_,
      nrow = H,
      ncol = length(nm$rows_g1),
      dimnames = list(as.character(tl_T$years.for), nm$rows_g1)
    )
    
    fit1 <- obs1 <- matrix(
      NA_real_,
      nrow = length(tl_T$years.fit),
      ncol = length(nm$rows_g1),
      dimnames = list(as.character(tl_T$years.fit), nm$rows_g1)
    )
    
    draws1[, , "Total"] <- extract_boot_slice(tl_T, a, 1)
    actual1[, "Total"] <- extract_actual_slice(tl_T, a, 1)
    fit1[, "Total"] <- extract_fit_slice(tl_T, a, 1)
    obs1[, "Total"] <- extract_obs_fit_slice(tl_T, a, 1)
    
    draws1[, , "Total - Males"] <- extract_boot_slice(ml_M, a, 1)
    actual1[, "Total - Males"] <- extract_actual_slice(ml_M, a, 1)
    fit1[, "Total - Males"] <- extract_fit_slice(ml_M, a, 1)
    obs1[, "Total - Males"] <- extract_obs_fit_slice(ml_M, a, 1)
    
    draws1[, , "Total - Females"] <- extract_boot_slice(ml_F, a, 1)
    actual1[, "Total - Females"] <- extract_actual_slice(ml_F, a, 1)
    fit1[, "Total - Females"] <- extract_fit_slice(ml_F, a, 1)
    obs1[, "Total - Females"] <- extract_obs_fit_slice(ml_F, a, 1)
    
    for (i in 1:8) {
      sM <- paste0("Cause ", i, " - Males")
      sF <- paste0("Cause ", i, " - Females")
      
      draws1[, , sM] <- extract_boot_slice(bl_M, a, i)
      actual1[, sM] <- extract_actual_slice(bl_M, a, i)
      fit1[, sM] <- extract_fit_slice(bl_M, a, i)
      obs1[, sM] <- extract_obs_fit_slice(bl_M, a, i)
      
      draws1[, , sF] <- extract_boot_slice(bl_F, a, i)
      actual1[, sF] <- extract_actual_slice(bl_F, a, i)
      fit1[, sF] <- extract_fit_slice(bl_F, a, i)
      obs1[, sF] <- extract_obs_fit_slice(bl_F, a, i)
    }
    
    out_g1[[a]] <- draws1
    actual_g1[[a]] <- actual1
    residuals_g1[[a]] <- fit1 - obs1
    
    
    ##################################################################
    ## g2: cause-based hierarchy
    ##################################################################
    
    draws2 <- array(
      NA_real_,
      dim = c(B, H, length(nm$rows_g2)),
      dimnames = list(NULL, as.character(tl_T$years.for), nm$rows_g2)
    )
    
    actual2 <- matrix(
      NA_real_,
      nrow = H,
      ncol = length(nm$rows_g2),
      dimnames = list(as.character(tl_T$years.for), nm$rows_g2)
    )
    
    fit2 <- obs2 <- matrix(
      NA_real_,
      nrow = length(tl_T$years.fit),
      ncol = length(nm$rows_g2),
      dimnames = list(as.character(tl_T$years.fit), nm$rows_g2)
    )
    
    draws2[, , "Total"] <- extract_boot_slice(tl_T, a, 1)
    actual2[, "Total"] <- extract_actual_slice(tl_T, a, 1)
    fit2[, "Total"] <- extract_fit_slice(tl_T, a, 1)
    obs2[, "Total"] <- extract_obs_fit_slice(tl_T, a, 1)
    
    for (i in 1:8) {
      sT <- paste0("Cause ", i, " - Total")
      sM <- paste0("Cause ", i, " - Males")
      sF <- paste0("Cause ", i, " - Females")
      
      draws2[, , sT] <- extract_boot_slice(ml_T, a, i)
      actual2[, sT] <- extract_actual_slice(ml_T, a, i)
      fit2[, sT] <- extract_fit_slice(ml_T, a, i)
      obs2[, sT] <- extract_obs_fit_slice(ml_T, a, i)
      
      draws2[, , sM] <- extract_boot_slice(bl_M, a, i)
      actual2[, sM] <- extract_actual_slice(bl_M, a, i)
      fit2[, sM] <- extract_fit_slice(bl_M, a, i)
      obs2[, sM] <- extract_obs_fit_slice(bl_M, a, i)
      
      draws2[, , sF] <- extract_boot_slice(bl_F, a, i)
      actual2[, sF] <- extract_actual_slice(bl_F, a, i)
      fit2[, sF] <- extract_fit_slice(bl_F, a, i)
      obs2[, sF] <- extract_obs_fit_slice(bl_F, a, i)
    }
    
    out_g2[[a]] <- draws2
    actual_g2[[a]] <- actual2
    residuals_g2[[a]] <- fit2 - obs2
    
    
    ##################################################################
    ## gT: full grouped structure
    ##################################################################
    
    drawsT <- array(
      NA_real_,
      dim = c(B, H, length(nm$rows_gT)),
      dimnames = list(NULL, as.character(tl_T$years.for), nm$rows_gT)
    )
    
    actualT <- matrix(
      NA_real_,
      nrow = H,
      ncol = length(nm$rows_gT),
      dimnames = list(as.character(tl_T$years.for), nm$rows_gT)
    )
    
    fitT <- obsT <- matrix(
      NA_real_,
      nrow = length(tl_T$years.fit),
      ncol = length(nm$rows_gT),
      dimnames = list(as.character(tl_T$years.fit), nm$rows_gT)
    )
    
    drawsT[, , "Total"] <- extract_boot_slice(tl_T, a, 1)
    actualT[, "Total"] <- extract_actual_slice(tl_T, a, 1)
    fitT[, "Total"] <- extract_fit_slice(tl_T, a, 1)
    obsT[, "Total"] <- extract_obs_fit_slice(tl_T, a, 1)
    
    drawsT[, , "Total - Males"] <- extract_boot_slice(ml_M, a, 1)
    actualT[, "Total - Males"] <- extract_actual_slice(ml_M, a, 1)
    fitT[, "Total - Males"] <- extract_fit_slice(ml_M, a, 1)
    obsT[, "Total - Males"] <- extract_obs_fit_slice(ml_M, a, 1)
    
    drawsT[, , "Total - Females"] <- extract_boot_slice(ml_F, a, 1)
    actualT[, "Total - Females"] <- extract_actual_slice(ml_F, a, 1)
    fitT[, "Total - Females"] <- extract_fit_slice(ml_F, a, 1)
    obsT[, "Total - Females"] <- extract_obs_fit_slice(ml_F, a, 1)
    
    for (i in 1:8) {
      sT <- paste0("Cause ", i, " - Total")
      sM <- paste0("Cause ", i, " - Males")
      sF <- paste0("Cause ", i, " - Females")
      
      drawsT[, , sT] <- extract_boot_slice(ml_T, a, i)
      actualT[, sT] <- extract_actual_slice(ml_T, a, i)
      fitT[, sT] <- extract_fit_slice(ml_T, a, i)
      obsT[, sT] <- extract_obs_fit_slice(ml_T, a, i)
      
      drawsT[, , sM] <- extract_boot_slice(bl_M, a, i)
      actualT[, sM] <- extract_actual_slice(bl_M, a, i)
      fitT[, sM] <- extract_fit_slice(bl_M, a, i)
      obsT[, sM] <- extract_obs_fit_slice(bl_M, a, i)
      
      drawsT[, , sF] <- extract_boot_slice(bl_F, a, i)
      actualT[, sF] <- extract_actual_slice(bl_F, a, i)
      fitT[, sF] <- extract_fit_slice(bl_F, a, i)
      obsT[, sF] <- extract_obs_fit_slice(bl_F, a, i)
    }
    
    out_gT[[a]] <- drawsT
    actual_gT[[a]] <- actualT
    residuals_gT[[a]] <- fitT - obsT
  }
  
  names(out_g1) <- names(out_g2) <- names(out_gT) <- as.character(ages)
  names(actual_g1) <- names(actual_g2) <- names(actual_gT) <- as.character(ages)
  names(residuals_g1) <- names(residuals_g2) <- names(residuals_gT) <- as.character(ages)
  
  list(
    names = nm,
    ages = ages,
    years.for = tl_T$years.for,
    years.fit = tl_T$years.fit,
    B = B,
    g1 = list(draws = out_g1, actual = actual_g1, residuals = residuals_g1),
    g2 = list(draws = out_g2, actual = actual_g2, residuals = residuals_g2),
    gT = list(draws = out_gT, actual = actual_gT, residuals = residuals_gT)
  )
}


######################################################################
## Main wrapper: run all six CoDA bootstraps and assemble structures
######################################################################

generate_coda_bootstrap_paths <- function(dt,
                                          ih = 15,
                                          years.fit = 1970:2000,
                                          lt_dir = "LT",
                                          nr_rank = 3,
                                          n_epsilon = 100,
                                          n_kappa = 100,
                                          seed = 1234,
                                          save_file = NULL) {
  
  message("Bootstrap CoDA: bottom-level M by cause")
  bl_M <- CoDa_boot_legacy(
    dt = dt, sex = "M", hstr = 1, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 1
  )
  
  message("Bootstrap CoDA: bottom-level F by cause")
  bl_F <- CoDa_boot_legacy(
    dt = dt, sex = "F", hstr = 1, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 2
  )
  
  message("Bootstrap CoDA: middle-level M all-cause")
  ml_M <- CoDa_boot_legacy(
    dt = dt, sex = "M", hstr = 2, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 3
  )
  
  message("Bootstrap CoDA: middle-level F all-cause")
  ml_F <- CoDa_boot_legacy(
    dt = dt, sex = "F", hstr = 2, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 4
  )
  
  message("Bootstrap CoDA: middle-level T by cause")
  ml_T <- CoDa_boot_legacy(
    dt = dt, sex = "T", hstr = 1, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 5
  )
  
  message("Bootstrap CoDA: top-level T all-cause")
  tl_T <- CoDa_boot_legacy(
    dt = dt, sex = "T", hstr = 2, ih = ih, years.fit = years.fit,
    lt_dir = lt_dir, nr_rank = nr_rank,
    n_epsilon = n_epsilon, n_kappa = n_kappa,
    seed = seed + 6
  )
  
  message("Assembling bootstrap paths for g1, g2, gT")
  
  assembled <- assemble_coda_bootstrap_structures(
    bl_M = bl_M,
    bl_F = bl_F,
    ml_M = ml_M,
    ml_F = ml_F,
    ml_T = ml_T,
    tl_T = tl_T
  )
  
  assembled$raw_boot <- list(
    bl_M = bl_M,
    bl_F = bl_F,
    ml_M = ml_M,
    ml_F = ml_F,
    ml_T = ml_T,
    tl_T = tl_T
  )
  
  if (!is.null(save_file)) {
    dir.create(dirname(save_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(assembled, save_file)
    message("Saved bootstrap paths to: ", save_file)
  }
  
  assembled
}