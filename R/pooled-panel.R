# Pooled panel cointegrating polynomial regression with fixed effects
# (de Jong & Wagner, 2016). Port of PanelEKC_indiv_eff_only.m ("oneway":
# individual fixed effects only) and PanelEKC_two_eff.m ("twoway":
# individual + time fixed effects).
#
# This is structurally different from pcpr(type = "mg"): mg estimates N
# separate slopes (one per unit) and averages them (full slope
# heterogeneity). This estimator assumes a *single common slope* shared by
# every unit, with only the intercepts (and, under "twoway", a common time
# effect) allowed to vary. Long-run variances are still estimated per unit
# (own bandwidth, own kernel-based Omega_i/Delta_i), but are then pooled
# (averaged across units) into one Omega_bar/Delta_bar used to build a
# single bias correction applied uniformly to every unit.
#
# Two point estimators come out of this: `beta_Mod` (a closed-form
# second-order bias-corrected LSDV/within estimator -- "Modified OLS", not
# the same construction as FM-OLS) and `beta_FM` (the actual Fully Modified
# analogue). The correction algebra uses theoretical constant matrices `M`,
# `Q` (moments of the limiting Brownian-motion functionals for a degree-q
# polynomial regressor) and `GT` (T-power scaling), tabulated in the
# original source only for q = 2 and q = 3 with a single integrated
# regressor -- this port keeps that restriction rather than attempting to
# generalize matrices that were not derived for the general case.
#
# Note this estimator does *not* drop the first time observation (unlike
# cpr()/pcpr(type="mg"), which difference and truncate to T-1): it sets
# v_{i,1} = x_{i,1} (as if x_{i,0} = 0) and uses the full T observations.
# That is a genuine methodological difference from the mg estimator, not
# an inconsistency to paper over.

#' Within-transformation (demeaning) for the pooled panel CPR estimator
#'
#' Port of simpledemean.m.
#'
#' @param y,x `T x N` matrices.
#' @param way `"oneway"` (individual fixed effects only) or `"twoway"`
#'   (individual + time fixed effects).
#' @param q Highest power of `x` to include (2 or 3).
#' @keywords internal
simpledemean <- function(y, x, way = c("oneway", "twoway"), q) {
  way <- match.arg(way)
  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn <- nrow(x)
  N <- ncol(x)

  x2 <- x^2
  x3 <- if (q == 3) x^3 else NULL

  if (way == "oneway") {
    xmeans <- colMeans(x)
    x2means <- colMeans(x2)
    ymeans <- colMeans(y)

    xtilde <- sweep(x, 2, xmeans, "-")
    xquad_tilde <- sweep(x2, 2, x2means, "-")
    ytilde <- sweep(y, 2, ymeans, "-")

    if (q == 3) {
      x3means <- colMeans(x3)
      xcub_tilde <- sweep(x3, 2, x3means, "-")
    } else {
      xcub_tilde <- NULL
    }
  } else { # twoway
    demean2 <- function(m) {
      col_mean <- matrix(colMeans(m), Tn, N, byrow = TRUE)
      row_mean <- matrix(rowMeans(m), Tn, N)
      m - col_mean - row_mean + sum(m) / (Tn * N)
    }
    xtilde <- demean2(x)
    xquad_tilde <- demean2(x2)
    ytilde <- demean2(y)
    xcub_tilde <- if (q == 3) demean2(x3) else NULL
  }

  list(ytilde = ytilde, xtilde = xtilde, xquad_tilde = xquad_tilde, xcub_tilde = xcub_tilde)
}

#' @keywords internal
.pooled_theory_matrices <- function(q) {
  if (q == 2) {
    list(
      M = matrix(c(1/6, 0, 0, 5/12), 2, 2, byrow = TRUE),
      Q = matrix(c(1/3, 0, 0, 59/60), 2, 2, byrow = TRUE)
    )
  } else if (q == 3) {
    list(
      M = matrix(c(1/6, 0, 3/8,
                   0,   5/12, 0,
                   3/8, 0, 39/20), 3, 3, byrow = TRUE),
      Q = matrix(c(1/3, 0, 9/10,
                   0,   59/60, 0,
                   9/10, 0, 101/20), 3, 3, byrow = TRUE)
    )
  } else {
    stop("The pooled panel estimator only supports q (max power of the single ",
         "integrated regressor) = 2 or 3 -- the theoretical bias-correction ",
         "matrices M, Q, GT are only tabulated for those cases in the original ",
         "de Jong & Wagner (2016) source.", call. = FALSE)
  }
}

#' Pooled panel cointegrating polynomial regression (de Jong & Wagner, 2016)
#'
#' A single common slope shared by all units (individual, and optionally
#' time, fixed effects removed by a within-transformation), with per-unit
#' long-run variances pooled into one bias correction. See the file-level
#' comment in `R/pooled-panel.R` for the full derivation and its known
#' restrictions (single regressor, `q` in `{2, 3}`).
#'
#' @param y,x `T x N` matrices (balanced panel: T time points, N units).
#' @param q Highest power of `x` to include; must be `2` or `3`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param bandwidth Bandwidth selection: `"And91"`, `"AM92"`, `"NW"`, or a
#'   fixed numeric value. Applied per unit (own Omega_i/Delta_i), as in
#'   [pcpr()]'s `"mg"` type.
#' @param effects `"oneway"` (individual fixed effects only) or `"twoway"`
#'   (individual + time fixed effects).
#' @return A list with `beta_lsdv`, `beta_Mod`, `beta_FM` (length-`q`
#'   coefficient vectors), `VCV_Mod`, `VCV_FM`, `VCV_FM_std` (`q x q`
#'   variance-covariance matrices), and `unit_info` (per-unit diagnostics:
#'   `Omega_i` etc.).
#' @keywords internal
fit_pooled_panel_cpr <- function(y, x, q, kernel, bandwidth, effects = "oneway") {
  theory <- .pooled_theory_matrices(q)
  M <- theory$M
  Q <- theory$Q

  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn <- nrow(x)
  N <- ncol(x)

  dm <- simpledemean(y, x, effects, q)
  yv <- as.numeric(dm$ytilde)
  Xtilde <- if (q == 2) {
    cbind(as.numeric(dm$xtilde), as.numeric(dm$xquad_tilde))
  } else {
    cbind(as.numeric(dm$xtilde), as.numeric(dm$xquad_tilde), as.numeric(dm$xcub_tilde))
  }

  XXtilde <- crossprod(Xtilde)
  invXXtilde <- solve(XXtilde)
  beta_lsdv <- as.numeric(invXXtilde %*% crossprod(Xtilde, yv))
  u_hat <- yv - Xtilde %*% beta_lsdv

  # v_{i,1} = x_{i,1} (as if x_{i,0} = 0): no truncation, keeps the full T.
  vt <- rbind(x[1, , drop = FALSE], diff(x))

  GT <- diag(Tn^(-seq(1, by = 0.5, length.out = q)), q)

  Eqn <- vector("list", N)
  Sum_Lr <- matrix(0, 2, 2)
  Sum_Dr <- matrix(0, 2, 2)
  Sum_Mi <- numeric(q)
  Sum_DMD <- matrix(0, q, q)
  Sum_Omega_udotv <- 0
  Sum_ODMD <- matrix(0, q, q)
  Sum_ODQD <- matrix(0, q, q)
  Sum_Sigma_13 <- matrix(0, q, q)
  Sum_OuuOvv <- 0
  Sum_OudotvOvv <- 0

  for (i in seq_len(N)) {
    idx <- ((i - 1) * Tn + 1):(i * Tn)
    u_hat_i <- u_hat[idx]
    vt_i <- vt[, i]
    X_tilde_i <- Xtilde[idx, , drop = FALSE]

    lv <- estimate_lr_var(cbind(u_hat_i, vt_i), kernel, bandwidth, demean = FALSE)
    Lr_i <- lv$Omega
    Dr_i <- lv$Delta
    Ovv_i <- Lr_i[2, 2]

    Sum_Lr <- Sum_Lr + Lr_i
    Sum_Dr <- Sum_Dr + Dr_i

    M_i <- gen_cpr_corr_vec(x[, i], seq_len(q))
    Sum_Mi <- Sum_Mi + M_i

    D_i <- diag(Ovv_i^(seq_len(q) / 2), q)
    Sum_DMD <- Sum_DMD + D_i %*% M %*% D_i

    Omega_udotv_i <- Lr_i[1, 1] - Lr_i[2, 1]^2 / Ovv_i
    Sum_Omega_udotv <- Sum_Omega_udotv + Omega_udotv_i
    Sum_ODMD <- Sum_ODMD + Omega_udotv_i * (D_i %*% M %*% D_i)

    Ouv2vv_i <- Lr_i[1, 2]^2 / Ovv_i
    Sum_ODQD <- Sum_ODQD + Ouv2vv_i * (D_i %*% Q %*% D_i)

    Sigma13_i <- matrix(0, q, q)
    Sigma13_i[1, 1] <- 0.25 * Lr_i[1, 2]^2
    if (q == 3) {
      Sigma13_i[1, 3] <- 0.5 * Ovv_i * Lr_i[1, 2]^2
      Sigma13_i[3, 1] <- Sigma13_i[1, 3]
      Sigma13_i[3, 3] <- Ovv_i^2 * Lr_i[1, 2]^2
    }
    Sum_Sigma_13 <- Sum_Sigma_13 + Sigma13_i

    if (effects == "twoway") {
      Sum_OuuOvv <- Sum_OuuOvv + Lr_i[1, 1] * Ovv_i
      Sum_OudotvOvv <- Sum_OudotvOvv + Omega_udotv_i * Ovv_i
    }

    Eqn[[i]] <- list(x = x[, i], u_hat = u_hat_i, vt = vt_i, X_tilde = X_tilde_i,
                      Omega = Lr_i, Delta = Dr_i, M = M_i, Omega_udotv = Omega_udotv_i)
  }

  Lr_mean <- Sum_Lr / N
  Dr_mean <- Sum_Dr / N

  # Modified OLS bias correction:
  base_vec <- numeric(q)
  base_vec[1] <- -0.5 * Tn * Lr_mean[1, 2]
  if (q == 3) base_vec[3] <- -(Tn^2) * Lr_mean[2, 2] * Lr_mean[1, 2]
  Sum_Ci_tilde_star <- Dr_mean[2, 1] * Sum_Mi + N * base_vec
  beta_Mod <- as.numeric(beta_lsdv - invXXtilde %*% Sum_Ci_tilde_star)

  # Fully Modified correction:
  Dr_vu_plus <- Dr_mean[2, 1] - Dr_mean[2, 2] / Lr_mean[2, 2] * Lr_mean[2, 1]
  Sum_FM_cor <- numeric(q)
  for (i in seq_len(N)) {
    C_plus_i <- Dr_vu_plus * Eqn[[i]]$M
    y_tildeplus_i <- dm$ytilde[, i] - Lr_mean[1, 2] / Lr_mean[2, 2] * Eqn[[i]]$vt
    FM_cor_i <- as.numeric(crossprod(Eqn[[i]]$X_tilde, y_tildeplus_i)) - C_plus_i
    Sum_FM_cor <- Sum_FM_cor + FM_cor_i
  }
  beta_FM <- as.numeric(invXXtilde %*% Sum_FM_cor)

  Sigma11 <- Sum_ODMD / N
  Sigma12 <- Sum_ODQD / N
  Sigma13 <- Sum_Sigma_13 / N
  Sigma1 <- Sigma11 + Sigma12 - Sigma13
  V1_hat <- Sum_DMD / N
  Omega_udotv_mean <- Sum_Omega_udotv / N

  if (effects == "oneway") {
    invV_hat <- solve(V1_hat)
    VCV_Mod <- (1 / N) * GT %*% (invV_hat %*% Sigma1 %*% invV_hat) %*% GT
    VCV_FM <- (1 / N) * GT %*% (invV_hat %*% Sigma11 %*% invV_hat) %*% GT
  } else {
    OuuOvv_mean <- Sum_OuuOvv / N
    OudotvOvv_mean <- Sum_OudotvOvv / N

    adj2 <- matrix(0, q, q); adj2[2, 2] <- Lr_mean[2, 2]^2 / 12
    V2_hat <- V1_hat - adj2
    invV_hat <- solve(V2_hat)

    adjA <- matrix(0, q, q); adjA[2, 2] <- OuuOvv_mean * Lr_mean[2, 2] / 6
    adjB <- matrix(0, q, q); adjB[2, 2] <- Lr_mean[1, 1] * Lr_mean[2, 2]^2 / 12
    Sigma2 <- Sigma1 - adjA + adjB
    VCV_Mod <- (1 / N) * GT %*% (invV_hat %*% Sigma2 %*% invV_hat) %*% GT

    adjC <- matrix(0, q, q); adjC[2, 2] <- Lr_mean[2, 2] * OudotvOvv_mean / 6
    adjD <- matrix(0, q, q); adjD[2, 2] <- Omega_udotv_mean * Lr_mean[2, 2]^2 / 12
    Sigma2plus <- Sigma11 - adjC + adjD
    VCV_FM <- (1 / N) * GT %*% (invV_hat %*% Sigma2plus %*% invV_hat) %*% GT
  }

  VCV_FM_std <- Omega_udotv_mean * invXXtilde

  list(beta_lsdv = beta_lsdv, beta_Mod = beta_Mod, beta_FM = beta_FM,
       VCV_Mod = VCV_Mod, VCV_FM = VCV_FM, VCV_FM_std = VCV_FM_std,
       unit_info = Eqn, n_units = N, n_time = Tn, q = q, effects = effects)
}

#' @keywords internal
fit_pmg_pcpr <- function(y_list, x_list, orders, w_list, deter_list,
                          estimator, bandwidth, kernel, unit_names,
                          effects = "oneway") {
  if (!identical(toupper(estimator), "FMOLS")) {
    stop("The pooled panel ('pmg') estimator currently only implements 'FMOLS' ",
         "(which provides both a Modified-OLS and a Fully-Modified-OLS point estimate).",
         call. = FALSE)
  }
  if (!is.null(w_list) && any(!vapply(w_list, is.null, logical(1)))) {
    stop("The pooled panel ('pmg') estimator does not support stationary regressors ",
         "(`w`) yet.", call. = FALSE)
  }
  m <- ncol(x_list[[1]])
  if (m != 1) {
    stop("The pooled panel ('pmg') estimator only supports a single integrated ",
         "regressor (ncol(x) == 1).", call. = FALSE)
  }
  if (!(is.numeric(orders) && length(orders) == 1 && orders %in% c(2, 3))) {
    stop("The pooled panel ('pmg') estimator requires `orders` to be a single ",
         "integer, 2 or 3 (the theoretical bias-correction matrices are only ",
         "tabulated for those cases).", call. = FALSE)
  }

  Y <- do.call(cbind, y_list)
  X <- do.call(cbind, lapply(x_list, function(xx) xx[, 1]))
  colnames(Y) <- colnames(X) <- unit_names

  fit <- fit_pooled_panel_cpr(Y, X, q = orders, kernel = kernel, bandwidth = bandwidth,
                               effects = effects)

  coef_names <- paste0("x1^", seq_len(orders))
  coefficients <- fit$beta_FM
  names(coefficients) <- coef_names
  se <- sqrt(diag(fit$VCV_FM))
  names(se) <- coef_names

  list(fits = fit, unit_coefficients = NULL, coefficients = coefficients, se = se)
}
