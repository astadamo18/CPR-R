# Fully Modified OLS estimator for the cointegrating polynomial regression.
# Direct port of FM_CPR.m (Hong & Wagner FM-OLS for CPR).
#
#   y_t = gamma'w_t + delta'deter_t + beta'X_t + u_t
#
# where X_t collects the requested powers of the integrated regressors x_t.

#' Resolve a bandwidth choice to a numeric value
#' @keywords internal
resolve_bandwidth <- function(u, kernel, bandwidth) {
  if (is.numeric(bandwidth)) return(bandwidth)
  switch(bandwidth,
    And91 = bw_andrews91(u, kernel),
    NW    = bw_newey_west(u, kernel, intercept = FALSE),
    stop("Unknown bandwidth selection '", bandwidth,
         "'. Use 'And91', 'AM92', 'NW', or a numeric value.", call. = FALSE)
  )
}

#' Long-run variance estimate for a given bandwidth/kernel choice
#'
#' Dispatches to VAR-pre-whitened Andrews & Monahan (1992) estimation when
#' `bandwidth = "AM92"`, otherwise resolves a scalar bandwidth (fixed,
#' Andrews 1991, or Newey-West) and calls [lr_var()] directly.
#' @keywords internal
estimate_lr_var <- function(u, kernel, bandwidth, demean = FALSE) {
  if (identical(bandwidth, "AM92")) {
    return(lr_var_prewhitened(u, kernel, pw_lag = 1, stab = TRUE, demean = demean))
  }
  bw <- resolve_bandwidth(u, kernel, bandwidth)
  lr_var(u, kernel, bw, demean)
}

#' Fully Modified OLS estimation of a cointegrating polynomial regression
#'
#' @param y Dependent variable, `T x 1`.
#' @param x Integrated (I(1)) regressors, `T x m`.
#' @param orders Powers of `x` to include; see [gen_var_poly_terms()].
#' @param w Optional stationary (I(0)) regressors, `T x s`, or `NULL`.
#' @param deter Deterministic regressors (including the constant, if any),
#'   `T x q`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param bandwidth Bandwidth selection: `"And91"`, `"AM92"`, `"NW"`, or a
#'   fixed numeric value.
#' @return A list of estimation results (coefficients, standard errors,
#'   t-values, residuals, variance-covariance matrices, ...); see
#'   `vignette` / package documentation for details. Intended to be called
#'   through [cpr()], not directly.
#' @keywords internal
fit_fmols_cpr <- function(y, x, orders, w = NULL, deter, kernel, bandwidth,
                           n_lag = NULL, n_lead = NULL) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn0 <- nrow(x)

  # Delta(x_t): computed from the untruncated series before dropping obs. 1.
  v <- diff(x)
  y <- y[2:Tn0, , drop = FALSE]
  x <- x[2:Tn0, , drop = FALSE]
  if (!is.null(w)) w <- as.matrix(w)[2:Tn0, , drop = FALSE]
  deter <- as.matrix(deter)[2:Tn0, , drop = FALSE]

  poly <- gen_var_poly_terms(x, orders, stochastic = TRUE)
  X <- poly$X
  P <- poly$P
  Mstar <- poly$Mstar

  Tn <- nrow(x)
  m <- ncol(x)
  kw <- if (is.null(w)) 0L else ncol(w)
  kd <- ncol(deter)

  idx_gamma <- seq_len(kw)
  idx_delta <- kw + seq_len(kd)
  idx_beta <- kw + kd + seq_len(ncol(X))

  J <- cbind(deter, X)
  Z <- if (kw > 0) cbind(w, J) else J
  ZZinv <- solve(crossprod(Z))

  # OLS
  b_ols <- ZZinv %*% crossprod(Z, y)
  u_ols <- as.numeric(y - Z %*% b_ols)

  if (kw > 0) iww <- solve(crossprod(w))

  # (1) Long-run variance estimators for [u_ols, Delta(x)]. Delta(x) is
  # always demeaned here (Wagner & Reichold 2023, Remark 5) -- this is not a
  # user-controlled option.
  v_dm <- sweep(v, 2, colMeans(v), "-")

  lv <- estimate_lr_var(cbind(u_ols, v_dm), kernel, bandwidth, demean = FALSE)
  Lr <- lv$Omega
  Dr <- lv$Delta

  Lr_vv <- Lr[-1, -1, drop = FALSE]
  Lr_vu <- Lr[-1, 1]
  Lr_vvvu <- solve(Lr_vv, Lr_vu)
  Omega_udotv <- as.numeric(Lr[1, 1] - Lr[1, -1, drop = FALSE] %*% Lr_vvvu)

  Dr_vu <- Dr[-1, 1]
  Dr_vv <- Dr[-1, -1, drop = FALSE]
  Lambda0 <- as.numeric(Dr_vu - Dr_vv %*% Lr_vvvu)

  # (2) Correction terms
  for (i in seq_len(m)) {
    ind <- (P[i] + 1):P[i + 1]
    Mstar[ind] <- Lambda0[i] * Mstar[ind]
  }

  Astar <- numeric(ncol(Z))
  if (kw > 0) {
    Astar[idx_gamma] <- (1 / Tn) * crossprod(w, u_ols) - (1 / Tn) * crossprod(w, v) %*% Lr_vvvu
  }
  Astar[idx_beta] <- Mstar

  # (3) FM estimator
  yplus <- y - v %*% Lr_vvvu
  bplus <- ZZinv %*% (crossprod(Z, yplus) - Astar)

  gamma_fm <- as.numeric(bplus[idx_gamma, , drop = FALSE])
  delta_fm <- as.numeric(bplus[idx_delta, , drop = FALSE])
  beta_fm <- as.numeric(bplus[idx_beta, , drop = FALSE])

  u_plus <- as.numeric(yplus - Z %*% bplus)
  fitted <- as.numeric(Z %*% bplus)

  # Inference for coefficients on the stationary regressors (HAC-type)
  if (kw > 0) {
    S <- w * matrix(u_plus, Tn, kw)
    SLr <- estimate_lr_var(S, kernel, bandwidth, demean = FALSE)$Omega
    varmat0 <- Tn * iww %*% SLr %*% iww
    se_gamma <- sqrt(diag(varmat0))
    t_gamma <- gamma_fm / se_gamma
  } else {
    varmat0 <- NULL
    se_gamma <- numeric(0)
    t_gamma <- numeric(0)
  }

  # Inference for coefficients on deterministic + integrated regressors
  Omega_udotv1 <- estimate_lr_var(u_plus, kernel, bandwidth, demean = FALSE)$Omega
  Omega_udotv1 <- as.numeric(Omega_udotv1)

  JJinv <- solve(crossprod(J))
  varmat <- Omega_udotv * JJinv
  varmat1 <- Omega_udotv1 * JJinv

  se_all <- sqrt(diag(varmat))
  se_delta <- se_all[seq_len(kd)]
  se_beta <- se_all[(kd + 1):length(se_all)]
  t_delta <- delta_fm / se_delta
  t_beta <- beta_fm / se_beta

  # "Naive" OLS-type HAC covariance matrix, for comparison only (ignores
  # cointegration -- inference on beta/delta should use varmat, not this).
  S_ols <- Z * matrix(u_ols, Tn, ncol(Z))
  SLr_ols <- estimate_lr_var(S_ols, kernel, bandwidth, demean = FALSE)$Omega
  varmatOLS <- Tn * ZZinv %*% SLr_ols %*% ZZinv

  list(
    coef_gamma = gamma_fm, coef_delta = delta_fm, coef_beta = beta_fm,
    se_gamma = se_gamma, se_delta = se_delta, se_beta = se_beta,
    t_gamma = t_gamma, t_delta = t_delta, t_beta = t_beta,
    coef_gamma_ols = as.numeric(b_ols[idx_gamma, , drop = FALSE]),
    coef_delta_ols = as.numeric(b_ols[idx_delta, , drop = FALSE]),
    coef_beta_ols = as.numeric(b_ols[idx_beta, , drop = FALSE]),
    fitted = fitted, residuals = u_plus, residuals_ols = u_ols,
    Astar = Astar, Lambda0 = Lambda0,
    Omega_udotv = Omega_udotv, Omega_udotv1 = Omega_udotv1,
    varmat = varmat, varmat1 = varmat1, varmat0 = varmat0, varmatOLS = varmatOLS,
    ZZinv = ZZinv,
    n_obs = Tn, kw = kw, kd = kd, m = m, P = P, powers = poly$powers
  )
}
