# Dynamic OLS estimator for the cointegrating polynomial regression.
# Direct port of DOLS_CPR.m (Schneeberger, based on Wagner's FM_CPR) and its
# GenLeadLag.m helper, both from MonitoringCPR_MatlabCode/MonitoringCPR/.
#
#   y_t = delta'deter_t + beta'X_t + chi'x_LeadLag_t + u_t
#
# where X_t collects the requested powers of the integrated regressors and
# x_LeadLag_t holds leads and lags of Delta(x_t) (levels of x only, not of
# its powers). Consistency is achieved by the lead/lag augmentation itself,
# via plain OLS on the augmented equation; the long-run variance step
# afterward is only needed for standard errors, not for a bias correction
# (unlike FM-OLS, there is no Astar/Mstar correction here).
#
# Notable differences from fit_fmols_cpr(), both faithful to the original
# source rather than smoothed over:
#  - DOLS_CPR.m does not support stationary regressors (`w`) at all.
#  - When n_lag = n_lead = 0, DOLS_CPR.m does NOT drop the first
#    observation (unlike FM_CPR.m's harmonized single-truncation
#    convention) -- it is plain OLS on the polynomial regression using the
#    full T-length sample. Truncation only happens, and only by exactly as
#    much as the lead/lag structure requires, when n_lag + n_lead > 0.
#  - The long-run variance for inference is the (univariate) long-run
#    variance of the DOLS residuals themselves (already "purified" of the
#    endogeneity by the lead/lag terms), not the FM-OLS Schur-complement
#    Omega_u.v built from [u_ols, Delta(x)] jointly.

#' Leads and lags of a (multivariate) series
#'
#' Port of GenLeadLag.m, in the natural T x n (time down rows) orientation.
#' Column blocks are ordered [contemporaneous, lag_1, ..., lag_{n_lag},
#' lead_1, ..., lead_{n_lead}], each block `ncol(v)` columns wide. Rows
#' where a given lag/lead is not available are filled with zero (matching
#' the original, which never trims for this reason -- the caller trims the
#' rows where any block would still be incomplete).
#'
#' @param v Matrix, `T x n`.
#' @param n_lag,n_lead Number of lags / leads to include.
#' @keywords internal
gen_lead_lag <- function(v, n_lag, n_lead) {
  v <- as.matrix(v)
  Tn <- nrow(v)
  n <- ncol(v)

  blocks <- list(v)
  for (j in seq_len(n_lag)) {
    lagblock <- matrix(0, Tn, n)
    if (Tn - j >= 1) lagblock[(j + 1):Tn, ] <- v[1:(Tn - j), ]
    blocks[[length(blocks) + 1]] <- lagblock
  }
  for (j in seq_len(n_lead)) {
    leadblock <- matrix(0, Tn, n)
    if (Tn - j >= 1) leadblock[1:(Tn - j), ] <- v[(j + 1):Tn, ]
    blocks[[length(blocks) + 1]] <- leadblock
  }
  do.call(cbind, blocks)
}

#' Dynamic OLS estimation of a cointegrating polynomial regression
#'
#' @inheritParams fit_fmols_cpr
#' @param n_lag Number of contemporaneous + lagging first differences of
#'   `x` to include (0 = no lags, not even the contemporaneous one).
#' @param n_lead Number of leading first differences of `x` to include.
#' @keywords internal
fit_dols_cpr <- function(y, x, orders, w = NULL, deter, kernel, bandwidth,
                          n_lag = 0, n_lead = 0) {
  if (!is.null(w)) {
    stop("The DOLS estimator does not support stationary regressors (`w`) -- ",
         "this matches the original DOLS_CPR.m, which only takes deterministic ",
         "components and (polynomial) integrated regressors.", call. = FALSE)
  }

  y <- as.matrix(y)
  x <- as.matrix(x)
  deter <- as.matrix(deter)
  Tn <- nrow(x)

  poly <- gen_var_poly_terms(x, orders, stochastic = FALSE)
  X <- poly$X

  kd <- ncol(deter)
  m_aug <- ncol(X)
  idx_delta <- seq_len(kd)
  idx_beta <- kd + seq_len(m_aug)

  Z <- cbind(deter, X)

  b_ols <- solve(crossprod(Z), crossprod(Z, y))
  u_ols <- as.numeric(y - Z %*% b_ols)

  if ((n_lag + n_lead) == 0) {
    Z_trunc <- Z
    y_trunc <- as.numeric(y)
    b_dols <- as.numeric(b_ols)
    u_dols <- u_ols
  } else {
    v <- diff(x)
    X_LeadLag <- gen_lead_lag(v, n_lag, n_lead)

    Z_untrunc <- cbind(Z[2:Tn, , drop = FALSE], X_LeadLag)
    keep <- (n_lag + 1):(nrow(Z_untrunc) - n_lead)
    Z_trunc <- Z_untrunc[keep, , drop = FALSE]
    y_trunc <- as.numeric(y[(n_lag + 2):(Tn - n_lead), , drop = FALSE])

    b_dols <- as.numeric(solve(crossprod(Z_trunc), crossprod(Z_trunc, y_trunc)))
    u_dols <- y_trunc - as.numeric(Z_trunc %*% b_dols)
  }

  # Long-run variance directly from the DOLS residuals (already purified of
  # the endogeneity by the lead/lag augmentation -- no Schur-complement
  # correction needed here, unlike FM-OLS).
  Omega_udotv <- as.numeric(estimate_lr_var(u_dols, kernel, bandwidth, demean = FALSE)$Omega)

  delta_dols <- b_dols[idx_delta]
  beta_dols <- b_dols[idx_beta]

  ZZinv_trunc <- solve(crossprod(Z_trunc))
  varmat_dols <- Omega_udotv * ZZinv_trunc
  se_all <- sqrt(diag(varmat_dols))
  se_delta <- se_all[idx_delta]
  se_beta <- se_all[idx_beta]

  # Fitted values from the structural part only (deter + X), dropping the
  # lead/lag nuisance coefficients, evaluated at the full (untruncated)
  # sample -- matches `result.Fitted = Z*b_dols(1:(kd+m_aug))` in the
  # original.
  fitted <- as.numeric(Z %*% b_dols[seq_len(kd + m_aug)])

  list(
    coef_gamma = numeric(0), coef_delta = delta_dols, coef_beta = beta_dols,
    se_gamma = numeric(0), se_delta = se_delta, se_beta = se_beta,
    t_gamma = numeric(0), t_delta = delta_dols / se_delta, t_beta = beta_dols / se_beta,
    coef_gamma_ols = numeric(0), coef_delta_ols = as.numeric(b_ols[idx_delta]),
    coef_beta_ols = as.numeric(b_ols[idx_beta]),
    fitted = fitted, residuals = u_dols, residuals_ols = u_ols,
    Omega_udotv = Omega_udotv, Omega_udotv1 = Omega_udotv,
    varmat = varmat_dols[c(idx_delta, idx_beta), c(idx_delta, idx_beta), drop = FALSE],
    varmat1 = NULL, varmat0 = NULL, varmatOLS = NULL,
    n_obs = length(y_trunc), kw = 0L, kd = kd, m = ncol(x), P = poly$P, powers = poly$powers
  )
}
