# Integrated Modified OLS (IM-OLS) estimator for the cointegrating
# polynomial regression. Port of IMOLS_NL.m (MonitoringCPR_MatlabCode) --
# its baseline "regression (1)" procedure -- plus Vhat_NEW.m for the
# sandwich covariance matrix.
#
#   y_t = delta'deter_t + beta'X_t + u_t
#
# Unlike FM-OLS, IM-OLS needs no Schur-complement bias correction and no
# truncation of the sample: it partial-sums (cumsum) both sides of the
# equation and adds the *un-summed* regressor levels back in as
# "augmentation" terms to absorb the endogeneity between X_t and u_t, then
# runs plain OLS on that transformed system (Vogelsang & Wagner, 2014).
#
# Vhat_NEW.m's own docstring says its output is the coefficient covariance
# "up to cond. LR-var." (its Lambda^2) -- i.e. a sandwich "meat" that still
# needs multiplying by the conditional long-run variance of u given
# Delta(x), Omega_udotv, to become the actual asymptotic covariance
# matrix. IMOLS_NL.m/Vhat_NEW.m don't spell out how to estimate that
# scaling factor themselves, but the newer, more general im_scmpr.m
# (IM-SCMPR-ExemplaryCode) completes the exact same sandwich formula with
# `coeff_var = kron(lrvar_cond, calM_hat)`, where lrvar_cond is estimated
# from [OLS residuals, Delta(x)]'s joint long-run variance -- precisely
# the same Omega_udotv construction fit_fmols_cpr() already uses. That
# construction is reused here rather than re-derived.
#
# Scope, faithful to what's actually implemented rather than padded out:
#  - Only IMOLS_NL.m's "regression (1)" (its baseline estimator) is
#    ported, not its optional "regression (2)" Z-augmented refinement
#    (selector = 2/3 in the original) -- that refinement's exact purpose
#    isn't spelled out well enough in the source alone to port faithfully.
#  - Only "full augmentation" (augtype = 2: every power of x doubles as
#    its own augmentation regressor) is ported, not "linear augmentation"
#    (augtype = 1: only the base level of x). The source leaves the choice
#    between them to the caller with no stated default; "full" is the one
#    that mirrors FM-OLS's own per-power correction (Mstar/Lambda0) most
#    closely, so it was chosen for consistency across estimators.
#  - Does not support stationary regressors (`w`) at all -- IMOLS_NL.m's
#    own signature never took one.

#' "Up to Lambda^2" sandwich covariance matrix for the IM-OLS regressors
#'
#' Port of Vhat_NEW.m. Still needs multiplying by the conditional
#' long-run variance of the residuals (Omega_udotv) to become the actual
#' asymptotic covariance matrix -- see the file-level comment above.
#'
#' @param g The regressor matrix of the partial-summed-and-augmented
#'   regression (`Xmat` below), `T x k`.
#' @keywords internal
imols_vhat <- function(g) {
  Tn <- nrow(g)
  S <- apply(g, 2, cumsum)
  SJ <- rbind(0, S[-Tn, , drop = FALSE])
  ST <- matrix(S[Tn, ], nrow = Tn, ncol = ncol(g), byrow = TRUE)
  DS <- ST - SJ
  gg <- solve(crossprod(g))
  Vh <- gg %*% t(DS)
  Vh %*% t(Vh)
}

#' Integrated Modified OLS (IM-OLS) estimation of a cointegrating
#' polynomial regression
#'
#' @inheritParams fit_fmols_cpr
#' @keywords internal
fit_imols_cpr <- function(y, x, orders, w = NULL, deter, kernel, bandwidth,
                           n_lag = NULL, n_lead = NULL) {
  if (!is.null(w)) {
    stop("The IMOLS estimator does not support stationary regressors (`w`) -- ",
         "this matches the original IMOLS_NL.m, which only takes deterministic ",
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

  # Plain OLS on levels -- only to get u_ols for the conditional long-run
  # variance step below (the same role fit_fmols_cpr() uses its own u_ols
  # for), not IM-OLS's point estimator itself.
  b_ols <- solve(crossprod(Z), crossprod(Z, y))
  u_ols <- as.numeric(y - Z %*% b_ols)

  v <- diff(x)
  v_dm <- sweep(v, 2, colMeans(v), "-")
  lv <- estimate_lr_var(cbind(u_ols[2:Tn], v_dm), kernel, bandwidth, demean = FALSE)
  Lr <- lv$Omega
  Lr_vv <- Lr[-1, -1, drop = FALSE]
  Lr_vu <- Lr[-1, 1]
  Omega_udotv <- as.numeric(Lr[1, 1] - Lr[1, -1, drop = FALSE] %*% solve(Lr_vv, Lr_vu))

  # IM-OLS point estimator: partial-sum both sides, augment with the
  # un-summed regressor levels, plain OLS on the transformed system.
  SD <- apply(deter, 2, cumsum)
  SX <- apply(X, 2, cumsum)
  Sy <- cumsum(y)
  Xmat <- cbind(SD, SX, X)

  theta_hat <- solve(crossprod(Xmat), crossprod(Xmat, Sy))
  delta_im <- as.numeric(theta_hat[idx_delta, ])
  beta_im <- as.numeric(theta_hat[idx_beta, ])

  fitted <- as.numeric(Z %*% theta_hat[seq_len(kd + m_aug), ])
  residuals <- as.numeric(y) - fitted

  varmat_full <- Omega_udotv * imols_vhat(Xmat)
  varmat <- varmat_full[seq_len(kd + m_aug), seq_len(kd + m_aug), drop = FALSE]
  se_all <- sqrt(diag(varmat))
  se_delta <- se_all[idx_delta]
  se_beta <- se_all[idx_beta]

  list(
    coef_gamma = numeric(0), coef_delta = delta_im, coef_beta = beta_im,
    se_gamma = numeric(0), se_delta = se_delta, se_beta = se_beta,
    t_gamma = numeric(0), t_delta = delta_im / se_delta, t_beta = beta_im / se_beta,
    coef_gamma_ols = numeric(0), coef_delta_ols = as.numeric(b_ols[idx_delta]),
    coef_beta_ols = as.numeric(b_ols[idx_beta]),
    fitted = fitted, residuals = residuals, residuals_ols = u_ols,
    Omega_udotv = Omega_udotv, Omega_udotv1 = Omega_udotv,
    varmat = varmat, varmat1 = NULL, varmat0 = NULL, varmatOLS = NULL,
    n_obs = Tn, kw = 0L, kd = kd, m = ncol(x), P = poly$P, powers = poly$powers
  )
}
