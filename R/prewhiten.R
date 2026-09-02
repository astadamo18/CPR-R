# VAR pre-whitening and Andrews & Monahan (1992) long-run variance.
# Ports of var_m.m, AndMon_Stab.m, AndMon_HAC92.m.

#' Lag a matrix by `n` periods
#'
#' @param x Matrix (or vector).
#' @param n Lag order.
#' @param v Fill value for the initial `n` rows.
#' @keywords internal
lag_matrix <- function(x, n = 1, v = 0) {
  x <- as.matrix(x)
  if (n < 1) return(matrix(numeric(0), 0, ncol(x)))
  Tn <- nrow(x)
  top <- matrix(v, n, ncol(x))
  rbind(top, x[1:(Tn - n), , drop = FALSE])
}

#' Estimate a VAR(p) by OLS (no deterministic terms)
#'
#' Normalization: `y_t = a1*y_{t-1} + ... + ap*y_{t-p} + e_t`.
#'
#' @param y Matrix of observations, `T x s`.
#' @param p VAR order.
#' @return A list with `coeffs` (`s x (s*p)`), `resids` (`(T-p) x s`), `aic`,
#'   `bic`.
#' @keywords internal
var_m <- function(y, p) {
  y <- as.matrix(y)
  Tn <- nrow(y)
  s <- ncol(y)

  regs <- matrix(0, Tn, p * s)
  for (i in seq_len(p)) {
    regs[, ((i - 1) * s + 1):(i * s)] <- lag_matrix(y, i)
  }

  y_eff <- y[(p + 1):Tn, , drop = FALSE]
  r_eff <- regs[(p + 1):Tn, , drop = FALSE]

  coeffs <- solve(crossprod(r_eff), crossprod(r_eff, y_eff))
  resids <- y_eff - r_eff %*% coeffs

  VCV <- crossprod(resids)
  rVCV <- VCV / (Tn - p * s)

  aic <- log(det(rVCV)) + (2 * p * s * s) / Tn
  bic <- log(det(rVCV)) + (p * s * s * log(Tn)) / Tn

  list(coeffs = t(coeffs), resids = resids, aic = aic, bic = bic)
}

#' Eigenvalue-stabilized AR(1)-polynomial estimate
#'
#' Andrews & Monahan (1992) shrink singular values of the estimated AR
#' polynomial (evaluated at 1) that exceed 0.97 in absolute value back to
#' 0.97. Implemented here via SVD, which is numerically equivalent to (and
#' more stable than) the original eigendecomposition-based construction.
#'
#' @param A Estimated AR polynomial matrix, evaluated at `z = 1`.
#' @param thresh Stabilization threshold (default `0.97`).
#' @keywords internal
stabilize_eigen <- function(A, thresh = 0.97) {
  s <- svd(A)
  d <- pmin(s$d, thresh)
  s$u %*% diag(d, nrow = length(d)) %*% t(s$v)
}

#' VAR-pre-whitened long-run variance estimation (Andrews & Monahan, 1992)
#'
#' Combines VAR pre-whitening with Andrews (1991) long-run variance
#' estimation.
#'
#' @param u Residual matrix, `T x m`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param pw_lag Pre-whitening VAR order (1, 2, ...).
#' @param stab Stabilize the eigenvalues of the AR(1) polynomial at 0.97.
#' @param demean Demean `u` (column-wise) before pre-whitening.
#' @return A list with `Omega`, `Delta`, `Sigma`.
#' @keywords internal
lr_var_prewhitened <- function(u, kernel, pw_lag = 1, stab = TRUE, demean = FALSE) {
  u <- as.matrix(u)
  Tn <- nrow(u)
  m <- ncol(u)

  if (demean) u <- sweep(u, 2, colMeans(u), "-")

  vm <- var_m(u, pw_lag)
  coeffs <- vm$coeffs
  resids <- vm$resids

  if (pw_lag == 1) {
    a1 <- diag(m) - coeffs
  } else {
    coeff_ext <- cbind(diag(m), coeffs)
    mult1 <- kronecker(matrix(-1, pw_lag, 1), diag(m))
    multmat <- rbind(diag(m), mult1)
    a1 <- coeff_ext %*% multmat
  }

  inva1 <- solve(a1)
  if (stab) inva1 <- stabilize_eigen(inva1)

  band_pw <- bw_andrews91(resids, kernel)

  Sigma <- crossprod(u) / Tn

  lv <- lr_var(resids, kernel, band_pw, demean)

  Omega <- inva1 %*% lv$Omega %*% t(inva1)

  Lambda_pw <- lv$Delta - lv$Sigma
  Delta <- Sigma + inva1 %*% Lambda_pw %*% t(inva1) + inva1 %*% coeffs %*% Sigma

  list(Omega = Omega, Delta = Delta, Sigma = Sigma)
}
