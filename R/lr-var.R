# Long-run variance estimation.
# Direct port of lr_var.m.

#' Kernel-based long-run variance estimation
#'
#' Computes the long-run variance `Omega`, the one-sided long-run variance
#' `Delta` (starting at lag 0), and the contemporaneous variance `Sigma` of
#' an input matrix of (residual) series.
#'
#' @param u Matrix (or vector) of series, `T x m`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param band Bandwidth, integer in `1, ..., T`.
#' @param demean Demean `u` (column-wise) before computing autocovariances.
#' @return A list with `Omega`, `Delta`, `Sigma` (each `m x m` matrices).
#' @keywords internal
lr_var <- function(u, kernel, band, demean = FALSE) {
  u <- as.matrix(u)
  Tn <- nrow(u)
  m <- ncol(u)

  if (demean) u <- sweep(u, 2, colMeans(u), "-")

  lw <- lr_weights(Tn, kernel, band)
  w <- lw$w
  j_max <- lw$upper

  Sigma <- crossprod(u) / Tn
  Omega <- matrix(0, m, m)
  Delta <- matrix(0, m, m)

  if (j_max >= 1) {
    for (j in seq_len(j_max)) {
      T1 <- crossprod(u[(j + 1):Tn, , drop = FALSE], u[1:(Tn - j), , drop = FALSE]) / Tn
      T2 <- crossprod(u[1:(Tn - j), , drop = FALSE], u[(j + 1):Tn, , drop = FALSE]) / Tn
      Omega <- Omega + w[j] * (T1 + T2)
      Delta <- Delta + w[j] * T2
    }
  }

  Omega <- Omega + Sigma
  Delta <- Delta + Sigma

  list(Omega = Omega, Delta = Delta, Sigma = Sigma)
}
