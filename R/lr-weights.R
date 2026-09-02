# Kernel weights for long-run variance estimation.
# Direct port of lr_weights.m (Wagner et al., MATLAB CTPUTests toolbox).

#' Kernel weights for long-run variance estimation
#'
#' Computes the vector of kernel weights at lags 1, ..., T-1 used to weight
#' sample autocovariances when building a long-run variance estimate.
#'
#' @param Tn Number of time-series observations.
#' @param kernel One of `"tr"` (truncated), `"ba"` (Bartlett), `"pa"`
#'   (Parzen), `"bo"` (Bohman), `"da"` (Daniell), `"qs"` (Quadratic
#'   Spectral).
#' @param band Bandwidth (need not be an integer, except for `"tr"`).
#' @return A list with `w` (weight vector of length `Tn - 1`) and `upper`
#'   (index of the largest lag with a possibly non-zero weight).
#' @keywords internal
lr_weights <- function(Tn, kernel, band) {
  w <- numeric(Tn - 1)
  M <- band

  upper <- if (kernel == "tr") {
    up <- min(M, Tn - 1)
    if (up >= 1) w[seq_len(up)] <- 1
    up
  } else if (kernel == "ba") {
    up <- ceiling(M) - 1
    if (up >= 1) {
      j <- seq_len(up)
      w[j] <- 1 - j / M
    }
    up
  } else if (kernel == "pa") {
    j1_end <- floor(M / 2)
    if (j1_end >= 1) {
      j <- seq_len(j1_end)
      jj <- j / M
      w[j] <- 1 - 6 * jj^2 + 6 * jj^3
    }
    j2_start <- j1_end + 1
    j2_end <- floor(M)
    if (j2_end >= j2_start) {
      j <- j2_start:j2_end
      jj <- j / M
      w[j] <- 2 * (1 - jj)^3
    }
    ceiling(M) - 1
  } else if (kernel == "bo") {
    up <- ceiling(M) - 1
    if (up >= 1) {
      j <- seq_len(up)
      jj <- j / M
      w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj) / pi
    }
    up
  } else if (kernel == "da") {
    up <- Tn - 1
    j <- seq_len(up)
    jj <- pi * j / M
    w[j] <- sin(jj) / jj
    up
  } else if (kernel == "qs") {
    up <- Tn - 1
    sc <- 6 * pi / 5
    j <- seq_len(up)
    jj <- j / M
    w[j] <- 25 / (12 * pi^2 * jj^2) * (sin(sc * jj) / (sc * jj) - cos(sc * jj))
    up
  } else {
    stop("Unknown kernel '", kernel, "'. Use one of 'tr', 'ba', 'pa', 'bo', 'da', 'qs'.", call. = FALSE)
  }

  list(w = w, upper = max(upper, 0))
}
