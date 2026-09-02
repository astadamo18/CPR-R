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
  Tmax <- Tn - 1L

  # MATLAB's `j+1:end` silently becomes an empty range when j+1 exceeds the
  # array length (its colon operator returns empty for reversed bounds); R's
  # `:` instead produces a *descending* sequence, which crashes downstream
  # indexing in lr_var(). A bandwidth from And91/NW can easily exceed T-1 for
  # a short, persistent series (e.g. T ~ 30). Since lags >= T-1 have no data
  # to form an autocovariance from anyway, capping `upper` at T-1 here is the
  # mathematically correct fix (equivalent to MATLAB's empty-range no-op),
  # not just a defensive patch.
  upper <- if (kernel == "tr") {
    up <- min(M, Tmax)
    if (up >= 1) w[seq_len(up)] <- 1
    up
  } else if (kernel == "ba") {
    up <- min(ceiling(M) - 1, Tmax)
    if (up >= 1) {
      j <- seq_len(up)
      w[j] <- 1 - j / M
    }
    up
  } else if (kernel == "pa") {
    j1_end <- min(floor(M / 2), Tmax)
    if (j1_end >= 1) {
      j <- seq_len(j1_end)
      jj <- j / M
      w[j] <- 1 - 6 * jj^2 + 6 * jj^3
    }
    j2_start <- floor(M / 2) + 1
    j2_end <- min(floor(M), Tmax)
    if (j2_end >= j2_start) {
      j <- j2_start:j2_end
      jj <- j / M
      w[j] <- 2 * (1 - jj)^3
    }
    min(ceiling(M) - 1, Tmax)
  } else if (kernel == "bo") {
    up <- min(ceiling(M) - 1, Tmax)
    if (up >= 1) {
      j <- seq_len(up)
      jj <- j / M
      w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj) / pi
    }
    up
  } else if (kernel == "da") {
    up <- Tmax
    j <- seq_len(up)
    jj <- pi * j / M
    w[j] <- sin(jj) / jj
    up
  } else if (kernel == "qs") {
    up <- Tmax
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
