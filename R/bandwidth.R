# Automatic bandwidth selection for long-run variance estimation.
# Ports of And_HAC91.m (Andrews, 1991) and bwNW.m (Newey & West, 1994).

#' Andrews (1991) automatic bandwidth selection
#'
#' AR(1)-based automatic bandwidth selection of Andrews (1991), applied
#' coordinate-wise to `v`.
#'
#' @param v Matrix of series, `T x dim(v)`.
#' @param kernel One of `"tr"`, `"ba"`, `"pa"`, `"th"` (Tukey-Hanning),
#'   `"qs"`.
#' @return Selected bandwidth (rounded up to the next integer).
#' @keywords internal
bw_andrews91 <- function(v, kernel) {
  v <- as.matrix(v)
  Tn <- nrow(v)
  dimv <- ncol(v)

  rhovec <- numeric(dimv)
  sigma2vec <- numeric(dimv)
  for (j in seq_len(dimv)) {
    y <- v[2:Tn, j]
    x <- v[1:(Tn - 1), j]
    rho <- sum(x * y) / sum(x * x)
    rhovec[j] <- rho
    sigma2vec[j] <- sum((y - x * rho)^2) / Tn
  }

  denom <- sum(sigma2vec^2 / (1 - rhovec)^4)
  numer2 <- sum(4 * rhovec^2 * sigma2vec^2 / (1 - rhovec)^8)
  a2 <- numer2 / denom
  numer1 <- sum(4 * rhovec^2 * sigma2vec^2 / ((1 - rhovec)^6 * (1 + rhovec)^2))
  a1 <- numer1 / denom

  bandwidth <- switch(kernel,
    tr = 0.6611 * (a2 * Tn)^(1 / 5),
    ba = 1.1447 * (a1 * Tn)^(1 / 3),
    pa = 2.6614 * (a2 * Tn)^(1 / 5),
    th = 1.7462 * (a2 * Tn)^(1 / 5),
    qs = 1.3221 * (a2 * Tn)^(1 / 5),
    stop("Andrews (1991) bandwidth selection is not defined for kernel '", kernel,
         "' (supported: 'tr', 'ba', 'pa', 'th', 'qs').", call. = FALSE)
  )

  ceiling(bandwidth)
}

#' Newey-West (1994) automatic bandwidth selection
#'
#' @param v Matrix of series, `T x m`.
#' @param kernel One of `"ba"`, `"pa"`, `"qs"`.
#' @param intercept If `TRUE`, the first column of `v` is given zero weight.
#' @param weights Optional weight vector of length `ncol(v)`; defaults to a
#'   vector of ones (with the first entry set to zero when `intercept =
#'   TRUE`).
#' @return Selected (real-valued) bandwidth.
#' @keywords internal
bw_newey_west <- function(v, kernel, intercept = FALSE, weights = NULL) {
  v <- as.matrix(v)
  Tn <- nrow(v)
  m <- ncol(v)
  vmat <- t(v)

  if (is.null(weights)) {
    weights <- rep(1, m)
    if (intercept) weights[1] <- 0
  }

  npower <- switch(kernel,
    ba = 2 / 9,
    pa = 4 / 25,
    qs = 2 / 25,
    stop("Newey-West bandwidth selection is only defined for kernel 'ba', 'pa', or 'qs'.", call. = FALSE)
  )
  n <- floor(4 * (Tn / 100)^npower)

  vmatw <- as.numeric(weights %*% vmat)

  sigma <- numeric(n + 1)
  for (j in seq_len(n + 1)) {
    sigma[j] <- sum(vmatw[j:Tn] * vmatw[1:(Tn - j + 1)]) / Tn
  }

  s0 <- sigma[1] + 2 * sum(sigma[-1])
  s1 <- 2 * sum(seq_len(n) * sigma[-1])
  s2 <- 2 * sum(seq_len(n)^2 * sigma[-1])

  q <- if (kernel == "ba") 1 else 2
  Tpower <- 1 / (2 * q + 1)

  gamma <- switch(kernel,
    ba = 1.1447 * ((s1 / s0)^2)^Tpower,
    pa = 2.6614 * ((s2 / s0)^2)^Tpower,
    qs = 1.3221 * ((s2 / s0)^2)^Tpower
  )

  gamma * Tn^Tpower
}
