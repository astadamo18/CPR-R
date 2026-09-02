# Polynomial regressor and FM-OLS correction-term construction.
# Ports of GenPowerReg.m, GenCPRCorrVec.m, GenVarPolyTerms.m.
#
# Simplification relative to the MATLAB source: the original GenPowerReg /
# GenCPRCorrVec each had an `all` ("yes"/"no") switch between "powers 1:p"
# and "explicit power vector". Passing powvec = 1:p to the explicit branch
# gives an identical result to the "yes" branch in both cases, so here there
# is a single explicit-power-vector code path.

#' Powers of a single regressor
#'
#' @param x Numeric vector, the regressor.
#' @param powvec Integer vector of powers to compute.
#' @return A `length(x) x length(powvec)` matrix.
#' @keywords internal
gen_power_reg <- function(x, powvec) {
  x <- as.numeric(x)
  powvec <- as.numeric(powvec)
  out <- vapply(powvec, function(p) x^p, numeric(length(x)))
  matrix(out, nrow = length(x), ncol = length(powvec))
}

#' FM-OLS correction terms for powers of an integrated regressor
#'
#' For power `p`, the correction term is `p * sum(x_t^(p-1))`.
#'
#' @param x Numeric vector, the (integrated) regressor.
#' @param powvec Integer vector of powers.
#' @return Numeric vector of correction terms, one per entry of `powvec`.
#' @keywords internal
gen_cpr_corr_vec <- function(x, powvec) {
  x <- as.numeric(x)
  Tn <- length(x)
  powvec <- as.numeric(powvec)
  if (length(powvec) == 0) return(numeric(0))
  maxpower <- max(powvec)

  sum_matrix <- matrix(1, Tn, maxpower)
  if (maxpower >= 2) sum_matrix[, 2] <- x
  if (maxpower >= 3) {
    for (j in 3:maxpower) sum_matrix[, j] <- sum_matrix[, j - 1] * x
  }

  sum_vec <- colSums(sum_matrix)
  full_vec <- seq_len(maxpower) * sum_vec
  full_vec[powvec]
}

#' Build polynomial regressors (and FM-OLS correction terms) for integrated
#' regressors
#'
#' @param x Matrix of (integrated) regressors, `T x m`.
#' @param orders Powers to include. Either a single integer (same max order
#'   `1:order` for every column), a numeric vector of length `ncol(x)`
#'   (per-column max order `1:order_i`), or a list of length `ncol(x)` giving
#'   explicit (possibly non-consecutive) powers per column.
#' @param stochastic If `TRUE` (default), also compute the FM-OLS correction
#'   terms `Mstar`. Set to `FALSE` for purely deterministic polynomial terms.
#' @return A list with `X` (the polynomial regressor matrix), `P` (index
#'   vector: column `P[i]+1` to `P[i+1]` of `X` belongs to regressor `i`),
#'   `Mstar` (correction terms, or `NULL` if `stochastic = FALSE`), and
#'   `powers` (the list of power vectors actually used, one per column of
#'   `x`, for labeling purposes).
#' @keywords internal
gen_var_poly_terms <- function(x, orders, stochastic = TRUE) {
  x <- as.matrix(x)
  m <- ncol(x)

  if (is.list(orders)) {
    if (length(orders) != m) {
      stop("`orders` as a list must have length ncol(x).", call. = FALSE)
    }
    powers <- lapply(orders, as.numeric)
  } else if (length(orders) == 1) {
    powers <- replicate(m, seq_len(orders), simplify = FALSE)
  } else if (length(orders) == m) {
    powers <- lapply(orders, seq_len)
  } else {
    stop("`orders` must be a single integer, a numeric vector of length ncol(x), ",
         "or a list of length ncol(x).", call. = FALSE)
  }

  X <- NULL
  Mstar <- if (stochastic) numeric(0) else NULL
  P <- 0
  for (i in seq_len(m)) {
    pw <- powers[[i]]
    X <- cbind(X, gen_power_reg(x[, i], pw))
    if (stochastic) Mstar <- c(Mstar, gen_cpr_corr_vec(x[, i], pw))
    P <- c(P, P[length(P)] + length(pw))
  }

  list(X = X, P = P, Mstar = Mstar, powers = powers)
}
