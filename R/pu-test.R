# Phillips-Ouliaris-type PU cointegration test for a cointegrating
# polynomial regression. Port of PU_test.m.
#
# Unlike ct_test() (which needs FM-OLS residuals and Omega_udotv1 from a
# fitted cpr() object), PU_test operates directly on the raw y/x series: it
# runs its own plain-OLS regression to get residuals for the denominator,
# and a separate local VAR(1)-with-deterministics on the stacked [y, x]
# system for the long-run variance in the numerator. So it is not simply
# "the same computation with different critical values" -- it is a
# genuinely different test statistic, reusing the same kernel/bandwidth/
# long-run-variance machinery as cpr() and ct_test().
#
# Critical value tables are extracted from the original CTPUTests MATLAB
# toolbox's PUcritval/*.mat files (same 9-percentile layout as CT's tables).
# Only d = 0, m = 1, p = 2 is bundled so far, matching the CEE example.

.pu_critval_table <- list(
  "d0_m1_p2" = c(
    `1%` = 2.0262593733019623, `2.5%` = 2.6460146295490197,
    `5%` = 3.4310662947470956, `10%` = 4.6943570569576085,
    `50%` = 13.436969556179534, `90%` = 30.982539242290752,
    `95%` = 37.87475517642149, `97.5%` = 44.432484707331454,
    `99%` = 52.94784857075556
  )
)

#' @keywords internal
pu_critval <- function(d, m, p) {
  key <- paste0("d", d, "_m", m, "_p", p)
  tab <- .pu_critval_table[[key]]
  if (is.null(tab)) {
    stop("No PU critical value table bundled for d=", d, ", m=", m, ", p=", p,
         ". Only d=0, m=1, p=2 is available so far; further tables can be added ",
         "to `.pu_critval_table` in R/pu-test.R from the original PUcritval/*.mat files.",
         call. = FALSE)
  }
  tab
}

#' Phillips-Ouliaris-type PU cointegration test for a cointegrating polynomial regression
#'
#' Tests the null hypothesis of *no* cointegration against the alternative
#' of cointegration (the reverse null from [ct_test()]), using the
#' Phillips-Ouliaris/Wagner-type ratio statistic of Wagner and co-authors
#' for cointegrating polynomial regressions. Port of `PU_test.m`.
#'
#' Operates directly on the raw series: an internal plain-OLS regression of
#' `y` on the deterministics, any linearly-entered regressors, and the
#' requested powers of the last column of `x`, plus a separate local
#' VAR(1)-with-deterministics on the stacked `[y, x]` system for the
#' long-run variance. It does not require (or take) a fitted [cpr()]
#' object.
#'
#' @param y Dependent variable, length `T`.
#' @param x Integrated regressors, `T x m`. If `m > 1`, only the *last*
#'   column gets the polynomial powers requested via `orders`; the
#'   remaining `m - 1` columns enter linearly (power 1).
#' @param d Deterministic specification: `-1` (none), `0` (intercept), `1`
#'   (intercept + trend). Only `0` is currently tabulated.
#' @param m Number of integrated regressors (`ncol(x)`). Only `1` is
#'   currently tabulated.
#' @param orders Powers of the last column of `x` to include; see
#'   [gen_var_poly_terms()]. Only a single max power of `2` is currently
#'   tabulated. If `orders` is exactly `1` (linear only), no critical
#'   values are available in the original toolbox either, and `critval`/
#'   `reject` come back as `NA`.
#' @param kernel Kernel function, see [lr_weights()].
#' @param bandwidth Bandwidth selection: `"And91"`, `"AM92"`, `"NW"`, or a
#'   fixed numeric value.
#' @param alpha Significance levels to test at; must be a subset of
#'   `c(0.1, 0.05, 0.01)`.
#'
#' @return A list with `statistic`, `alpha`, `critval`, and `reject`
#'   (logical per `alpha`; `TRUE` means reject the null of no
#'   cointegration, i.e. evidence *for* cointegration).
#' @export
pu_test <- function(y, x, d, m, orders, kernel, bandwidth, alpha = c(0.1, 0.05, 0.01)) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn <- nrow(y)

  const <- rep(1, Tn)
  trend <- seq_len(Tn)
  deter <- switch(as.character(d),
    "-1" = matrix(numeric(0), Tn, 0),
    "0"  = matrix(const, ncol = 1, dimnames = list(NULL, "const")),
    "1"  = cbind(const = const, trend = trend),
    stop("`d` must be -1, 0, or 1.", call. = FALSE)
  )

  z <- cbind(y, x)

  linear_only <- is.numeric(orders) && length(orders) == 1 && orders == 1
  xpower <- x[, m, drop = FALSE]

  if (linear_only) {
    powerreg <- NULL
  } else {
    powerreg <- gen_var_poly_terms(xpower, orders, stochastic = FALSE)$X
  }

  if (linear_only) {
    regmat4uhat <- cbind(deter, x)
  } else {
    x_rest <- if (m > 1) x[, seq_len(m - 1), drop = FALSE] else NULL
    regmat4uhat <- cbind(deter, x_rest, powerreg)
  }

  coeff4uhat <- solve(crossprod(regmat4uhat), crossprod(regmat4uhat, y))
  uhat <- y - regmat4uhat %*% coeff4uhat

  # Local VAR(1)-with-deterministics on the stacked [y, x] system:
  depvar <- z[2:Tn, , drop = FALSE]
  lagz <- lag_matrix(z, 1)[2:Tn, , drop = FALSE]
  indepvar <- cbind(deter[2:Tn, , drop = FALSE], lagz)

  varcoeff <- solve(crossprod(indepvar), crossprod(indepvar, depvar))
  varresid <- depvar - indepvar %*% varcoeff

  Lr <- estimate_lr_var(varresid, kernel, bandwidth, demean = FALSE)$Omega
  omega_udotv <- as.numeric(Lr[1, 1] - Lr[1, -1, drop = FALSE] %*% solve(Lr[-1, -1, drop = FALSE], Lr[-1, 1]))

  statistic <- omega_udotv / (Tn^(-2) * sum(uhat^2))

  if (linear_only) {
    return(list(statistic = statistic, alpha = alpha,
                critval = rep(NA_real_, length(alpha)), reject = rep(NA, length(alpha))))
  }

  crit_table <- pu_critval(d, m, orders)
  pct_label <- c(`0.1` = "90%", `0.05` = "95%", `0.01` = "99%")
  lab <- pct_label[as.character(alpha)]
  if (any(is.na(lab))) {
    stop("`alpha` must be a subset of c(0.1, 0.05, 0.01).", call. = FALSE)
  }
  critval <- as.numeric(crit_table[lab])

  list(statistic = statistic, alpha = alpha, critval = critval, reject = statistic > critval)
}
