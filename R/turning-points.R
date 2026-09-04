# Turning-point analysis for cointegrating polynomial regressions.
#
# For a fitted CPR y_t = const + beta_1*x_t + beta_2*x_t^2 + ... (the
# stationary regressors `w` and any trend, if present, are deliberately left
# out of this curve -- there is no single canonical value to hold them at),
# a "turning point" is a value x* where d(prediction)/dx = 0: the standard
# EKC-style interpretation of a quadratic/cubic cointegrating polynomial
# relationship (e.g. de Jong & Wagner's income/emissions curve). Turning
# point *location* only depends on the slope coefficients (beta_1, beta_2,
# ...); the curve's *level* additionally needs the constant, which never
# moves where the turning point sits on the x-axis but is still needed to
# plot (and label) the right y-value there -- hence including it explicitly
# everywhere below, even though it "does not matter" for x* itself.

#' Derivative-root turning points of a polynomial in a single regressor
#'
#' @param beta Numeric vector of slope coefficients for the powers of x, in
#'   the same order as `powers` (e.g. `c(beta1, beta2)` for `x^1`, `x^2`).
#' @param powers Integer vector of the powers `beta` corresponds to.
#' @param const Constant (intercept) to add to the curve's level; `0` if
#'   there is none to add.
#' @param x_range Optional length-2 numeric vector; turning points outside
#'   this range are dropped (interior turning points only). `NULL` keeps
#'   every real root.
#' @return A data frame with columns `x`, `y` (curve value at `x`, including
#'   `const`), and `type` (`"maximum"`, `"minimum"`, or `"inflection"`),
#'   sorted by `x`. Zero rows if there is no turning point (e.g. a purely
#'   linear relationship, `max(powers) < 2`).
#' @keywords internal
poly_turning_points <- function(beta, powers, const = 0, x_range = NULL) {
  stopifnot(length(beta) == length(powers))
  none <- data.frame(x = numeric(0), y = numeric(0), type = character(0))
  max_p <- max(powers)
  if (max_p < 2) return(none)

  # Coefficients of the derivative polynomial in x^0, x^1, ..., x^(max_p-1):
  # d/dx(beta_k * x^p_k) = p_k * beta_k * x^(p_k - 1), which lands in slot
  # p_k of this (1-indexed) vector.
  deriv_coef <- numeric(max_p)
  for (k in seq_along(powers)) {
    p <- powers[k]
    deriv_coef[p] <- deriv_coef[p] + p * beta[k]
  }
  if (all(deriv_coef == 0)) return(none)

  roots <- polyroot(deriv_coef)
  is_real <- abs(Im(roots)) < 1e-6 * pmax(1, abs(Re(roots)))
  real_roots <- sort(unique(Re(roots)[is_real]))
  if (length(real_roots) == 0) return(none)

  if (!is.null(x_range)) {
    real_roots <- real_roots[real_roots >= min(x_range) & real_roots <= max(x_range)]
  }
  if (length(real_roots) == 0) return(none)

  curve_value <- function(xv) const + sum(beta * xv^powers)
  second_deriv <- function(xv) {
    sum(ifelse(powers < 2, 0, powers * (powers - 1) * beta * xv^pmax(powers - 2, 0)))
  }

  d2 <- vapply(real_roots, second_deriv, numeric(1))
  type <- ifelse(d2 > 1e-8, "minimum", ifelse(d2 < -1e-8, "maximum", "inflection"))
  y <- vapply(real_roots, curve_value, numeric(1))

  data.frame(x = real_roots, y = y, type = type, row.names = NULL)
}

#' Constant (intercept) coefficient of a fit, or 0 if there is none
#' @keywords internal
get_const_coef <- function(coefficients) {
  if ("const" %in% names(coefficients)) return(unname(coefficients[["const"]]))
  0
}

#' Turning point(s) of a fitted cointegrating polynomial regression
#'
#' Where the fitted curve's slope with respect to the (single) integrated
#' regressor crosses zero -- the EKC-style turning point of a
#' quadratic/cubic cointegrating polynomial relationship.
#'
#' @param object A fitted `"cpr"` or `"pcpr"` object.
#' @param ... Passed on to methods.
#' @return A data frame; see the method-specific documentation.
#' @export
turning_points <- function(object, ...) {
  UseMethod("turning_points")
}

#' @param x_range Restrict to turning points inside this x-range (a length-2
#'   vector); `"data"` (default) uses the observed range of the fitted
#'   regressor, i.e. interior turning points only. Pass `NULL` to keep every
#'   real root, including extrapolated ones outside the observed data.
#' @rdname turning_points
#' @export
turning_points.cpr <- function(object, x_range = "data", ...) {
  if (ncol(object$x) != 1) {
    stop("turning_points() only supports a cpr fit with a single integrated ",
         "regressor (ncol(x) == 1).", call. = FALSE)
  }
  xname <- colnames(object$x)[1]
  powers1 <- object$fit$powers[[1]]
  beta <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const <- get_const_coef(object$coefficients)

  if (identical(x_range, "data")) x_range <- range(object$x[, 1])
  poly_turning_points(beta, powers1, const = const, x_range = x_range)
}

#' @details
#' For `type = "mg"`, the turning point of each per-unit fit is computed
#' first (each restricted to that unit's own observed x-range), then
#' averaged across units *by type* (all `"maximum"` points averaged
#' together, all `"minimum"` points averaged together) -- the mean-group
#' philosophy applied to a nonlinear function of the estimates, the same way
#' [pcpr()]'s coefficients themselves are a mean-group average of per-unit
#' coefficients. The reported `y` at each averaged `x` uses the panel's own
#' group-mean curve (`object$coefficients`, which already averages the
#' constant across units along with the slopes).
#'
#' For `type = "pmg"`, there is a single common slope, so at most one
#' turning point of each type. The pooled model has no single estimated
#' constant (individual, and possibly time, fixed effects absorb it
#' instead), so the constant used for the curve's level is the average,
#' across units, of that unit's own implied fixed effect `alpha_i =
#' mean(y_i) - beta_FM' * mean(x_i^powers)` -- reconstructed from each
#' unit's own *raw* (not demeaned) `y`/`x`, so it is a real, data-scale
#' level, not an artifact of the within-transformed estimation (verified
#' against an independent dummy-variable regression: exact to floating-point
#' precision when this formula is paired with `beta_lsdv`, for both
#' `effects = "oneway"` and `"twoway"`). Pairing it with `beta_FM` instead
#' (as here, for consistency with the slope actually reported/plotted)
#' makes each *individual* `alpha_i` only an approximation of that unit's
#' true fixed effect; what stays exact regardless of which beta is used is
#' the single number this function actually reports: `mean(alpha_i) +
#' beta' * mean_i(x_i^powers-bar)` always reproduces the panel's true
#' grand-mean `y` exactly, by construction.
#' @rdname turning_points
#' @export
turning_points.pcpr <- function(object, ...) {
  if (identical(object$type, "PMG")) {
    return(pmg_turning_points(object))
  }
  mg <- mg_unit_turning_points(object)
  if (nrow(mg$unit_turning_points) == 0) {
    return(mg$unit_turning_points)
  }
  mg$average
}

#' @keywords internal
mg_unit_turning_points <- function(object) {
  unit_fits <- object$unit_fits
  N <- length(unit_fits)
  xname <- colnames(unit_fits[[1]]$x)[1]
  powers1 <- unit_fits[[1]]$fit$powers[[1]]

  unit_tp <- vector("list", N)
  for (i in seq_len(N)) {
    df_i <- turning_points.cpr(unit_fits[[i]], x_range = "data")
    df_i$unit <- rep(object$units[i], nrow(df_i))
    unit_tp[[i]] <- df_i
  }
  unit_tp <- do.call(rbind, unit_tp)

  if (nrow(unit_tp) == 0) {
    return(list(unit_turning_points = unit_tp, average = unit_tp))
  }

  beta_mg <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const_mg <- get_const_coef(object$coefficients)

  counts <- table(unit_tp$type)
  avg_x <- stats::aggregate(x ~ type, data = unit_tp, FUN = mean)
  avg_x$y <- vapply(avg_x$x, function(xv) const_mg + sum(beta_mg * xv^powers1), numeric(1))
  avg_x$n_units <- as.integer(counts[avg_x$type])
  avg_x <- avg_x[order(avg_x$x), c("x", "y", "type", "n_units")]

  list(unit_turning_points = unit_tp, average = avg_x)
}

#' @keywords internal
pmg_turning_points <- function(object) {
  fit <- object$unit_fits
  powers1 <- seq_len(fit$q)
  beta <- unname(object$coefficients[paste0("x1^", powers1)])
  const <- pmg_average_const(fit, beta, powers1)

  x_range <- range(unlist(lapply(fit$unit_info, function(u) u$x)))
  poly_turning_points(beta, powers1, const = const, x_range = x_range)
}

#' Average, across units, of each unit's own implied fixed effect
#'
#' `alpha_i = mean(y_i) - beta' * mean(x_i^powers)`, using each unit's own
#' *raw* `y`/`x` (`unit_info[[i]]$y`/`$x`, stored by [fit_pooled_panel_cpr()]
#' before any within/demeaning transformation) -- not anything read off the
#' demeaned estimation itself, which would just be ~0 by construction (a
#' relative position, not a real level). See the `@details` on
#' [turning_points.pcpr()] for the exactness this formula does (and does
#' not) guarantee.
#' @keywords internal
pmg_average_const <- function(fit, beta, powers1) {
  alpha_i <- vapply(fit$unit_info, function(u) {
    xp <- gen_power_reg(u$x, powers1)
    mean(u$y) - as.numeric(colMeans(xp) %*% beta)
  }, numeric(1))
  mean(alpha_i)
}
