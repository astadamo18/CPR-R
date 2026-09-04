# plot() methods for turning-point analysis: the fitted curve of a
# cointegrating polynomial regression against its (single) integrated
# regressor, with turning point(s) marked and labeled. See the file-level
# comment in R/turning-points.R for why the constant is always included in
# the curve even though it never moves a turning point's x-location.

#' Plot the fitted curve and turning point(s) of a cointegrating polynomial
#' regression
#'
#' @param x A fitted `"cpr"` object with a single integrated regressor.
#' @param y Ignored (required by the [plot()] generic's signature).
#' @param n Number of points in the smooth curve grid.
#' @param x_range Range to draw the curve over and restrict turning points
#'   to; defaults to the observed range of the regressor.
#' @param show_data Add the observed `(x, y)` points for reference.
#' @param digits Rounding used in the turning-point labels.
#' @param xlab,ylab,main Plot labels; `xlab`/`main` default sensibly if left
#'   `NULL`.
#' @param ... Passed on to the underlying [plot()] call.
#' @return Invisibly, the turning-point data frame (see [turning_points()]).
#' @export
plot.cpr <- function(x, y = NULL, n = 200, x_range = NULL, show_data = TRUE,
                      digits = 3, xlab = NULL, ylab = "prediction",
                      main = "Turning point analysis", ...) {
  object <- x
  if (ncol(object$x) != 1) {
    stop("plot.cpr() only supports a fit with a single integrated regressor.", call. = FALSE)
  }
  xname <- colnames(object$x)[1]
  powers1 <- object$fit$powers[[1]]
  beta <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const <- get_const_coef(object$coefficients)

  xr <- if (is.null(x_range)) range(object$x[, 1]) else x_range
  grid <- seq(xr[1], xr[2], length.out = n)
  curve_y <- const + as.numeric(gen_power_reg(grid, powers1) %*% beta)
  tp <- poly_turning_points(beta, powers1, const = const, x_range = xr)

  if (is.null(xlab)) xlab <- xname
  graphics::plot(grid, curve_y, type = "l", lwd = 2, xlab = xlab, ylab = ylab, main = main, ...)
  if (show_data) {
    graphics::points(object$x[, 1], object$y, pch = 16,
                      col = grDevices::adjustcolor("black", 0.35))
  }
  draw_turning_points(tp, digits = digits)
  invisible(tp)
}

#' Plot the fitted curve and (averaged) turning point(s) of a panel
#' cointegrating polynomial regression
#'
#' For `type = "mg"`: the group-mean curve (using [pcpr()]'s own
#' group-mean coefficients, constant included) with each unit's own curve
#' shown faintly for context, and turning points averaged across units by
#' type (see [turning_points.pcpr()]). For `type = "pmg"`: the single pooled
#' curve, using the average implied fixed effect as its constant (see
#' [turning_points.pcpr()]).
#'
#' @param x A fitted `"pcpr"` object.
#' @param y Ignored (required by the [plot()] generic's signature).
#' @param n Number of points in the smooth curve grid(s).
#' @param show_units For `type = "mg"`, also draw each unit's own curve
#'   faintly. Ignored for `type = "pmg"` (there is only one curve).
#' @param digits Rounding used in the turning-point labels.
#' @param xlab,ylab,main Plot labels; default sensibly if left `NULL`.
#' @param ... Passed on to the underlying [plot()] call.
#' @return Invisibly, the turning-point data (see [turning_points.pcpr()]).
#' @export
plot.pcpr <- function(x, y = NULL, n = 200, show_units = TRUE, digits = 3,
                       xlab = NULL, ylab = "prediction", main = NULL, ...) {
  object <- x
  if (identical(object$type, "PMG")) {
    plot_pcpr_pmg(object, n = n, digits = digits, xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    plot_pcpr_mg(object, n = n, show_units = show_units, digits = digits,
                  xlab = xlab, ylab = ylab, main = main, ...)
  }
}

#' @keywords internal
plot_pcpr_mg <- function(object, n, show_units, digits, xlab, ylab, main, ...) {
  unit_fits <- object$unit_fits
  N <- length(unit_fits)
  xname <- colnames(unit_fits[[1]]$x)[1]
  powers1 <- unit_fits[[1]]$fit$powers[[1]]

  xr <- range(unlist(lapply(unit_fits, function(f) f$x[, 1])))
  beta_mg <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const_mg <- get_const_coef(object$coefficients)

  grid <- seq(xr[1], xr[2], length.out = n)
  curve_y <- const_mg + as.numeric(gen_power_reg(grid, powers1) %*% beta_mg)

  if (is.null(xlab)) xlab <- xname
  if (is.null(main)) main <- "Turning point analysis (mean group)"
  graphics::plot(grid, curve_y, type = "l", lwd = 2, xlab = xlab, ylab = ylab, main = main, ...)

  if (show_units) {
    for (i in seq_len(N)) {
      f <- unit_fits[[i]]
      xr_i <- range(f$x[, 1])
      beta_i <- unname(f$coefficients[paste0(xname, "^", powers1)])
      const_i <- get_const_coef(f$coefficients)
      grid_i <- seq(xr_i[1], xr_i[2], length.out = n)
      curve_i <- const_i + as.numeric(gen_power_reg(grid_i, powers1) %*% beta_i)
      graphics::lines(grid_i, curve_i, col = grDevices::adjustcolor("steelblue", 0.35))
    }
  }

  avg <- mg_unit_turning_points(object)$average
  if (nrow(avg) > 0) {
    labels <- paste0("avg. ", avg$type, " (N=", avg$n_units, ")\n(",
                      round(avg$x, digits), ", ", round(avg$y, digits), ")")
    draw_turning_points(avg, digits = digits, labels = labels)
  }
  invisible(avg)
}

#' @keywords internal
plot_pcpr_pmg <- function(object, n, digits, xlab, ylab, main, ...) {
  fit <- object$unit_fits
  powers1 <- seq_len(fit$q)
  beta <- unname(object$coefficients[paste0("x1^", powers1)])
  const <- pmg_average_const(fit, beta, powers1)

  xr <- range(unlist(lapply(fit$unit_info, function(u) u$x)))
  grid <- seq(xr[1], xr[2], length.out = n)
  curve_y <- const + as.numeric(gen_power_reg(grid, powers1) %*% beta)
  tp <- poly_turning_points(beta, powers1, const = const, x_range = xr)

  if (is.null(xlab)) xlab <- "x1"
  if (is.null(main)) main <- "Turning point analysis (pooled panel)"
  graphics::plot(grid, curve_y, type = "l", lwd = 2, xlab = xlab, ylab = ylab, main = main, ...)
  draw_turning_points(tp, digits = digits)
  invisible(tp)
}

#' Mark and label turning points on the current plot
#' @keywords internal
draw_turning_points <- function(tp, digits = 3, labels = NULL) {
  if (nrow(tp) == 0) return(invisible())
  if (is.null(labels)) {
    labels <- paste0(tp$type, "\n(", round(tp$x, digits), ", ", round(tp$y, digits), ")")
  }
  graphics::abline(v = tp$x, lty = 2, col = "red")
  graphics::points(tp$x, tp$y, pch = 19, col = "red", cex = 1.3)
  graphics::text(tp$x, tp$y, labels = labels, pos = 3, offset = 1, col = "red", cex = 0.8, xpd = TRUE)
  invisible()
}
