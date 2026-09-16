# plot() methods for turning-point analysis: the fitted curve of a
# cointegrating polynomial regression against its (single) integrated
# regressor, with turning point(s) marked and labeled. See the file-level
# comment in R/turning-points.R for why the constant -- and any stationary
# regressor w's contribution at w's own mean -- is always included in the
# curve even though neither ever moves a turning point's x-location.
#
# A turning point can fall outside the observed range of x (e.g. pooling a
# common slope across panel units easily pushes the implied vertex beyond
# what any single unit's data covers) -- turning_points() reports it either
# way, flagged via `interior` rather than dropped (see R/turning-points.R),
# and plot() extends the drawn curve just far enough to actually show it,
# dashing the extrapolated segment so it reads as projection rather than
# observed relationship.

#' Extend a data range just enough to include some extra points
#'
#' Used so a turning point outside the observed x-range still gets drawn
#' on the curve instead of being cropped off the plot -- with a little
#' padding so it isn't sitting right on the axis edge. Returns `base_range`
#' unchanged if every point in `extra_points` already falls inside it.
#' @keywords internal
extend_range_for_points <- function(base_range, extra_points) {
  full <- range(c(base_range, extra_points))
  if (isTRUE(all.equal(full, base_range))) return(base_range)
  pad <- 0.05 * diff(full)
  full + c(-pad, pad)
}

#' Draw a curve, dashing whatever part of it falls outside the observed
#' data range
#'
#' @param grid,curve_y The (already computed) curve to draw.
#' @param xr_data The observed x-range; `grid` values inside it are drawn
#'   solid, values outside it (extrapolated) dashed.
#' @keywords internal
draw_curve <- function(grid, curve_y, xr_data) {
  mid <- grid >= xr_data[1] & grid <= xr_data[2]
  graphics::lines(grid[mid], curve_y[mid], lwd = 2, lty = 1)
  left <- grid <= xr_data[1]
  if (sum(left) > 1) graphics::lines(grid[left], curve_y[left], lwd = 2, lty = 2)
  right <- grid >= xr_data[2]
  if (sum(right) > 1) graphics::lines(grid[right], curve_y[right], lwd = 2, lty = 2)
}

#' Plot the fitted curve and turning point(s) of a cointegrating polynomial
#' regression
#'
#' @param x A fitted `"cpr"` object with a single integrated regressor.
#' @param y Ignored (required by the [plot()] generic's signature).
#' @param n Number of points in the smooth curve grid.
#' @param x_range Observed range to draw the curve solid over and classify
#'   turning points against; defaults to the observed range of the
#'   regressor. The curve is still extended (dashed) beyond this range
#'   whenever a turning point falls outside it, so it's never hidden.
#' @param show_data Add the observed `(x, y)` points for reference.
#' @param digits Rounding used in the turning-point labels.
#' @param xlab,ylab,main Plot labels; `xlab`/`main` default sensibly if left
#'   `NULL`.
#' @param id Optional label identifying this fit (e.g. a country/unit name),
#'   appended to `main` as `"<main> - <id>"`. Purely cosmetic -- has no
#'   effect on the fit or the returned turning-point data. `NULL` (default)
#'   leaves `main` as-is. Uses a plain hyphen rather than an em dash: some
#'   graphics devices (e.g. the default bitmap `png()`) fall back to "..."
#'   for characters their font doesn't cover.
#' @param ... Passed on to the underlying [plot()] call.
#' @return Invisibly, the turning-point data frame (see [turning_points()]).
#' @export
plot.cpr <- function(x, y = NULL, n = 200, x_range = NULL, show_data = TRUE,
                      digits = 3, xlab = NULL, ylab = "prediction",
                      main = "Turning point analysis", id = NULL, ...) {
  object <- x
  if (ncol(object$x) != 1) {
    stop("plot.cpr() only supports a fit with a single integrated regressor.", call. = FALSE)
  }
  xname <- colnames(object$x)[1]
  powers1 <- object$fit$powers[[1]]
  beta <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const <- get_level_offset(object$coefficients, object$w)

  xr_data <- if (is.null(x_range)) range(object$x[, 1]) else x_range
  tp <- poly_turning_points(beta, powers1, const = const, x_range = xr_data)
  xr_full <- extend_range_for_points(xr_data, tp$x)

  grid <- sort(unique(c(seq(xr_full[1], xr_full[2], length.out = n), xr_data)))
  curve_y <- const + as.numeric(gen_power_reg(grid, powers1) %*% beta)

  if (is.null(xlab)) xlab <- xname
  if (!is.null(id)) main <- paste0(main, " - ", id)
  dots <- list(...)
  # The axis limits must cover the actual data too, not just the fitted
  # curve -- an observation's residual can easily put it outside the
  # curve's own range (e.g. right at the edge of x, where a single noisy
  # point pulls y well below/above the smooth fit), and plot() sets the
  # visible region from its first two arguments only; points() added
  # afterward are silently clipped if they fall outside it.
  if (show_data) {
    if (is.null(dots$ylim)) dots$ylim <- range(c(curve_y, object$y))
    if (is.null(dots$xlim)) dots$xlim <- range(c(grid, object$x[, 1]))
  }
  do.call(graphics::plot, c(list(grid, curve_y, type = "n",
                                  xlab = xlab, ylab = ylab, main = main), dots))
  draw_curve(grid, curve_y, xr_data)
  if (show_data) {
    graphics::points(object$x[, 1], object$y, pch = 16,
                      col = grDevices::adjustcolor("black", 0.35))
  }
  draw_turning_points(tp, digits = digits)
  invisible(tp)
}

#' Plot the fitted curve and turning point(s) of a panel cointegrating
#' polynomial regression
#'
#' For `type = "mg"`: the group-mean curve (using [pcpr()]'s own
#' group-mean coefficients, constant and any `w` included), with its own
#' turning point marked -- always exactly on the drawn curve (see
#' [turning_points.pcpr()]), extended (dashed) beyond the observed data if
#' the turning point falls outside it. For `type = "pmg"`: the single
#' pooled curve, using the average implied fixed effect as its constant
#' (see [turning_points.pcpr()]).
#'
#' @param x A fitted `"pcpr"` object.
#' @param y Ignored (required by the [plot()] generic's signature).
#' @param n Number of points in the smooth curve grid.
#' @param digits Rounding used in the turning-point labels.
#' @param xlab,ylab,main Plot labels; default sensibly if left `NULL`.
#' @param ... Passed on to the underlying [plot()] call.
#' @return Invisibly, the turning-point data (see [turning_points.pcpr()]).
#' @export
plot.pcpr <- function(x, y = NULL, n = 200, digits = 3,
                       xlab = NULL, ylab = "prediction", main = NULL, ...) {
  object <- x
  if (identical(object$type, "PMG")) {
    plot_pcpr_pmg(object, n = n, digits = digits, xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    plot_pcpr_mg(object, n = n, digits = digits, xlab = xlab, ylab = ylab, main = main, ...)
  }
}

#' @keywords internal
plot_pcpr_mg <- function(object, n, digits, xlab, ylab, main, ...) {
  unit_fits <- object$unit_fits
  xname <- colnames(unit_fits[[1]]$x)[1]
  powers1 <- unit_fits[[1]]$fit$powers[[1]]

  xr_data <- range(unlist(lapply(unit_fits, function(f) f$x[, 1])))
  beta_mg <- unname(object$coefficients[paste0(xname, "^", powers1)])
  const_mg <- get_level_offset(object$coefficients, pooled_w(unit_fits))

  tp <- mg_turning_points(object)
  xr_full <- extend_range_for_points(xr_data, tp$x)

  grid <- sort(unique(c(seq(xr_full[1], xr_full[2], length.out = n), xr_data)))
  curve_y <- const_mg + as.numeric(gen_power_reg(grid, powers1) %*% beta_mg)

  if (is.null(xlab)) xlab <- xname
  if (is.null(main)) main <- "Turning point analysis (mean group)"
  graphics::plot(grid, curve_y, type = "n", xlab = xlab, ylab = ylab, main = main, ...)
  draw_curve(grid, curve_y, xr_data)
  draw_turning_points(tp, digits = digits)
  invisible(tp)
}

#' @keywords internal
plot_pcpr_pmg <- function(object, n, digits, xlab, ylab, main, ...) {
  fit <- object$unit_fits
  powers1 <- seq_len(fit$q)
  beta <- unname(object$coefficients[paste0("x1^", powers1)])
  const <- pmg_average_const(fit, beta, powers1)

  xr_data <- range(unlist(lapply(fit$unit_info, function(u) u$x)))
  tp <- pmg_turning_points(object)
  xr_full <- extend_range_for_points(xr_data, tp$x)

  grid <- sort(unique(c(seq(xr_full[1], xr_full[2], length.out = n), xr_data)))
  curve_y <- const + as.numeric(gen_power_reg(grid, powers1) %*% beta)

  if (is.null(xlab)) xlab <- "x1"
  if (is.null(main)) main <- "Turning point analysis (pooled panel)"
  graphics::plot(grid, curve_y, type = "n", xlab = xlab, ylab = ylab, main = main, ...)
  draw_curve(grid, curve_y, xr_data)
  draw_turning_points(tp, digits = digits)
  invisible(tp)
}

#' Mark and label turning points on the current plot
#'
#' An extrapolated turning point (`interior = FALSE`) is drawn in orange
#' instead of red, with "(extrapolated)" appended to its label -- shown,
#' not hidden, but visually marked as projection rather than something the
#' data actually covers.
#' @keywords internal
draw_turning_points <- function(tp, digits = 3, labels = NULL) {
  if (nrow(tp) == 0) return(invisible())
  extrapolated <- !is.na(tp$interior) & !tp$interior
  if (is.null(labels)) {
    suffix <- ifelse(extrapolated, " (extrapolated)", "")
    labels <- paste0(tp$type, suffix, "\n(", round(tp$x, digits), ", ", round(tp$y, digits), ")")
  }
  col <- ifelse(extrapolated, "darkorange3", "red")
  graphics::abline(v = tp$x, lty = 2, col = col)
  graphics::points(tp$x, tp$y, pch = 19, col = col, cex = 1.3)
  graphics::text(tp$x, tp$y, labels = labels, pos = 3, offset = 1, col = col, cex = 0.8, xpd = TRUE)
  invisible()
}
