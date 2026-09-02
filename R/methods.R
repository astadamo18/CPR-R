# S3 methods for "cpr" objects.

#' @export
print.cpr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Cointegrating Polynomial Regression (", x$estimator, ")\n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' @export
summary.cpr <- function(object, ...) {
  structure(object, class = c("summary.cpr", class(object)))
}

#' @export
print.summary.cpr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Cointegrating Polynomial Regression (", x$estimator, ")\n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  bw_label <- if (is.character(x$bandwidth)) x$bandwidth else format(x$bandwidth, digits = digits)
  cat("\nKernel:", x$kernel, " | Bandwidth:", bw_label, "\n")
  cat("Observations used:", x$n_obs, "\n\n")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coef_table, digits = digits, signif.stars = TRUE)
  invisible(x)
}

#' @export
coef.cpr <- function(object, ...) object$coefficients

#' @export
fitted.cpr <- function(object, ...) object$fitted.values

#' @export
residuals.cpr <- function(object, ...) object$residuals

# S3 methods for "pcpr" objects.

#' @export
print.pcpr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Panel Cointegrating Polynomial Regression (", x$estimator, ", type = ", x$type, ")\n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  cat("\nGroup-mean coefficients:\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' @export
summary.pcpr <- function(object, ...) {
  structure(object, class = c("summary.pcpr", class(object)))
}

#' @export
print.summary.pcpr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Panel Cointegrating Polynomial Regression (", x$estimator, ", type = ", x$type, ")\n\n", sep = "")
  cat("Call:\n")
  print(x$call)
  bw_label <- if (is.character(x$bandwidth)) x$bandwidth else format(x$bandwidth, digits = digits)
  cat("\nKernel:", x$kernel, " | Bandwidth:", bw_label, "\n")
  cat("Units (N):", x$n_units, " | Time points per unit (T):", x$n_time, "\n\n")
  cat("Group-mean coefficients:\n")
  stats::printCoefmat(x$coef_table, digits = digits, signif.stars = TRUE)
  invisible(x)
}

#' @export
coef.pcpr <- function(object, ...) object$coefficients
