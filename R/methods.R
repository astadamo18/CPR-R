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
