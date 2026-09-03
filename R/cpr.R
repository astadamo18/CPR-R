# User-facing entry point for single-equation cointegrating polynomial
# regression.

#' Build a matrix of deterministic regressors
#'
#' Convenience helper for the `deter` argument of [cpr()].
#'
#' @param Tn Number of observations.
#' @param const Include a constant column.
#' @param trend Include a linear trend column (`1, 2, ..., Tn`).
#' @return A `Tn x q` matrix with columns named `"const"` / `"trend"` as
#'   applicable.
#' @export
make_deterministics <- function(Tn, const = TRUE, trend = FALSE) {
  cols <- list()
  if (const) cols$const <- rep(1, Tn)
  if (trend) cols$trend <- seq_len(Tn)
  if (length(cols) == 0) {
    stop("At least one of `const` or `trend` must be TRUE.", call. = FALSE)
  }
  do.call(cbind, cols)
}

#' Fit a cointegrating polynomial regression (CPR)
#'
#' Single-equation estimation of
#' \deqn{y_t = \gamma' w_t + \delta' deter_t + \beta' X_t + u_t}
#' where `X_t` collects the requested powers of the integrated (I(1))
#' regressors `x_t`, `w_t` are stationary (I(0)) regressors, and `deter_t`
#' are deterministic components (constant, trend, ...).
#'
#' @param y Either a numeric response vector (or single-column matrix), the
#'   I(1) dependent variable -- or, lm()-like, a two-sided formula
#'   `y ~ x1 + x2` naming columns of `data` (in which case `data` is
#'   required and `x` is ignored). A formula's right-hand side names
#'   columns verbatim; see the file-level comment in `R/formula-data.R` for
#'   why transformed terms like `log(x1)` are not supported directly.
#' @param x Numeric matrix (or vector) of I(1) regressors, one column per
#'   variable, when `y` is not a formula. At least one column is required.
#' @param data A data frame to look up columns in, for the `y ~ x1 + x2`
#'   formula form (required for it) and/or for `w`/`deter` given as
#'   one-sided formulas (e.g. `w = ~ w1 + w2`).
#' @param orders Powers of `x` to include. Either a single integer (same max
#'   order `1:order` for every column of `x`), a numeric vector of length
#'   `ncol(x)` (per-column max order), or a list of length `ncol(x)` giving
#'   explicit (possibly non-consecutive) powers per column, e.g.
#'   `list(c(1, 2), c(1, 3))`.
#' @param w Optional matrix (or vector) of stationary (I(0)) regressors.
#' @param deter Optional matrix of deterministic regressors. Defaults to a
#'   constant only (via [make_deterministics()]); pass e.g.
#'   `make_deterministics(length(y), trend = TRUE)` to add a linear trend.
#' @param estimator Estimation method: `"FMOLS"` (fully modified OLS) or
#'   `"DOLS"` (dynamic OLS: the polynomial regression augmented with leads
#'   and lags of `Delta(x)`, see `n_lag`/`n_lead`; does not support `w`).
#'   `"MOLS"`, `"IMOLS"` are reserved for future extensions and currently
#'   raise an informative error.
#' @param bandwidth Bandwidth selection for the long-run variance estimator:
#'   `"And91"` (Andrews, 1991; default), `"AM92"` (Andrews & Monahan, 1992,
#'   VAR(1) pre-whitened), `"NW"` (Newey & West, 1994), or a fixed numeric
#'   bandwidth.
#' @param kernel Kernel function for the long-run variance estimator:
#'   `"ba"` (Bartlett; default), `"tr"` (truncated), `"pa"` (Parzen), `"bo"`
#'   (Bohman), `"da"` (Daniell), or `"qs"` (Quadratic Spectral). Note that
#'   `bandwidth = "And91"` only supports `"tr"`, `"ba"`, `"pa"`, `"qs"`
#'   (plus `"th"`, not offered here as a kernel for the LR-variance weights
#'   themselves), and `bandwidth = "NW"` only supports `"ba"`, `"pa"`,
#'   `"qs"`.
#' @param n_lag,n_lead Only used when `estimator = "DOLS"`: number of
#'   lagging / leading first differences of `x` to include as additional
#'   (nuisance) regressors. `0, 0` (the default) means plain OLS on the
#'   polynomial regression, with no lead/lag augmentation and, unlike
#'   `"FMOLS"`, no dropped first observation.
#'
#' @details
#' `Delta(x_t) = x_t - x_{t-1}` is always demeaned before it enters the
#' long-run variance estimation step (Wagner & Reichold, 2023, Remark 5);
#' this demeaning is hardwired and not a user-controlled option.
#'
#' @return An object of class `"cpr"`, with `print()` and `summary()`
#'   methods. Also carries the resolved `y`/`x` (post `data`/formula
#'   lookup, pre any estimator-specific truncation) as `$y`/`$x`, so other
#'   functions needing the original series (e.g. [pu_test()]) can be
#'   dispatched straight off the fit.
#' @export
cpr <- function(y, x = NULL, orders, w = NULL, deter = NULL,
                 estimator = "FMOLS",
                 bandwidth = "And91",
                 kernel = "ba",
                 n_lag = 0, n_lead = 0,
                 data = NULL) {

  cl <- match.call()

  if (inherits(y, "formula")) {
    fd <- extract_formula_xy(y, data)
    y <- fd$y
    x <- fd$x
  }
  if (!is.null(w) && inherits(w, "formula")) w <- resolve_formula_vars(w, data, "w")
  if (!is.null(deter) && inherits(deter, "formula")) deter <- resolve_formula_vars(deter, data, "deter")

  if (is.null(x)) {
    stop("`x` must be supplied, either directly or implicitly via the right-hand ",
         "side of a `y ~ x1 + x2` formula (with `data`).", call. = FALSE)
  }

  y <- as.matrix(y)
  x <- as.matrix(x)
  Tn <- nrow(x)

  if (nrow(y) != Tn) stop("`y` and `x` must have the same number of observations.", call. = FALSE)
  if (ncol(x) < 1) stop("`x` must have at least one column.", call. = FALSE)

  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  xnames <- colnames(x)

  if (is.null(deter)) {
    deter <- make_deterministics(Tn, const = TRUE, trend = FALSE)
  } else {
    deter <- as.matrix(deter)
    if (is.null(colnames(deter))) colnames(deter) <- paste0("deter", seq_len(ncol(deter)))
  }

  if (!is.null(w)) {
    w <- as.matrix(w)
    if (ncol(w) == 0) {
      w <- NULL
    } else if (is.null(colnames(w))) {
      colnames(w) <- paste0("w", seq_len(ncol(w)))
    }
  }

  estimator <- match.arg(toupper(estimator), names(.cpr_estimators))
  kernel <- match.arg(kernel, c("tr", "ba", "pa", "bo", "da", "qs"))
  if (is.character(bandwidth)) bandwidth <- match.arg(bandwidth, c("And91", "AM92", "NW"))

  fit_fun <- .cpr_estimators[[estimator]]
  fit <- fit_fun(y = y, x = x, orders = orders, w = w, deter = deter,
                 kernel = kernel, bandwidth = bandwidth, n_lag = n_lag, n_lead = n_lead)

  beta_names <- unlist(mapply(function(nm, pw) paste0(nm, "^", pw),
                               xnames, fit$powers, SIMPLIFY = FALSE))

  coef_names <- c(colnames(w), colnames(deter), beta_names)
  coefficients <- c(fit$coef_gamma, fit$coef_delta, fit$coef_beta)
  se <- c(fit$se_gamma, fit$se_delta, fit$se_beta)
  names(coefficients) <- coef_names
  names(se) <- coef_names

  tval <- coefficients / se
  pval <- 2 * stats::pnorm(-abs(tval))

  coef_table <- cbind(Estimate = coefficients, `Std. Error` = se,
                       `z value` = tval, `Pr(>|z|)` = pval)

  structure(
    list(
      call = cl,
      estimator = estimator,
      kernel = kernel,
      bandwidth = bandwidth,
      coefficients = coefficients,
      coef_table = coef_table,
      fitted.values = fit$fitted,
      residuals = fit$residuals,
      residuals_ols = fit$residuals_ols,
      n_obs = fit$n_obs,
      kw = fit$kw, kd = fit$kd, m = fit$m,
      y = as.numeric(y), x = x,
      fit = fit
    ),
    class = "cpr"
  )
}
