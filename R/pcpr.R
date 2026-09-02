# Panel cointegrating polynomial regression.
#
# type = "mg" (mean group): the *same* single-unit FM-OLS estimation used by
# cpr() is run once per cross-sectional unit -- literally calling cpr() N
# times, not a separate re-implementation -- and the unit-specific
# coefficients are averaged (Pesaran & Smith, 1995, mean-group estimator).
# This mirrors GroupMeanFMOLS.m conceptually (mean of N individual FM-OLS
# estimators) but, per design, reuses cpr() itself as the unit-level engine
# rather than GroupMeanFMOLS.m's separate Frisch-Waugh-demeaned code path,
# so unit-level results are guaranteed identical to calling cpr() on that
# unit by hand.
#
# type = "pmg" (panel / pooled mean group, common slope across units) is
# registered as a placeholder for now; see .pcpr_types below.
#
# Consistency rule: pcpr() requires a *balanced* panel (every unit observed
# at the same set of time points). This guarantees that the long-run
# variance for every unit's model is estimated over the same number of time
# points T -- the "how many timepoints was the variance estimated over"
# question is answered identically for every unit. The auto-selected
# bandwidth itself is still allowed to differ across units (that is the
# point of allowing full slope heterogeneity in "mg"); only the sample
# length T feeding into it is held fixed across units. A fixed numeric
# `bandwidth` applies identically to all units either way.

fit_mg_pcpr <- function(y_list, x_list, orders, w_list, deter_list,
                         estimator, bandwidth, kernel, unit_names) {
  N <- length(y_list)
  fits <- vector("list", N)
  names(fits) <- unit_names

  for (i in seq_len(N)) {
    fits[[i]] <- cpr(y_list[[i]], x_list[[i]], orders = orders,
                      w = w_list[[i]], deter = deter_list[[i]],
                      estimator = estimator, bandwidth = bandwidth, kernel = kernel)
  }

  coef_names <- names(fits[[1]]$coefficients)
  same_names <- vapply(fits, function(f) identical(names(f$coefficients), coef_names), logical(1))
  if (!all(same_names)) {
    stop("Inconsistent coefficients across units: ", paste(unit_names[!same_names], collapse = ", "),
         ". Check that `orders`, `w`, and `deter` give the same specification for every unit.",
         call. = FALSE)
  }

  B <- do.call(rbind, lapply(fits, function(f) f$coefficients))
  rownames(B) <- unit_names

  betaMG <- colMeans(B)
  if (N > 1) {
    # Pesaran & Smith (1995) mean-group variance: the between-unit sample
    # variance of the individual estimates, scaled by 1/N. Nonparametric in
    # the sense that it does not require a model for cross-sectional
    # dependence; it is the standard MG inference approach.
    covB <- stats::cov(B)
    varMG <- diag(covB) / N
  } else {
    varMG <- rep(NA_real_, ncol(B))
  }

  list(fits = fits, unit_coefficients = B, coefficients = betaMG, se = sqrt(varMG))
}

fit_pmg_pcpr <- function(...) {
  stop("The 'pmg' (panel / pooled mean group) type is not implemented yet. ",
       "pcpr()'s type interface is a registry of functions with signature ",
       "(y_list, x_list, orders, w_list, deter_list, estimator, bandwidth, kernel, ",
       "unit_names); 'pmg' can be added as a new entry in `.pcpr_types` without ",
       "changing pcpr() itself. Only 'mg' is available for now.", call. = FALSE)
}

.pcpr_types <- list(MG = fit_mg_pcpr, PMG = fit_pmg_pcpr)

#' Fit a panel cointegrating polynomial regression (panel CPR)
#'
#' @param y Numeric vector, the (stacked, long-format) dependent variable.
#' @param x Numeric matrix (or vector) of I(1) regressors, stacked long-format.
#' @param id Vector identifying the cross-sectional unit of each row of
#'   `y`/`x`. The panel must be balanced: every unit must have the same
#'   number of time observations.
#' @param time Optional time index per row, used to sort observations within
#'   each unit before differencing/estimation. If `NULL`, rows are assumed
#'   to already be sorted by time within each unit.
#' @param orders,w,deter,estimator,bandwidth,kernel Passed through to
#'   [cpr()] for each unit; see its documentation. `w` and `deter`, if
#'   supplied, must be stacked long-format like `y`/`x`.
#' @param type Panel estimator type. `"mg"` (mean group; default): average
#'   of N unit-specific [cpr()] fits, with Pesaran & Smith (1995) between-unit
#'   inference. `"pmg"` (panel / pooled mean group, common slope across
#'   units) is reserved for future use and currently raises an informative
#'   error.
#'
#' @return An object of class `"pcpr"`, with `print()` and `summary()`
#'   methods. `object$unit_fits` holds the individual per-unit `"cpr"`
#'   objects (identical to calling [cpr()] on that unit directly), and
#'   `object$unit_coefficients` the matrix of unit-specific coefficients
#'   that were averaged.
#' @export
pcpr <- function(y, x, id, time = NULL, orders, w = NULL, deter = NULL,
                  estimator = "FMOLS", bandwidth = "And91", kernel = "ba",
                  type = "mg") {

  cl <- match.call()

  y <- as.numeric(y)
  x <- as.matrix(x)
  Tn_total <- length(y)

  if (nrow(x) != Tn_total) stop("`y` and `x` must have the same number of observations.", call. = FALSE)
  if (length(id) != Tn_total) stop("`id` must have the same length as `y`.", call. = FALSE)
  if (!is.null(time) && length(time) != Tn_total) stop("`time` must have the same length as `y`.", call. = FALSE)

  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))

  if (!is.null(w)) {
    w <- as.matrix(w)
    if (nrow(w) != Tn_total) stop("`w` must have the same number of observations as `y`.", call. = FALSE)
  }
  if (!is.null(deter)) {
    deter <- as.matrix(deter)
    if (nrow(deter) != Tn_total) stop("`deter` must have the same number of observations as `y`.", call. = FALSE)
  }

  id <- as.character(id)
  units <- unique(id)
  N <- length(units)

  y_list <- vector("list", N)
  x_list <- vector("list", N)
  w_list <- vector("list", N)
  deter_list <- vector("list", N)
  Tn_unit <- integer(N)

  for (i in seq_len(N)) {
    idx <- which(id == units[i])
    if (!is.null(time)) idx <- idx[order(time[idx])]
    y_list[[i]] <- y[idx]
    x_list[[i]] <- x[idx, , drop = FALSE]
    # Note: `w_list[[i]] <- NULL` would *delete* element i (an R list
    # gotcha), so only assign when there is something to assign; w_list /
    # deter_list start out as N-element lists of NULL from vector("list", N).
    if (!is.null(w)) w_list[[i]] <- w[idx, , drop = FALSE]
    if (!is.null(deter)) deter_list[[i]] <- deter[idx, , drop = FALSE]
    Tn_unit[i] <- length(idx)
  }

  if (length(unique(Tn_unit)) > 1) {
    stop("pcpr() requires a balanced panel: every unit must have the same number of ",
         "time observations, but found ", paste(unique(Tn_unit), collapse = ", "),
         ". This keeps the long-run variance estimation for every unit's model based ",
         "on the same number of time points.", call. = FALSE)
  }
  Tn <- Tn_unit[1]

  if (is.null(deter)) {
    const_deter <- make_deterministics(Tn, const = TRUE, trend = FALSE)
    deter_list <- replicate(N, const_deter, simplify = FALSE)
  }

  type <- match.arg(toupper(type), names(.pcpr_types))
  fit_fun <- .pcpr_types[[type]]
  fit <- fit_fun(y_list, x_list, orders = orders, w_list = w_list, deter_list = deter_list,
                 estimator = estimator, bandwidth = bandwidth, kernel = kernel, unit_names = units)

  tval <- fit$coefficients / fit$se
  pval <- 2 * stats::pnorm(-abs(tval))
  coef_table <- cbind(Estimate = fit$coefficients, `Std. Error` = fit$se,
                       `z value` = tval, `Pr(>|z|)` = pval)

  structure(
    list(
      call = cl,
      estimator = estimator, kernel = kernel, bandwidth = bandwidth, type = type,
      units = units, n_units = N, n_time = Tn,
      coefficients = fit$coefficients, coef_table = coef_table,
      unit_coefficients = fit$unit_coefficients, unit_fits = fit$fits
    ),
    class = "pcpr"
  )
}
