# KPSS/Shin-type CT cointegration test for a cointegrating polynomial
# regression. Port of CT_test.m.
#
# Critical value tables are extracted from the original CTPUTests MATLAB
# toolbox's CTcritval/*.mat files (asymptotic critical values, 9 tabulated
# percentiles: 1%, 2.5%, 5%, 10%, 50%, 90%, 95%, 97.5%, 99%). Only the
# d = 0 (intercept), m = 1 (one integrated regressor), p = 2 (max power 2)
# table is bundled so far, matching the deJongWagner2022 NOIP ~ GNIPC +
# GNIPC^2 example. Additional (d, m, p) combinations can be added to
# `.ct_critval_table` the same way: load the corresponding CT_d_*_m_*_p_*.mat
# file and extract its single 1x9 row.

.ct_critval_table <- list(
  "d0_m1_p2" = c(
    `1%` = 0.01879925362135678, `2.5%` = 0.02240143320219408,
    `5%` = 0.02633577319256316, `10%` = 0.032087121877752825,
    `50%` = 0.07422383271386807, `90%` = 0.21328227677587056,
    `95%` = 0.29265553319890514, `97.5%` = 0.37864620187920706,
    `99%` = 0.5044477684037733
  )
)

#' @keywords internal
ct_critval <- function(d, m, p) {
  key <- paste0("d", d, "_m", m, "_p", p)
  tab <- .ct_critval_table[[key]]
  if (is.null(tab)) {
    stop("No CT critical value table bundled for d=", d, ", m=", m, ", p=", p,
         ". Only d=0, m=1, p=2 is available so far; further tables can be added ",
         "to `.ct_critval_table` in R/ct-test.R from the original CTcritval/*.mat files.",
         call. = FALSE)
  }
  tab
}

#' KPSS/Shin-type CT cointegration test for a cointegrating polynomial regression
#'
#' Tests the null hypothesis of cointegration (no unit root in the FM-OLS
#' residuals) against the alternative of no cointegration, using the
#' KPSS/Shin-type statistic of Wagner and co-authors for cointegrating
#' polynomial regressions. Port of `CT_test.m`.
#'
#' An S3 generic: call it either on the raw ingredients (`ct_test(uplus,
#' omega, d, m, p, alpha)`, the [`ct_test.default()`] method) or directly on
#' a fitted [cpr()] object (`ct_test(fit, alpha = ...)`, the
#' [`ct_test.cpr()`] method), which pulls `uplus`, `omega`, `m`, `p`, and
#' `d` out of the fit for you -- `d` is exactly the deterministic choice
#' already made via `cpr(..., deter = ...)`, so there's no need to state
#' it again unless you want to override the inference.
#'
#' @param x Either a numeric vector of residuals (`uplus`, dispatching to
#'   [ct_test.default()]) or a fitted [cpr()] object (dispatching to
#'   [ct_test.cpr()]).
#' @param ... Passed on to the method.
#'
#' @return A list with `statistic` (the CT test statistic), `alpha`,
#'   `critval` (critical value per `alpha`), and `reject` (logical per
#'   `alpha`; `TRUE` means reject the null of cointegration, i.e. evidence
#'   *against* cointegration).
#' @export
ct_test <- function(x, ...) {
  UseMethod("ct_test")
}

#' @describeIn ct_test Default method: supply the residuals and long-run
#'   variance directly.
#' @param omega Estimated long-run variance of the residuals in `x` (e.g.
#'   `fit$fit$Omega_udotv1`).
#' @param d Deterministic specification: `-1` (none), `0` (intercept), `1`
#'   (intercept + trend). Only `0` is currently tabulated.
#' @param m Number of integrated regressors (critical values tabulated up to
#'   4). Only `1` is currently tabulated.
#' @param p Highest power included (critical values tabulated up to 4). Only
#'   `2` is currently tabulated.
#' @param alpha Significance levels to test at; must be a subset of
#'   `c(0.1, 0.05, 0.01)` (the levels available in the tabulated critical
#'   values).
#' @export
ct_test.default <- function(x, omega, d, m, p, alpha = c(0.1, 0.05, 0.01), ...) {
  uplus <- as.numeric(x)
  Tn <- length(uplus)

  partsum <- cumsum(uplus) / sqrt(Tn)
  statistic <- sum(partsum^2) / (Tn * omega)

  crit_table <- ct_critval(d, m, p)
  pct_label <- c(`0.1` = "90%", `0.05` = "95%", `0.01` = "99%")
  lab <- pct_label[as.character(alpha)]
  if (any(is.na(lab))) {
    stop("`alpha` must be a subset of c(0.1, 0.05, 0.01).", call. = FALSE)
  }
  critval <- as.numeric(crit_table[lab])

  list(statistic = statistic, alpha = alpha, critval = critval,
       reject = statistic > critval)
}

#' Infer CT_test.m's `d` (-1/0/1) from a fitted cpr object's `deter`
#'
#' Only succeeds for the standard specifications [make_deterministics()]
#' produces: none, `"const"` only, or `"const"` + `"trend"` (in that
#' order). Anything else (custom `deter`, unusual column names/order) has
#' no well-defined `d` and must be passed explicitly to [ct_test.cpr()].
#' @keywords internal
infer_ct_d <- function(x) {
  kd <- x$kd
  deter_names <- names(x$coefficients)[x$kw + seq_len(kd)]

  if (kd == 0) return(-1L)
  if (kd == 1 && identical(deter_names, "const")) return(0L)
  if (kd == 2 && identical(deter_names, c("const", "trend"))) return(1L)

  stop("Cannot automatically infer `d` from this fit's deterministic ",
       "specification (", kd, " column(s): ",
       if (kd > 0) paste(deter_names, collapse = ", ") else "none",
       "). `d` only has a well-defined meaning here for the standard cases ",
       "make_deterministics() produces: none (d = -1), const only (d = 0), ",
       "or const + trend (d = 1). Pass `d` explicitly.", call. = FALSE)
}

#' @describeIn ct_test `cpr` method: `uplus`, `omega`, `m`, and `p` are read
#'   straight off the fit (`x$fit$residuals`, `x$fit$Omega_udotv1`,
#'   `x$fit$m`, and the highest power in `x$fit$powers`). `d` defaults to
#'   `NULL`, which infers it from the fit's own `deter` (see
#'   [infer_ct_d()]) -- since `d` is exactly the choice you already made
#'   via `cpr(..., deter = ...)`, there is no need to state it twice.
#'   Pass `d` explicitly to override the inference (e.g. for a custom
#'   `deter` the auto-detection can't classify). Works for any estimator
#'   whose fit provides residuals and a long-run variance (currently
#'   `"FMOLS"` and `"DOLS"`).
#' @export
ct_test.cpr <- function(x, d = NULL, alpha = c(0.1, 0.05, 0.01), ...) {
  uplus <- x$fit$residuals
  omega <- x$fit$Omega_udotv1
  if (is.null(uplus) || is.null(omega)) {
    stop("This cpr fit (estimator = '", x$estimator, "') does not provide the ",
         "residuals/long-run variance ct_test() needs.", call. = FALSE)
  }
  if (is.null(d)) d <- infer_ct_d(x)
  m <- x$fit$m
  p <- max(unlist(x$fit$powers))
  ct_test.default(uplus, omega, d = d, m = m, p = p, alpha = alpha)
}
