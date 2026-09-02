# Estimator registry for cpr().
#
# Only FMOLS is implemented for now. DOLS, MOLS, and IMOLS are registered as
# named placeholders so that cpr(estimator = "DOLS") fails with a clear,
# actionable message rather than "unknown estimator" -- and so that adding a
# real implementation later is a matter of replacing the placeholder
# function with a real `fit_*_cpr()` of the same signature
# (y, x, orders, w, deter, kernel, bandwidth) -> list(...).

fit_dols_cpr <- function(...) {
  stop("The DOLS estimator is not implemented yet. cpr()'s estimator interface ",
       "is a registry of functions with signature (y, x, orders, w, deter, kernel, ",
       "bandwidth); DOLS can be added as a new entry in `.cpr_estimators` without ",
       "changing cpr() itself. Only 'FMOLS' is available for now.", call. = FALSE)
}

fit_mols_cpr <- function(...) {
  stop("The MOLS estimator is not implemented yet. See fit_dols_cpr() for how ",
       "to add it via the `.cpr_estimators` registry. Only 'FMOLS' is available for now.",
       call. = FALSE)
}

fit_imols_cpr <- function(...) {
  stop("The IMOLS estimator is not implemented yet. See fit_dols_cpr() for how ",
       "to add it via the `.cpr_estimators` registry. Only 'FMOLS' is available for now.",
       call. = FALSE)
}

.cpr_estimators <- list(
  FMOLS = fit_fmols_cpr,
  DOLS = fit_dols_cpr,
  MOLS = fit_mols_cpr,
  IMOLS = fit_imols_cpr
)
