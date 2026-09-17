# Estimator registry for cpr().
#
# FMOLS (fmols.R), DOLS (dols.R), and IMOLS (imols.R) are implemented. MOLS
# remains a named placeholder so that cpr(estimator = "MOLS") fails with a
# clear, actionable message rather than "unknown estimator" -- and so that
# adding a real implementation later is a matter of replacing the
# placeholder function with a real `fit_*_cpr()` of the same signature
# (y, x, orders, w, deter, kernel, bandwidth, n_lag, n_lead) -> list(...).
# (No standalone single-equation MATLAB source was found for MOLS; only a
# panel-embedded bias-correction inside PanelEKC_*.m, already ported as
# pcpr(type = "pmg")'s beta_Mod.)

fit_mols_cpr <- function(...) {
  stop("The MOLS estimator is not implemented. No standalone single-equation ",
       "MOLS MATLAB source was found to port faithfully (only a panel-embedded ",
       "bias correction inside PanelEKC_*.m, ported as pcpr(type = 'pmg')'s ",
       "beta_Mod). See fit_dols_cpr() in R/dols.R for how to add a new estimator ",
       "via the `.cpr_estimators` registry if one is derived later. ",
       "'FMOLS', 'DOLS', and 'IMOLS' are available for now.", call. = FALSE)
}

.cpr_estimators <- list(
  FMOLS = fit_fmols_cpr,
  DOLS = fit_dols_cpr,
  MOLS = fit_mols_cpr,
  IMOLS = fit_imols_cpr
)
