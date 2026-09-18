# Example: pooled panel CPR with pcpr(type = "pmg").
#
# Unlike pcpr(type = "mg") (N separate slopes, averaged -- full slope
# heterogeneity), type = "pmg" assumes a single common slope shared by all
# N countries, with only the intercepts (oneway) or intercepts + time
# effects (twoway) allowed to vary. Ported from de Jong & Wagner (2016),
# MATLAB source deJongWagner2022/PanelEKC_indiv_eff_only.m (oneway) and
# PanelEKC_two_eff.m (twoway).
#
# Restrictions inherited from the original source: `orders` must be exactly
# 2 or 3 (the bias-correction matrices are only tabulated for those cases,
# for a single polynomial regressor). Additional integrated regressors
# ARE supported (see the second example below), but only as an ad hoc,
# not theoretically derived, extension -- see the file-level comment in
# R/pooled-panel.R.
#
# Run from the package root with: Rscript examples/example_pcpr_pmg.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "imols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)

fit_oneway <- pcpr(
  panel$NOIP / 1000, panel$GNIPC / 1000,
  id = panel$COUNTRY, time = panel$YEAR,
  orders = 2, kernel = "ba", bandwidth = "And91",
  type = "pmg", effects = "oneway"
)
cat("=== pmg, oneway (individual fixed effects only) ===\n")
print(summary(fit_oneway))

fit_twoway <- pcpr(
  panel$NOIP / 1000, panel$GNIPC / 1000,
  id = panel$COUNTRY, time = panel$YEAR,
  orders = 2, kernel = "ba", bandwidth = "And91",
  type = "pmg", effects = "twoway"
)
cat("\n=== pmg, twoway (individual + time fixed effects) ===\n")
print(summary(fit_twoway))

## Three point estimates come out of the pooled model, not just FM-OLS:
cat("\nAll three pooled point estimates (oneway):\n")
cat("beta_lsdv (within/LSDV):    ", round(fit_oneway$unit_fits$beta_lsdv, 4), "\n")
cat("beta_Mod  (bias-corrected):  ", round(fit_oneway$unit_fits$beta_Mod, 4), "\n")
cat("beta_FM   (fully modified): ", round(fit_oneway$unit_fits$beta_FM, 4), "\n")

## For comparison: pcpr(type = "mg") allows the slope to differ by country
## instead of assuming one common slope.
fit_mg <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
               orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")
cat("\nFor comparison, pcpr(type = 'mg') group-mean coefficients:\n")
print(round(fit_mg$coefficients, 4))

## ---- An additional integrated regressor ----
##
## The bundled CEE panel only has NOIP/GNIPC, so a second integrated
## series is fabricated here (a random-walk-like column per country) just
## to demonstrate the mechanism -- not a real macroeconomic variable.
## `orders` still only describes the polynomial regressor (GNIPC, first
## column of x); the extra column (FAKEZ) automatically enters linearly
## only, with zero bias correction and its own (block-diagonal) HC0
## standard error -- see R/pooled-panel.R's file-level comment for why.
set.seed(1)
panel$FAKEZ <- unlist(lapply(split(seq_len(nrow(panel)), panel$COUNTRY), function(idx) {
  cumsum(rnorm(length(idx)))
})) / 1000

fit_pmg_extra <- pcpr(
  panel$NOIP / 1000, cbind(GNIPC = panel$GNIPC / 1000, FAKEZ = panel$FAKEZ),
  id = panel$COUNTRY, time = panel$YEAR,
  orders = 2, kernel = "ba", bandwidth = "And91",
  type = "pmg", effects = "oneway"
)
cat("\n=== pmg with an additional integrated regressor (FAKEZ) ===\n")
print(summary(fit_pmg_extra))
