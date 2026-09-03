# Example: pooled panel CPR with pcpr(type = "pmg").
#
# Unlike pcpr(type = "mg") (N separate slopes, averaged -- full slope
# heterogeneity), type = "pmg" assumes a single common slope shared by all
# N countries, with only the intercepts (oneway) or intercepts + time
# effects (twoway) allowed to vary. Ported from de Jong & Wagner (2016),
# MATLAB source deJongWagner2022/PanelEKC_indiv_eff_only.m (oneway) and
# PanelEKC_two_eff.m (twoway).
#
# Restrictions inherited from the original source: a single integrated
# regressor, and `orders` must be exactly 2 or 3 (the bias-correction
# matrices are only tabulated for those cases).
#
# Run from the package root with: Rscript examples/example_pcpr_pmg.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "methods.R"
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
