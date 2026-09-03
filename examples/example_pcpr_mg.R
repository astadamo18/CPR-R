# Example: mean-group panel CPR with pcpr(type = "mg").
#
# Runs the same NOIP ~ const + GNIPC + GNIPC^2 specification as
# examples/example_panel_cee.R, but as a single pcpr() call instead of a
# manual per-country loop: pcpr(type = "mg") calls cpr() once per country
# (identical estimation to the manual loop) and reports the group-mean
# (Pesaran & Smith, 1995) coefficient across all 13 countries.
#
# Run from the package root with: Rscript examples/example_pcpr_mg.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)

fit <- pcpr(
  panel$NOIP / 1000, panel$GNIPC / 1000,
  id   = panel$COUNTRY,
  time = panel$YEAR,
  orders    = 2,
  estimator = "FMOLS",
  bandwidth = "And91",
  kernel    = "ba",
  type      = "mg"
)

print(summary(fit))

cat("\nUnit-specific (per-country) FM-OLS estimates that were averaged:\n")
print(round(fit$unit_coefficients, 3))

## Each row above is exactly what cpr() would return if called on that
## country by hand -- see fit$unit_fits[["Czechia"]], a full "cpr" object.
