# Example: ct_test()'s print output -- test statistic, critical values,
# H0/H1, an approximate (interpolated) p-value, and significance stars.
#
# Run from the package root with: Rscript examples/example_ct_test_print.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]
y <- cz$NOIP / 1000
x <- cz$GNIPC / 1000

fit <- cpr(y, x, orders = 2, kernel = "ba", bandwidth = "And91")

cat("=== d = 0 (intercept only, inferred from fit's deter) ===\n")
print(ct_test(fit))

## A fit with a trend uses the d = 1 critical value table automatically --
## a genuinely different table, not the same one reused.
fit_trend <- cpr(y, x, orders = 2, deter = make_deterministics(nrow(cz), trend = TRUE),
                  kernel = "ba", bandwidth = "And91")
cat("\n\n=== d = 1 (intercept + trend, also inferred) ===\n")
print(ct_test(fit_trend))

## The whole (d, m, p) critical value table is picked automatically -- here
## from a cubic (p = 3) specification, no table management needed by hand.
fit_cubic <- cpr(y, x, orders = 3, kernel = "ba", bandwidth = "And91")
cat("\n\n=== p = 3 (cubic polynomial regressor, also inferred) ===\n")
print(ct_test(fit_cubic))
