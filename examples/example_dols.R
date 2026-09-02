# Example: Dynamic OLS (DOLS) estimation with cpr(estimator = "DOLS"),
# compared against FM-OLS on the same CEE data (Czechia: NOIP ~ GNIPC +
# GNIPC^2).
#
# Run from the package root with: Rscript examples/example_dols.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "cpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]
y <- cz$NOIP / 1000
x <- cz$GNIPC / 1000

cat("=== FMOLS ===\n")
print(summary(cpr(y, x, orders = 2, kernel = "ba", bandwidth = "And91")))

## DOLS with n_lag = n_lead = 0: plain OLS on the polynomial regression.
## Unlike FMOLS, DOLS does not drop the first observation in this case
## (28 observations used here vs. 27 for FMOLS).
cat("\n=== DOLS, n_lag = n_lead = 0 (plain OLS, full sample) ===\n")
print(summary(cpr(y, x, orders = 2, estimator = "DOLS", kernel = "ba", bandwidth = "And91")))

## DOLS with one lead and one lag of Delta(GNIPC) as nuisance regressors:
## loses 1 (differencing) + n_lag + n_lead = 3 observations relative to the
## full sample.
cat("\n=== DOLS, n_lag = 1, n_lead = 1 ===\n")
print(summary(cpr(y, x, orders = 2, estimator = "DOLS", kernel = "ba", bandwidth = "And91",
                   n_lag = 1, n_lead = 1)))

## DOLS does not support stationary regressors (`w`) -- this matches the
## original DOLS_CPR.m, not a limitation introduced by the port:
tryCatch(
  cpr(y, x, orders = 2, estimator = "DOLS", w = matrix(rnorm(length(y)), ncol = 1)),
  error = function(e) cat("\nExpected error:", conditionMessage(e), "\n")
)
