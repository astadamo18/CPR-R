# Example: lm()-like formula/data interface for cpr() and pcpr().
#
# Both functions still accept the original vector-based calling convention
# unchanged (y, x, orders, ...) -- this is an additional, not a
# replacement, interface: pass a data frame via `data` and refer to its
# columns either through a formula (`y ~ x1 + x2`) or, for pcpr()'s
# `id`/`time`, a *bare* column name (id = country, no quotes needed --
# lm()-like); quoted strings and raw vectors work too.
#
# Run from the package root with: Rscript examples/example_formula_data.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R",
  "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
panel$noip1000 <- panel$NOIP / 1000
panel$gnipc1000 <- panel$GNIPC / 1000

## ---- cpr(): a single country, formula + data ----
cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]

fit <- cpr(noip1000 ~ gnipc1000, data = cz, orders = 2, kernel = "ba", bandwidth = "And91")
cat("=== cpr(noip1000 ~ gnipc1000, data = cz, orders = 2) ===\n")
print(summary(fit))

## Coefficient names now come from the formula's variable names
## ("gnipc1000^1"/"^2") instead of the generic "x1^1"/"x1^2" the vector
## interface falls back to when `x` has no column names.

## `w` and `deter` can be one-sided formulas against the same `data`:
set.seed(1)
cz$z <- rnorm(nrow(cz))     # an unrelated stationary regressor, for illustration
fit_w <- cpr(noip1000 ~ gnipc1000, data = cz, w = ~z, orders = 2, kernel = "ba", bandwidth = "And91")
cat("\n=== with a stationary regressor via w = ~z ===\n")
print(summary(fit_w))

## ---- pcpr(): the full panel, formula + data, id/time as bare column names ----
fit_mg <- pcpr(noip1000 ~ gnipc1000, data = panel, id = COUNTRY, time = YEAR,
               orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")
cat("\n=== pcpr(noip1000 ~ gnipc1000, data = panel, id = COUNTRY, time = YEAR) ===\n")
print(summary(fit_mg))

## ct_test()/pu_test() work the same way whether the fit came from the
## formula/data interface or the vector interface -- they just dispatch on
## the fitted cpr object either way.
cat("\n=== ct_test(fit) / pu_test(fit) ===\n")
print(ct_test(fit))
print(pu_test(fit))
