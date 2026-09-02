# Example: fitting a cointegrating polynomial regression with cpr().
#
# Run from the package root with: Rscript examples/example_cpr.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "estimators.R", "cpr.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

set.seed(1)

## Simulate a quadratic cointegrating relationship:
##   y_t = 2 + 0.5*x_t + 0.1*x_t^2 + u_t,  x_t integrated of order 1,
##   u_t stationary and autocorrelated (so OLS t-stats are invalid, but
##   FM-OLS ones are asymptotically valid).
Tn <- 250
x <- cumsum(rnorm(Tn))
u <- as.numeric(arima.sim(list(ar = 0.5), n = Tn))
y <- 2 + 0.5 * x + 0.1 * x^2 + u

## Default: constant only, powers 1 and 2 of x, Andrews (1991) bandwidth
## with a Bartlett kernel.
fit <- cpr(y, x, orders = 2)
print(summary(fit))

## Same model with a linear trend added and the VAR(1) pre-whitened
## Andrews & Monahan (1992) bandwidth instead:
fit_trend <- cpr(y, x, orders = 2,
                  deter = make_deterministics(Tn, const = TRUE, trend = TRUE),
                  bandwidth = "AM92")
print(summary(fit_trend))

## A model with an additional stationary regressor:
w <- matrix(rnorm(Tn), ncol = 1, dimnames = list(NULL, "z"))
fit_w <- cpr(y, x, orders = 2, w = w)
print(summary(fit_w))
