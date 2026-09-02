# Base-R sanity tests for the CPR package (no external dependencies).
# Run with: Rscript tests/test-cpr.R

# Assumes this script is run from the package root, e.g. `Rscript tests/test-cpr.R`.
# Sourced in dependency order (estimators.R references fit_fmols_cpr() at
# load time, so fmols.R must be sourced first).
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "estimators.R", "cpr.R", "pcpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

set.seed(42)

## ---- Simulate a simple cointegrated series with a quadratic relationship ----
Tn <- 300
e_x <- rnorm(Tn)
x <- cumsum(e_x)                       # I(1) regressor
u <- arima.sim(list(ar = 0.4), n = Tn) # stationary, autocorrelated errors
y <- 2 + 0.5 * x + 0.1 * x^2 + as.numeric(u)

## ---- 1. Basic fit runs and returns a "cpr" object ----
fit <- cpr(y, x, orders = 2)
stopifnot(inherits(fit, "cpr"))
stopifnot(is.numeric(fit$coefficients))
stopifnot(length(fit$coefficients) == 3)   # const, x^1, x^2
stopifnot(all(is.finite(fit$coefficients)))
stopifnot(all(is.finite(fit$coef_table)))
cat("[OK] basic FMOLS fit runs, coefficients finite\n")

## ---- 2. Coefficient recovery is in the right ballpark ----
stopifnot(abs(fit$coefficients["const"] - 2) < 1.5)
stopifnot(abs(fit$coefficients["x1^1"] - 0.5) < 0.5)
stopifnot(abs(fit$coefficients["x1^2"] - 0.1) < 0.2)
cat("[OK] coefficients close to DGP values\n")

## ---- 3. print/summary methods work ----
out <- capture.output(print(fit))
stopifnot(any(grepl("Cointegrating Polynomial Regression", out)))
out2 <- capture.output(print(summary(fit)))
stopifnot(any(grepl("Coefficients", out2)))
cat("[OK] print/summary methods work\n")

## ---- 4. All kernel x bandwidth combinations that should be valid run without error ----
combos <- list(
  list(bandwidth = "And91", kernel = "ba"),
  list(bandwidth = "And91", kernel = "tr"),
  list(bandwidth = "And91", kernel = "pa"),
  list(bandwidth = "And91", kernel = "qs"),
  list(bandwidth = "NW",    kernel = "ba"),
  list(bandwidth = "NW",    kernel = "pa"),
  list(bandwidth = "NW",    kernel = "qs"),
  list(bandwidth = "AM92",  kernel = "ba"),
  list(bandwidth = 5,       kernel = "ba"),
  list(bandwidth = 5,       kernel = "bo"),
  list(bandwidth = 5,       kernel = "da")
)
for (cmb in combos) {
  f <- cpr(y, x, orders = 2, bandwidth = cmb$bandwidth, kernel = cmb$kernel)
  stopifnot(all(is.finite(f$coefficients)))
}
cat("[OK] all valid bandwidth/kernel combinations run\n")

## ---- 5. Invalid combination (And91 + Bohman) errors informatively ----
err <- tryCatch({
  cpr(y, x, orders = 2, bandwidth = "And91", kernel = "bo")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err))
stopifnot(grepl("not defined", conditionMessage(err)))
cat("[OK] invalid bandwidth/kernel combination errors informatively\n")

## ---- 6. Unimplemented estimators fail with a clear, structured message ----
for (est in c("DOLS", "MOLS", "IMOLS")) {
  err <- tryCatch({
    cpr(y, x, orders = 2, estimator = est)
    NULL
  }, error = function(e) e)
  stopifnot(!is.null(err))
  stopifnot(grepl("not implemented", conditionMessage(err)))
}
cat("[OK] DOLS/MOLS/IMOLS raise clear 'not implemented' errors\n")

## ---- 7. Stationary regressor (w) and explicit trend deterministic ----
w <- matrix(rnorm(Tn), ncol = 1)
deter <- make_deterministics(Tn, const = TRUE, trend = TRUE)
fit_w <- cpr(y, x, orders = 2, w = w, deter = deter)
stopifnot(inherits(fit_w, "cpr"))
stopifnot(all(c("const", "trend") %in% names(fit_w$coefficients)))
stopifnot("w1" %in% names(fit_w$coefficients))
cat("[OK] stationary regressors + trend deterministic work\n")

## ---- 8. Multiple integrated regressors with per-column orders (list form) ----
x2 <- cumsum(rnorm(Tn))
X2 <- cbind(x1 = x, x2 = x2)
fit_multi <- cpr(y, X2, orders = list(c(1, 2), c(1)))
stopifnot(length(fit_multi$fit$powers) == 2)
stopifnot(identical(fit_multi$fit$powers[[1]], c(1, 2)))
stopifnot(identical(fit_multi$fit$powers[[2]], 1))
cat("[OK] multiple integrated regressors with explicit per-column orders work\n")

## ---- 9. Low-level building blocks: sanity checks against known facts ----

# lr_weights: Bartlett with integer bandwidth M gives linearly decaying
# weights 1 - j/M for j = 1..M-1, then zero.
lw <- lr_weights(50, "ba", 5)
stopifnot(isTRUE(all.equal(lw$w[1:4], 1 - (1:4) / 5)))
stopifnot(all(lw$w[5:49] == 0))

# lr_var: for iid data (no autocorrelation), Omega should be close to Sigma.
set.seed(1)
iid_u <- matrix(rnorm(5000), ncol = 1)
lv <- lr_var(iid_u, "ba", 5, demean = FALSE)
stopifnot(abs(lv$Omega[1, 1] - lv$Sigma[1, 1]) < 0.2)

# gen_cpr_corr_vec: corr term for power 1 is T, for power 2 is 2*sum(x).
xx <- rnorm(20)
cv <- gen_cpr_corr_vec(xx, c(1, 2))
stopifnot(isTRUE(all.equal(cv[1], length(xx))))
stopifnot(isTRUE(all.equal(cv[2], 2 * sum(xx))))

cat("[OK] low-level building blocks match known closed forms\n")

## ---- 10. Regression test against the original MATLAB FM_OLS_panel.m /
## CT_test.m output (deJongWagner2022 CEE panel: NOIP ~ GNIPC + GNIPC^2) ----
panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]
fit_cz <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2, kernel = "ba", bandwidth = "And91")
ct_cz <- ct_test(fit_cz$fit$residuals, fit_cz$fit$Omega_udotv1, d = 0, m = 1, p = 2)

# Reference values from the original MATLAB output (screenshot / FM_OLS_panel.m):
# const=13.327 (p=0.000), GNIPC=-1.219 (p=0.000), GNIPC^2=0.014 (p=0.000), CT=0.101, no rejection at 10/5/1%.
stopifnot(isTRUE(all.equal(round(fit_cz$coef_table["const", "Estimate"], 3), 13.327)))
stopifnot(isTRUE(all.equal(round(fit_cz$coef_table["x1^1", "Estimate"], 3), -1.219)))
stopifnot(isTRUE(all.equal(round(fit_cz$coef_table["x1^2", "Estimate"], 3), 0.014)))
stopifnot(isTRUE(all.equal(round(ct_cz$statistic, 3), 0.101)))
stopifnot(!any(ct_cz$reject))
cat("[OK] matches original MATLAB FM_OLS_panel.m / CT_test.m output (Czechia)\n")

## ---- 11. pcpr(): mean-group panel estimator ----

fit_mg <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
               orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")
stopifnot(inherits(fit_mg, "pcpr"))
stopifnot(fit_mg$n_units == 13L)
stopifnot(fit_mg$n_time == 28L)

# Design contract: the per-unit estimation inside pcpr(mg) must be *identical*
# to calling cpr() on that unit directly (same function, not a parallel
# reimplementation) -- checked here bit-for-bit for two units.
for (cname in c("Czechia", "Slovenia")) {
  sub <- panel[panel$COUNTRY == cname, ]
  sub <- sub[order(sub$YEAR), ]
  fit_solo <- cpr(sub$NOIP / 1000, sub$GNIPC / 1000, orders = 2, kernel = "ba", bandwidth = "And91")
  stopifnot(identical(fit_solo$coefficients, fit_mg$unit_fits[[cname]]$coefficients))
}
cat("[OK] pcpr(type='mg') unit-level fits are identical to standalone cpr() calls\n")

# Group-mean coefficient is exactly the column mean of the unit coefficients.
stopifnot(isTRUE(all.equal(unname(fit_mg$coefficients), unname(colMeans(fit_mg$unit_coefficients)))))
cat("[OK] group-mean coefficient equals the mean of the unit-specific estimates\n")

# print/summary work.
out3 <- capture.output(print(summary(fit_mg)))
stopifnot(any(grepl("Group-mean coefficients", out3)))
cat("[OK] pcpr print/summary methods work\n")

# Unbalanced panel is rejected with an informative error.
panel_unbalanced <- panel[!(panel$COUNTRY == "Czechia" & panel$YEAR == max(panel$YEAR)), ]
err_bal <- tryCatch({
  pcpr(panel_unbalanced$NOIP / 1000, panel_unbalanced$GNIPC / 1000,
       id = panel_unbalanced$COUNTRY, time = panel_unbalanced$YEAR, orders = 2)
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_bal))
stopifnot(grepl("balanced panel", conditionMessage(err_bal)))
cat("[OK] pcpr() rejects an unbalanced panel with a clear error\n")

# type = "pmg" is a registered-but-unimplemented placeholder.
err_pmg <- tryCatch({
  pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
       orders = 2, type = "pmg")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_pmg))
stopifnot(grepl("not implemented", conditionMessage(err_pmg)))
cat("[OK] pcpr(type='pmg') raises a clear 'not implemented' error\n")

cat("\nAll tests passed.\n")
