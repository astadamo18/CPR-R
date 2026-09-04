# Base-R sanity tests for the CPR package (no external dependencies).
# Run with: Rscript tests/test-cpr.R

# Assumes this script is run from the package root, e.g. `Rscript tests/test-cpr.R`.
# Sourced in dependency order (estimators.R references fit_fmols_cpr() at
# load time, so fmols.R must be sourced first).
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R",
  "turning-points.R", "plot.R", "methods.R"
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
for (est in c("MOLS", "IMOLS")) {
  err <- tryCatch({
    cpr(y, x, orders = 2, estimator = est)
    NULL
  }, error = function(e) e)
  stopifnot(!is.null(err))
  stopifnot(grepl("not implemented", conditionMessage(err)))
}
cat("[OK] MOLS/IMOLS raise clear 'not implemented' errors\n")

## ---- 6b. DOLS: plain OLS (n_lag = n_lead = 0) matches OLS on the same design ----
fit_dols0 <- cpr(y, x, orders = 2, estimator = "DOLS", n_lag = 0, n_lead = 0)
stopifnot(inherits(fit_dols0, "cpr"))
stopifnot(all(is.finite(fit_dols0$coefficients)))
# With no lead/lag augmentation, DOLS's beta/delta should equal plain OLS
# on [deter, X] using the FULL (untruncated) sample -- unlike FM-OLS, DOLS
# does not drop the first observation in this case.
X_full <- cbind(x, x^2)
Z_full <- cbind(1, X_full)
b_ols_manual <- solve(crossprod(Z_full), crossprod(Z_full, y))
stopifnot(isTRUE(all.equal(unname(fit_dols0$coefficients["const"]), unname(b_ols_manual[1]), tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(unname(fit_dols0$coefficients["x1^1"]), unname(b_ols_manual[2]), tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(unname(fit_dols0$coefficients["x1^2"]), unname(b_ols_manual[3]), tolerance = 1e-8)))
stopifnot(fit_dols0$fit$n_obs == Tn)  # full sample, no truncation
cat("[OK] DOLS with n_lag=n_lead=0 matches plain OLS on the full (untruncated) sample\n")

## ---- 6c. DOLS with leads/lags: runs, truncates correctly, improves over naive OLS ----
fit_dols2 <- cpr(y, x, orders = 2, estimator = "DOLS", n_lag = 2, n_lead = 2)
stopifnot(all(is.finite(fit_dols2$coefficients)))
stopifnot(fit_dols2$fit$n_obs == Tn - 2 - 2 - 1)  # drop 1 (differencing) + n_lag + n_lead
cat("[OK] DOLS with leads/lags runs and truncates the sample correctly\n")

## ---- 6d. DOLS does not support `w` ----
err_dols_w <- tryCatch({
  cpr(y, x, orders = 2, estimator = "DOLS", w = matrix(rnorm(Tn), ncol = 1))
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_dols_w))
stopifnot(grepl("does not support stationary regressors", conditionMessage(err_dols_w)))
cat("[OK] DOLS rejects `w` with a clear error\n")

## ---- 6e. gen_lead_lag: closed-form check ----
v_ll <- matrix(1:10, ncol = 1)
ll <- gen_lead_lag(v_ll, n_lag = 1, n_lead = 1)
# columns: [contemporaneous, lag1, lead1]
stopifnot(identical(as.numeric(ll[, 1]), as.numeric(v_ll)))
stopifnot(identical(as.numeric(ll[2:10, 2]), as.numeric(v_ll[1:9, 1])))  # lag1
stopifnot(ll[1, 2] == 0)
stopifnot(identical(as.numeric(ll[1:9, 3]), as.numeric(v_ll[2:10, 1])))  # lead1
stopifnot(ll[10, 3] == 0)
cat("[OK] gen_lead_lag matches the expected lag/lead alignment\n")

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

## ---- 10b. ct_test() dispatches on a fitted cpr object, `d` auto-inferred ----
ct_cz_direct <- ct_test(fit_cz)  # no manual fit$fit$residuals/Omega_udotv1/d needed
stopifnot(isTRUE(all.equal(ct_cz_direct$statistic, ct_cz$statistic)))
stopifnot(identical(ct_cz_direct$reject, ct_cz$reject))
# Explicit d still overrides the inference, with an identical result here:
ct_cz_explicit <- ct_test(fit_cz, d = 0)
stopifnot(isTRUE(all.equal(ct_cz_explicit$statistic, ct_cz_direct$statistic)))
# Also works on a DOLS fit (any estimator whose fit exposes residuals + Omega):
fit_cz_dols <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2, estimator = "DOLS",
                    kernel = "ba", bandwidth = "And91")
ct_cz_dols <- ct_test(fit_cz_dols)
stopifnot(is.finite(ct_cz_dols$statistic))
# A fit with a trend infers d = 1, using the real (now-bundled) d=1 table --
# its critical values differ from the d=0 fit's, proving the inference
# actually changed which table was looked up, not just defaulted:
fit_cz_trend <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2,
                     deter = make_deterministics(nrow(cz), trend = TRUE),
                     kernel = "ba", bandwidth = "And91")
ct_trend <- ct_test(fit_cz_trend)
stopifnot(ct_trend$d == 1L)
stopifnot(!isTRUE(all.equal(ct_trend$critval, ct_cz_direct$critval)))
stopifnot(isTRUE(all.equal(unname(ct_trend$critval), unname(ct_critval(1, 1, 2)[c("90%", "95%", "99%")]))))
# A non-standard deter (not const-only or const+trend) can't be classified
# and asks for `d` explicitly rather than guessing:
fit_cz_custom <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2,
                      deter = matrix(rnorm(nrow(cz)), ncol = 1, dimnames = list(NULL, "z")),
                      kernel = "ba", bandwidth = "And91")
err_ambig <- tryCatch({ ct_test(fit_cz_custom); NULL }, error = function(e) e)
stopifnot(!is.null(err_ambig))
stopifnot(grepl("Cannot automatically infer", conditionMessage(err_ambig)))
cat("[OK] ct_test() infers `d` from the fit's deter (or errors clearly when it can't)\n")
cat("[OK] ct_test() works directly on a cpr object (FMOLS and DOLS fits)\n")

## ---- 10c. Full CT critical-value grid (all 48 (d,m,p) combos) loads ----
n_ok <- 0
for (dd in c(-1, 0, 1)) for (mm in 1:4) for (pp in 1:4) {
  tab <- ct_critval(dd, mm, pp)
  stopifnot(length(tab) == 9)
  stopifnot(all(diff(tab) > 0))  # percentiles must be strictly increasing
  n_ok <- n_ok + 1
}
stopifnot(n_ok == 48)
cat("[OK] all 48 bundled CT critical-value tables (d in {-1,0,1}, m,p in {1..4}) load and are monotone\n")

## ---- 10d. print.ct_test() output ----
ct_out <- capture.output(print(ct_cz_direct))
stopifnot(any(grepl("H0: cointegration", ct_out)))
stopifnot(any(grepl("H1: no cointegration", ct_out)))
stopifnot(any(grepl("Test statistic", ct_out)))
stopifnot(any(grepl("Critical values", ct_out)))
stopifnot(any(grepl("Decision", ct_out)))
stopifnot(!any(grepl("p-value|Signif. codes", ct_out)))  # dropped for now: only 9 tabulated percentiles
cat("[OK] print.ct_test() shows statistic, critical values, decisions, and hypotheses (no p-value/stars for now)\n")

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

# type = "pmg": pooled panel estimator (de Jong & Wagner 2016), common
# slope shared by all 13 countries. Run both oneway (individual fixed
# effects) and twoway (individual + time fixed effects), both q = 2 and 3.
for (eff in c("oneway", "twoway")) {
  for (q in c(2, 3)) {
    fit_pmg <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                     orders = q, kernel = "ba", bandwidth = "And91", type = "pmg", effects = eff)
    stopifnot(inherits(fit_pmg, "pcpr"))
    stopifnot(length(fit_pmg$coefficients) == q)
    stopifnot(all(is.finite(fit_pmg$coefficients)))
    stopifnot(all(is.finite(fit_pmg$coef_table)))
    stopifnot(is.null(fit_pmg$unit_coefficients))  # one pooled slope, not one per unit
    stopifnot(all(c("beta_lsdv", "beta_Mod", "beta_FM") %in% names(fit_pmg$unit_fits)))
  }
}
cat("[OK] pcpr(type='pmg') runs for oneway/twoway effects and q=2/3\n")

# pmg restrictions are rejected with informative errors: q outside {2,3},
# a stationary regressor `w`, and more than one integrated regressor.
err_pmg_q <- tryCatch({
  pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
       orders = 4, type = "pmg")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_pmg_q))
stopifnot(grepl("2 or 3", conditionMessage(err_pmg_q)))

err_pmg_w <- tryCatch({
  pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
       w = rnorm(nrow(panel)), orders = 2, type = "pmg")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_pmg_w))
stopifnot(grepl("stationary regressors", conditionMessage(err_pmg_w)))
cat("[OK] pcpr(type='pmg') rejects unsupported orders/`w` with clear errors\n")

# Internal consistency: with a single unit (N=1), the "oneway" pooled model
# reduces to a demeaned-intercept single-series FM-OLS. It should be in the
# same ballpark as cpr() on that unit (not identical: the pooled estimator
# does not truncate the first observation the way cpr() does), as a sanity
# check that the port is not wildly wrong.
cz_only <- panel[panel$COUNTRY == "Czechia", ]
cz_only <- cz_only[order(cz_only$YEAR), ]
fit_pmg_n1 <- pcpr(cz_only$NOIP / 1000, cz_only$GNIPC / 1000, id = cz_only$COUNTRY, time = cz_only$YEAR,
                    orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")
fit_cpr_n1 <- cpr(cz_only$NOIP / 1000, cz_only$GNIPC / 1000, orders = 2, kernel = "ba", bandwidth = "And91")
stopifnot(sign(fit_pmg_n1$coefficients["x1^1"]) == sign(fit_cpr_n1$coefficients["x1^1"]))
stopifnot(sign(fit_pmg_n1$coefficients["x1^2"]) == sign(fit_cpr_n1$coefficients["x1^2"]))
cat("[OK] pmg with N=1 is qualitatively consistent with standalone cpr()\n")

## ---- 12. pu_test(): Phillips-Ouliaris-type PU test ----

cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]
pu_cz <- pu_test(cz$NOIP / 1000, cz$GNIPC / 1000, d = 0, m = 1, orders = 2,
                  kernel = "ba", bandwidth = "And91")
stopifnot(inherits(pu_cz, "pu_test"))
stopifnot(is.finite(pu_cz$statistic))
stopifnot(length(pu_cz$reject) == 3)
cat("[OK] pu_test() runs and returns a finite statistic\n")

# Bundled critical values are the genuine ones extracted from the original
# PUcritval/PU_d_0_m_1_p_2.mat (not invented): spot-check the 5% (95th
# percentile) critical value used for decisions.
stopifnot(isTRUE(all.equal(pu_cz$critval[2], 37.87475517642149)))
cat("[OK] PU critical values match the original PUcritval/*.mat table\n")

# All 48 (d,m,p) combos load (full grid, same as ct_critval()'s).
n_ok_pu <- 0
for (dd in c(-1, 0, 1)) for (mm in 1:4) for (pp in 1:4) {
  tab <- pu_critval(dd, mm, pp)
  stopifnot(length(tab) == 9)
  stopifnot(all(diff(tab) > 0))
  n_ok_pu <- n_ok_pu + 1
}
stopifnot(n_ok_pu == 48)
cat("[OK] all 48 bundled PU critical-value tables load and are monotone\n")

# Out-of-range (d,m,p) still errors informatively.
err_pu <- tryCatch({ pu_critval(0, 5, 2); NULL }, error = function(e) e)
stopifnot(!is.null(err_pu))
stopifnot(grepl("No PU critical value table", conditionMessage(err_pu)))
cat("[OK] pu_test() errors informatively for an untabulated (d, m, p)\n")

# pu_test() dispatches on a fitted cpr object, same as ct_test().
fit_cz_pu <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2, kernel = "ba", bandwidth = "And91")
pu_cz_direct <- pu_test(fit_cz_pu)
stopifnot(isTRUE(all.equal(pu_cz_direct$statistic, pu_cz$statistic)))
stopifnot(pu_cz_direct$d == 0L)
# Works off a DOLS fit too (pu_test doesn't use FM-OLS residuals at all).
fit_cz_pu_dols <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2, estimator = "DOLS",
                       kernel = "ba", bandwidth = "And91")
pu_cz_dols <- pu_test(fit_cz_pu_dols)
stopifnot(is.finite(pu_cz_dols$statistic))
cat("[OK] pu_test() dispatches on a fitted cpr object (FMOLS and DOLS fits)\n")

# print.pu_test() output.
pu_out <- capture.output(print(pu_cz_direct))
stopifnot(any(grepl("H0: no cointegration", pu_out)))
stopifnot(any(grepl("H1: cointegration", pu_out)))
stopifnot(any(grepl("Test statistic", pu_out)))
stopifnot(any(grepl("Critical values", pu_out)))
stopifnot(any(grepl("Decision", pu_out)))
cat("[OK] print.pu_test() shows statistic, critical values, decisions, and hypotheses\n")

## ---- 13. lm()-like formula/data interface for cpr() and pcpr() ----

# cpr(): formula + data gives identical results to the vector interface,
# with nicer coefficient names (from the formula's RHS, not "x1").
fit_cz_formula <- cpr(NOIP1000 ~ GNIPC1000,
                       data = data.frame(NOIP1000 = cz$NOIP / 1000, GNIPC1000 = cz$GNIPC / 1000),
                       orders = 2, kernel = "ba", bandwidth = "And91")
stopifnot(isTRUE(all.equal(unname(fit_cz_formula$coefficients), unname(fit_cz$coefficients))))
stopifnot(identical(names(fit_cz_formula$coefficients), c("const", "GNIPC1000^1", "GNIPC1000^2")))
cat("[OK] cpr(formula, data = ...) matches the vector interface, with formula-derived names\n")

# Formula without `data`, and a formula naming a missing column, both error clearly.
err_no_data <- tryCatch({ cpr(y ~ x, orders = 2); NULL }, error = function(e) e)
stopifnot(!is.null(err_no_data))
stopifnot(grepl("`data` must be supplied", conditionMessage(err_no_data)))

err_missing_col <- tryCatch({
  cpr(NOIP1000 ~ nosuchcolumn, data = data.frame(NOIP1000 = 1:5, x = 1:5), orders = 2)
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_missing_col))
stopifnot(grepl("Column\\(s\\) not found", conditionMessage(err_missing_col)))
cat("[OK] cpr() formula interface errors clearly when `data` is missing or a column isn't found\n")

# `w` and `deter` can also be given as one-sided formulas against `data`.
set.seed(7)
Tn2 <- 150
df_wz <- data.frame(y = 0, x = cumsum(rnorm(Tn2)), z = rnorm(Tn2), trend = seq_len(Tn2))
u2 <- as.numeric(arima.sim(list(ar = 0.4), n = Tn2))
df_wz$y <- 2 + 0.5 * df_wz$x + 0.1 * df_wz$x^2 + 0.8 * df_wz$z + u2
fit_wz <- cpr(y ~ x, data = df_wz, w = ~z, deter = ~trend, orders = 2, kernel = "ba", bandwidth = "And91")
stopifnot(all(c("z", "trend", "x^1", "x^2") %in% names(fit_wz$coefficients)))
stopifnot(abs(fit_wz$coefficients["z"] - 0.8) < 0.3)
cat("[OK] cpr()'s `w`/`deter` accept one-sided formulas against `data`\n")

# pcpr(): formula + data, with id/time as column-name strings, matches the
# vector interface exactly.
panel$noip1000 <- panel$NOIP / 1000
panel$gnipc1000 <- panel$GNIPC / 1000
fit_mg_formula <- pcpr(noip1000 ~ gnipc1000, data = panel, id = "COUNTRY", time = "YEAR",
                        orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")
stopifnot(isTRUE(all.equal(unname(fit_mg_formula$coefficients), unname(fit_mg$coefficients))))
stopifnot(fit_mg_formula$n_units == 13L && fit_mg_formula$n_time == 28L)
cat("[OK] pcpr(formula, data = ..., id = \"...\", time = \"...\") matches the vector interface\n")

# id/time also accept *bare* (unquoted) column names, lm()-like -- not just
# quoted strings.
fit_mg_bare <- pcpr(noip1000 ~ gnipc1000, data = panel, id = COUNTRY, time = YEAR,
                     orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")
stopifnot(isTRUE(all.equal(unname(fit_mg_bare$coefficients), unname(fit_mg$coefficients))))

# Missing `id` altogether, and a bare id naming a column that exists
# nowhere, both error clearly rather than with a confusing NSE internal
# name or a silent wrong answer.
err_no_id <- tryCatch({ pcpr(noip1000 ~ gnipc1000, data = panel, orders = 2); NULL },
                       error = function(e) e)
stopifnot(!is.null(err_no_id))
stopifnot(grepl("`id` must be supplied", conditionMessage(err_no_id)))

err_bad_id <- tryCatch({
  pcpr(noip1000 ~ gnipc1000, data = panel, id = NOSUCHCOLUMN, orders = 2)
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_bad_id))
cat("[OK] pcpr()'s `id`/`time` accept bare (unquoted) column names, lm()-like\n")

# cpr's returned object carries the resolved raw y/x (post data/formula
# lookup, pre estimator truncation) for reuse by other functions.
stopifnot(length(fit_cz_formula$y) == nrow(cz))
stopifnot(isTRUE(all.equal(fit_cz_formula$x[, 1], cz$GNIPC / 1000, check.attributes = FALSE)))
cat("[OK] cpr() stores the resolved raw y/x on the returned object\n")

## ---- 16. turning_points() / plot(): EKC-style turning point analysis ----

# Single fit (Czechia, quadratic): the vertex of a pure quadratic
# const + b1*x + b2*x^2 is at x* = -b1/(2*b2) in closed form.
b1 <- fit_cz$coefficients["x1^1"]
b2 <- fit_cz$coefficients["x1^2"]
expected_x <- unname(-b1 / (2 * b2))
tp_cz <- turning_points(fit_cz)
stopifnot(nrow(tp_cz) == 1)
stopifnot(isTRUE(all.equal(tp_cz$x, expected_x)))
expected_y <- unname(fit_cz$coefficients["const"] + b1 * expected_x + b2 * expected_x^2)
stopifnot(isTRUE(all.equal(tp_cz$y, expected_y)))
stopifnot(tp_cz$type == (if (b2 > 0) "minimum" else "maximum"))

# A purely linear fit (orders = 1) has no turning point.
fit_cz_linear <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 1, kernel = "ba", bandwidth = "And91")
tp_linear <- turning_points(fit_cz_linear)
stopifnot(nrow(tp_linear) == 0)
stopifnot(identical(names(tp_linear), c("x", "y", "type")))

# turning_points() only supports a single integrated regressor.
err_tp_multi <- tryCatch({ turning_points(fit_multi); NULL }, error = function(e) e)
stopifnot(!is.null(err_tp_multi))
stopifnot(grepl("single integrated regressor", conditionMessage(err_tp_multi)))
cat("[OK] turning_points.cpr() matches the closed-form quadratic vertex, and handles the linear/multi-regressor edge cases\n")

# Panel mean-group: the reported average is, by construction, the mean (by
# type) of each unit's own turning point -- checked by recomputing it
# independently here, plus checked that the "proper constant" used for the
# labeled y-value is the panel's own group-mean constant, not zero or an
# unweighted per-unit average of y* values.
unit_tp_list <- lapply(fit_mg$unit_fits, turning_points)
unit_tp_x <- vapply(unit_tp_list, function(d) if (nrow(d) == 1) d$x else NA_real_, numeric(1))
expected_avg_x <- mean(unit_tp_x, na.rm = TRUE)
tp_mg <- turning_points(fit_mg)
stopifnot(nrow(tp_mg) == 1)  # every unit that has one turning point here has a "minimum"
stopifnot(isTRUE(all.equal(tp_mg$x, expected_avg_x)))
stopifnot(tp_mg$n_units == sum(!is.na(unit_tp_x)))
const_mg <- fit_mg$coefficients["const"]
b1_mg <- fit_mg$coefficients["x1^1"]
b2_mg <- fit_mg$coefficients["x1^2"]
expected_y_mg <- unname(const_mg + b1_mg * tp_mg$x + b2_mg * tp_mg$x^2)
stopifnot(isTRUE(all.equal(tp_mg$y, expected_y_mg)))
cat("[OK] turning_points.pcpr(type='mg') averages per-unit turning points and labels them with the group-mean curve\n")

# Panel pooled (pmg): a single common slope, so at most one turning point
# per type; here it happens to fall outside the observed x-range for both
# effects specifications, which is a real (if unexciting) finding, not a
# bug -- checked directly against the unrestricted root.
fit_pmg2 <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                  orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")
beta_pmg <- unname(fit_pmg2$coefficients[c("x1^1", "x1^2")])
expected_root_pmg <- -beta_pmg[1] / (2 * beta_pmg[2])
tp_pmg_unrestricted <- poly_turning_points(beta_pmg, c(1, 2), const = 0, x_range = NULL)
stopifnot(isTRUE(all.equal(tp_pmg_unrestricted$x, expected_root_pmg)))
tp_pmg <- turning_points(fit_pmg2)
stopifnot(identical(names(tp_pmg), c("x", "y", "type")))
stopifnot(nrow(tp_pmg) == 0)  # outside the observed range for this data
cat("[OK] turning_points.pcpr(type='pmg') uses the single common-slope root, restricted to the observed range\n")

# The pmg constant reconstruction (pmg_average_const()) must use each
# unit's own *raw* y/x, not anything derived from the demeaned/within
# estimation (which would just be ~0, a relative position rather than a
# real level). Checked two ways against an independent ground truth:
# (a) per-unit alpha_i, paired with beta_lsdv (not the reported beta_FM),
# must match a genuine dummy-variable (LSDV) lm() regression exactly; (b)
# the single averaged constant this function actually reports must, when
# paired with *any* beta (including the reported beta_FM), exactly
# reproduce the panel's true grand-mean y -- both properties should hold
# for oneway and twoway effects alike.
for (eff in c("oneway", "twoway")) {
  fit_pmg_eff <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                       orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg", effects = eff)
  beta_lsdv_eff <- fit_pmg_eff$unit_fits$beta_lsdv
  beta_fm_eff <- unname(fit_pmg_eff$coefficients[c("x1^1", "x1^2")])

  alpha_i_lsdv <- vapply(fit_pmg_eff$unit_fits$unit_info, function(u) {
    mean(u$y) - as.numeric(colMeans(gen_power_reg(u$x, c(1, 2))) %*% beta_lsdv_eff)
  }, numeric(1))

  if (eff == "oneway") {
    lm_ground_truth <- lm(NOIP1000 ~ factor(COUNTRY) + GNIPC1000 + I(GNIPC1000^2) - 1,
                           data = transform(panel, NOIP1000 = NOIP / 1000, GNIPC1000 = GNIPC / 1000))
    lm_alpha <- coef(lm_ground_truth)[paste0("factor(COUNTRY)", fit_pmg_eff$units)]
  } else {
    lm_ground_truth <- lm(NOIP1000 ~ GNIPC1000 + I(GNIPC1000^2) + factor(COUNTRY) + factor(YEAR),
                           data = transform(panel, NOIP1000 = NOIP / 1000, GNIPC1000 = GNIPC / 1000),
                           contrasts = list(`factor(COUNTRY)` = "contr.sum", `factor(YEAR)` = "contr.sum"))
    cc <- coef(lm_ground_truth)
    alpha_sum <- cc[grepl("factor\\(COUNTRY\\)", names(cc))]
    alpha_all <- c(alpha_sum, -sum(alpha_sum))
    names(alpha_all) <- sort(unique(panel$COUNTRY))
    lm_alpha <- cc[["(Intercept)"]] + alpha_all[fit_pmg_eff$units]
  }
  stopifnot(isTRUE(all.equal(unname(alpha_i_lsdv), unname(lm_alpha), tolerance = 1e-8)))

  grand_ybar <- mean(panel$NOIP / 1000)
  mean_xbar_powers <- colMeans(t(vapply(fit_pmg_eff$unit_fits$unit_info, function(u) {
    colMeans(gen_power_reg(u$x, c(1, 2)))
  }, numeric(2))))
  const_reported <- pmg_average_const(fit_pmg_eff$unit_fits, beta_fm_eff, c(1, 2))
  stopifnot(isTRUE(all.equal(const_reported + as.numeric(mean_xbar_powers %*% beta_fm_eff), grand_ybar)))
}
cat("[OK] pmg's reconstructed constant uses raw (not demeaned) data: matches an independent LSDV regression exactly, and always reproduces the panel's grand-mean y\n")

# plot() methods run without error (redirected to a throwaway pdf() device,
# no display needed) and return the same turning-point data invisibly.
plot_dev_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_file)
invisible_tp_cz <- plot(fit_cz)
invisible_tp_mg <- plot(fit_mg)
invisible_tp_pmg <- plot(fit_pmg2)
grDevices::dev.off()
unlink(plot_dev_file)
stopifnot(isTRUE(all.equal(invisible_tp_cz, tp_cz)))
stopifnot(isTRUE(all.equal(invisible_tp_mg, tp_mg)))
stopifnot(isTRUE(all.equal(invisible_tp_pmg, tp_pmg)))

err_plot_multi <- tryCatch({ plot(fit_multi); NULL }, error = function(e) e)
stopifnot(!is.null(err_plot_multi))
stopifnot(grepl("single integrated regressor", conditionMessage(err_plot_multi)))
cat("[OK] plot.cpr()/plot.pcpr() run without error and return the same turning-point data as turning_points()\n")

cat("\nAll tests passed.\n")
