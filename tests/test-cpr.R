# Base-R sanity tests for the CPR package (no external dependencies).
# Run with: Rscript tests/test-cpr.R

# Assumes this script is run from the package root, e.g. `Rscript tests/test-cpr.R`.
# Sourced in dependency order (estimators.R references fit_fmols_cpr() at
# load time, so fmols.R must be sourced first).
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "imols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R",
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
err <- tryCatch({
  cpr(y, x, orders = 2, estimator = "MOLS")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err))
stopifnot(grepl("not implemented", conditionMessage(err)))
cat("[OK] MOLS raises a clear 'not implemented' error\n")

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

## ---- 6f. IMOLS: coefficient recovery, `w` rejection, and an independent
## re-derivation of the partial-sum-and-augment point estimator plus its
## "up to Lambda^2" sandwich covariance matrix (Vhat_NEW.m), matching
## fit_imols_cpr()'s own computation exactly but written independently
## (a plain loop here vs. apply()/rbind() in production) so a
## transcription bug wouldn't just be replicated by both. ----
fit_imols <- cpr(y, x, orders = 2, estimator = "IMOLS", kernel = "ba", bandwidth = "And91")
stopifnot(inherits(fit_imols, "cpr"))
stopifnot(all(is.finite(fit_imols$coefficients)))
stopifnot(fit_imols$fit$n_obs == Tn)  # no truncation, unlike FM-OLS
stopifnot(abs(fit_imols$coefficients["const"] - 2) < 2)
stopifnot(abs(fit_imols$coefficients["x1^1"] - 0.5) < 0.5)
stopifnot(abs(fit_imols$coefficients["x1^2"] - 0.1) < 0.2)

err_imols_w <- tryCatch({
  cpr(y, x, orders = 2, estimator = "IMOLS", w = matrix(rnorm(Tn), ncol = 1))
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_imols_w))
stopifnot(grepl("does not support stationary regressors", conditionMessage(err_imols_w)))

# Independent re-derivation:
deter_im <- make_deterministics(Tn, const = TRUE, trend = FALSE)
poly_im <- gen_var_poly_terms(matrix(x, ncol = 1), 2, stochastic = FALSE)
X_im <- poly_im$X
Z_im <- cbind(deter_im, X_im)
b_ols_im <- solve(crossprod(Z_im), crossprod(Z_im, y))
u_ols_im <- as.numeric(y - Z_im %*% b_ols_im)
v_im <- diff(x)
v_dm_im <- v_im - mean(v_im)
lv_im <- lr_var(cbind(u_ols_im[2:Tn], v_dm_im), "ba",
                resolve_bandwidth(cbind(u_ols_im[2:Tn], v_dm_im), "ba", "And91"), demean = FALSE)
Lr_im <- lv_im$Omega
Omega_udotv_expected <- as.numeric(Lr_im[1, 1] - Lr_im[1, 2] * Lr_im[2, 1] / Lr_im[2, 2])
stopifnot(isTRUE(all.equal(fit_imols$fit$Omega_udotv, Omega_udotv_expected)))

Sy_im <- cumsum(y)
SD_im <- apply(deter_im, 2, cumsum)
SX_im <- apply(X_im, 2, cumsum)
# Column layout of Xmat_im: [1: cumsum(const), 2: cumsum(x^1), 3: cumsum(x^2),
# 4: x^1 (augmentation), 5: x^2 (augmentation)] -- delta/beta are rows 1-3
# of theta_expected, the augmentation ("gamma") nuisance coefficients rows 4-5.
Xmat_im <- cbind(SD_im, SX_im, X_im)  # "full augmentation"
theta_expected <- solve(crossprod(Xmat_im), crossprod(Xmat_im, Sy_im))
stopifnot(isTRUE(all.equal(unname(fit_imols$coefficients["const"]), unname(theta_expected[1, ]))))
stopifnot(isTRUE(all.equal(unname(fit_imols$coefficients["x1^1"]), unname(theta_expected[2, ]))))
stopifnot(isTRUE(all.equal(unname(fit_imols$coefficients["x1^2"]), unname(theta_expected[3, ]))))

# Independent (loop-based, not apply()/rbind()) re-derivation of Vhat_NEW.m:
vhat_check <- function(g) {
  Tg <- nrow(g); k <- ncol(g)
  S <- matrix(0, Tg, k)
  for (j in seq_len(k)) S[, j] <- cumsum(g[, j])
  DS <- matrix(0, Tg, k)
  total <- S[Tg, ]
  for (t in seq_len(Tg)) {
    prev <- if (t == 1) rep(0, k) else S[t - 1, ]
    DS[t, ] <- total - prev
  }
  gg <- solve(crossprod(g))
  Vh <- gg %*% t(DS)
  Vh %*% t(Vh)
}
varmat_expected <- Omega_udotv_expected * vhat_check(Xmat_im)
se_expected <- sqrt(diag(varmat_expected)[1:3])  # const, x1^1, x1^2 -- rows 4-5 are the augmentation block
stopifnot(isTRUE(all.equal(unname(fit_imols$coef_table["const", "Std. Error"]), unname(se_expected[1]))))
stopifnot(isTRUE(all.equal(unname(fit_imols$coef_table["x1^1", "Std. Error"]), unname(se_expected[2]))))
stopifnot(isTRUE(all.equal(unname(fit_imols$coef_table["x1^2", "Std. Error"]), unname(se_expected[3]))))
cat("[OK] IMOLS matches an independent re-derivation of the partial-sum-and-augment estimator and its sandwich covariance, recovers DGP coefficients, and rejects `w`\n")

## ---- 7. Stationary regressor (w) and explicit trend deterministic ----
w <- matrix(rnorm(Tn), ncol = 1)
deter <- make_deterministics(Tn, const = TRUE, trend = TRUE)
fit_w <- cpr(y, x, orders = 2, w = w, deter = deter)
stopifnot(inherits(fit_w, "cpr"))
stopifnot(all(c("const", "trend") %in% names(fit_w$coefficients)))
stopifnot("w1" %in% names(fit_w$coefficients))
cat("[OK] stationary regressors + trend deterministic work\n")

## ---- 7b. w's HAC covariance reuses the single bandwidth FM_CPR.m
## resolves from [u_ols, Delta(x)] (its `bandw`), rather than re-selecting a
## fresh one from S = w*u_plus. Found by tracing the original source: it
## explicitly reuses `bandw` for every later HAC step, it never re-runs
## And91/NW fresh per series. Using a persistent (AR) w makes the two
## resolved bandwidths differ materially, so this is a real behavior
## difference to lock in, not a coincidental match either way would give. ----
set.seed(11)
w_persistent <- matrix(as.numeric(arima.sim(list(ar = 0.7), n = Tn)), ncol = 1)
fit_w2 <- cpr(y, x, orders = 2, w = w_persistent, kernel = "ba", bandwidth = "And91")

x7b_mat <- as.matrix(x)
v7b <- diff(x7b_mat)
y7b <- y[2:Tn]; x7b <- x7b_mat[2:Tn, , drop = FALSE]; w7b <- w_persistent[2:Tn, , drop = FALSE]
deter7b <- make_deterministics(Tn - 1, const = TRUE, trend = FALSE)
poly7b <- gen_var_poly_terms(x7b, 2, stochastic = TRUE)
J7b <- cbind(deter7b, poly7b$X)
Z7b <- cbind(w7b, J7b)
b_ols7b <- solve(crossprod(Z7b), crossprod(Z7b, y7b))
u_ols7b <- as.numeric(y7b - Z7b %*% b_ols7b)
v_dm7b <- sweep(v7b, 2, colMeans(v7b), "-")
bw_shared_expected <- resolve_bandwidth(cbind(u_ols7b, v_dm7b), "ba", "And91")

u_plus7b <- fit_w2$fit$residuals
S7b <- w7b * matrix(u_plus7b, length(u_plus7b), 1)
SLr_expected <- lr_var(S7b, "ba", bw_shared_expected, demean = FALSE)$Omega
iww7b <- solve(crossprod(w7b))
varmat0_expected <- (Tn - 1) * iww7b %*% SLr_expected %*% iww7b
se_gamma_expected <- sqrt(diag(varmat0_expected))
stopifnot(isTRUE(all.equal(unname(fit_w2$fit$se_gamma), unname(se_gamma_expected))))

# ... and confirm the shared bandwidth genuinely differs from what a fresh
# resolution on S would give -- i.e. the discrepancy this fix closes is
# real, not a case where both approaches happen to agree anyway:
bw_fresh <- resolve_bandwidth(S7b, "ba", "And91")
stopifnot(!isTRUE(all.equal(bw_shared_expected, bw_fresh)))
cat("[OK] w's HAC standard errors reuse the shared [u_ols, Delta(x)] bandwidth, not a freshly-resolved one from S\n")

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
# ct_test()'s bundled critical values were simulated specifically for
# FM-OLS residuals (CT_test.m's own docstring says so) -- there is no
# result establishing a DOLS fit's residuals follow the same null
# distribution, so ct_test() must refuse a non-FMOLS fit rather than
# silently reusing the FM-OLS table for it.
fit_cz_dols <- cpr(cz$NOIP / 1000, cz$GNIPC / 1000, orders = 2, estimator = "DOLS",
                    kernel = "ba", bandwidth = "And91")
err_ct_dols <- tryCatch({ ct_test(fit_cz_dols); NULL }, error = function(e) e)
stopifnot(!is.null(err_ct_dols))
stopifnot(grepl("only supports estimator = 'FMOLS'", conditionMessage(err_ct_dols)))
# pu_test(), by contrast, never touches the fit's estimator-specific
# residuals (only its raw y/x), so it has no such restriction -- it should
# run identically well on a DOLS fit as on an FMOLS one:
pu_cz_dols <- pu_test(fit_cz_dols)
stopifnot(is.finite(pu_cz_dols$statistic))
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
cat("[OK] ct_test() works directly on a cpr object, and rejects non-FMOLS fits (pu_test() does not)\n")

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

# pcpr(type = "mg") is estimator-agnostic by construction (it just calls
# cpr(..., estimator = ...) per unit and averages), so IMOLS works through
# it with no panel-specific code at all -- checked the same way as the
# FMOLS case above (unit-level fits identical to standalone cpr() calls).
fit_mg_imols <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                      orders = 2, kernel = "ba", bandwidth = "And91", type = "mg", estimator = "IMOLS")
stopifnot(inherits(fit_mg_imols, "pcpr"))
for (cname in c("Czechia", "Slovenia")) {
  sub <- panel[panel$COUNTRY == cname, ]
  sub <- sub[order(sub$YEAR), ]
  fit_solo_imols <- cpr(sub$NOIP / 1000, sub$GNIPC / 1000, orders = 2, estimator = "IMOLS",
                         kernel = "ba", bandwidth = "And91")
  stopifnot(identical(fit_solo_imols$coefficients, fit_mg_imols$unit_fits[[cname]]$coefficients))
}
stopifnot(isTRUE(all.equal(unname(fit_mg_imols$coefficients), unname(colMeans(fit_mg_imols$unit_coefficients)))))
# pmg remains FM-OLS-only (its bias-correction algebra is FM-OLS-specific,
# not something any per-unit estimator can be swapped into).
err_pmg_imols <- tryCatch({
  pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
       orders = 2, type = "pmg", estimator = "IMOLS")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_pmg_imols))
stopifnot(grepl("only implements 'FMOLS'", conditionMessage(err_pmg_imols)))
cat("[OK] pcpr(type='mg') works with estimator='IMOLS' for free; pcpr(type='pmg') still rejects it\n")

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

# pmg with an additional *integrated* regressor (ncol(x) > 1): ad hoc
# extension following the same convention as the third-party xtpcmg.ado
# Stata port (R/pooled-panel.R's file-level comment) -- the extra
# regressor enters linearly only (order 1), gets zero bias correction, and
# a separate block-diagonal HC0 standard error.
set.seed(20260918)
panel_z <- panel
panel_z$FAKEZ <- unlist(lapply(split(seq_len(nrow(panel_z)), panel_z$COUNTRY), function(idx) {
  cumsum(rnorm(length(idx)))
})) / 1000

fit_pmg_single <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                        orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")
fit_pmg_multi <- pcpr(panel_z$NOIP / 1000, cbind(GNIPC = panel_z$GNIPC / 1000, FAKEZ = panel_z$FAKEZ),
                       id = panel_z$COUNTRY, time = panel_z$YEAR,
                       orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")
stopifnot(identical(names(fit_pmg_multi$coefficients), c("x1^1", "x1^2", "FAKEZ^1")))
stopifnot(all(is.finite(fit_pmg_multi$coefficients)))
stopifnot(all(is.finite(fit_pmg_multi$coef_table)))
# Adding an (uncorrelated-with-y-by-construction, but still routed through
# the same demeaned regression) extra regressor perturbs the polynomial
# coefficients somewhat, through the shared XX' inverse, but not wildly --
# same sign, same order of magnitude as the single-regressor fit.
stopifnot(sign(fit_pmg_multi$coefficients["x1^1"]) == sign(fit_pmg_single$coefficients["x1^1"]))
stopifnot(sign(fit_pmg_multi$coefficients["x1^2"]) == sign(fit_pmg_single$coefficients["x1^2"]))
# VCV is block-diagonal: no cross-covariance estimated between the
# polynomial block and the additional-regressor block (see file-level
# comment in R/pooled-panel.R -- no derived formula for it exists).
V_multi <- fit_pmg_multi$unit_fits$VCV_FM
stopifnot(identical(dim(V_multi), c(3L, 3L)))
stopifnot(all(V_multi[1:2, 3] == 0))
stopifnot(all(V_multi[3, 1:2] == 0))
# `orders` still only accepts a single integer 2/3 -- it applies to the
# first (polynomial) regressor; a vector/list is rejected with a clear
# error pointing at that restriction.
err_pmg_orders_vec <- tryCatch({
  pcpr(panel_z$NOIP / 1000, cbind(GNIPC = panel_z$GNIPC / 1000, FAKEZ = panel_z$FAKEZ),
       id = panel_z$COUNTRY, time = panel_z$YEAR, orders = c(2, 1), type = "pmg")
  NULL
}, error = function(e) e)
stopifnot(!is.null(err_pmg_orders_vec))
stopifnot(grepl("single integer, 2 or 3", conditionMessage(err_pmg_orders_vec)))
# effects = "twoway" and q = 3 also work with the extra regressor.
fit_pmg_multi_tw <- pcpr(panel_z$NOIP / 1000, cbind(GNIPC = panel_z$GNIPC / 1000, FAKEZ = panel_z$FAKEZ),
                          id = panel_z$COUNTRY, time = panel_z$YEAR,
                          orders = 3, kernel = "ba", bandwidth = "And91", type = "pmg", effects = "twoway")
stopifnot(identical(names(fit_pmg_multi_tw$coefficients), c("x1^1", "x1^2", "x1^3", "FAKEZ^1")))
stopifnot(all(is.finite(fit_pmg_multi_tw$coefficients)))
# turning_points()/plot() still work unchanged (they only ever read off
# the "x1^" powers, ignoring any additional-regressor coefficients).
tp_pmg_multi <- turning_points(fit_pmg_multi)
stopifnot(identical(names(tp_pmg_multi), c("x", "y", "type", "interior")))
cat("[OK] pcpr(type='pmg') supports an additional integrated regressor (linear, zero-corrected, block-diagonal VCV)\n")

# pmg restrictions are rejected with informative errors: q outside {2,3}
# and a stationary regressor `w`. (Additional *integrated* regressors are
# supported -- see the dedicated section below -- just not `w`.)
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
stopifnot(identical(names(tp_linear), c("x", "y", "type", "interior")))

# turning_points() only supports a single integrated regressor.
err_tp_multi <- tryCatch({ turning_points(fit_multi); NULL }, error = function(e) e)
stopifnot(!is.null(err_tp_multi))
stopifnot(grepl("single integrated regressor", conditionMessage(err_tp_multi)))
cat("[OK] turning_points.cpr() matches the closed-form quadratic vertex, and handles the linear/multi-regressor edge cases\n")

# A stationary regressor w's contribution is included in the curve's level,
# evaluated at w's own mean -- not implicitly at w = 0, which can be a
# wild extrapolation whenever 0 falls outside w's actually observed range
# (found via a real GNIPC ~ NOIP + REER fit, where REER never comes near
# zero). The turning point's *location* is still unaffected by w, only
# its level.
set.seed(23)
w_tp <- matrix(rnorm(Tn, mean = 50, sd = 5), ncol = 1)  # never near zero
fit_w_tp <- cpr(y, x, orders = 2, w = w_tp, kernel = "ba", bandwidth = "And91")
tp_w <- turning_points(fit_w_tp)
b1_w <- fit_w_tp$coefficients[["x1^1"]]
b2_w <- fit_w_tp$coefficients[["x1^2"]]
expected_x_w <- -b1_w / (2 * b2_w)
stopifnot(nrow(tp_w) == 1)
stopifnot(isTRUE(all.equal(tp_w$x, expected_x_w)))  # location: unaffected by w
expected_level_w <- fit_w_tp$coefficients[["const"]] + fit_w_tp$coefficients[["w1"]] * mean(w_tp)
expected_y_w <- expected_level_w + b1_w * tp_w$x + b2_w * tp_w$x^2
stopifnot(isTRUE(all.equal(tp_w$y, expected_y_w)))
# ... and this genuinely differs from the old (implicit w = 0) level:
stopifnot(!isTRUE(all.equal(expected_level_w, fit_w_tp$coefficients[["const"]])))
# plot.cpr() must use the exact same level (checked via its invisible return):
plot_dev_w <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_w)
tp_w_plot <- plot(fit_w_tp)
grDevices::dev.off()
unlink(plot_dev_w)
stopifnot(isTRUE(all.equal(tp_w_plot, tp_w)))
cat("[OK] turning_points()/plot() include a stationary regressor w's contribution at w's own mean, not implicitly at w = 0\n")

# Panel mean-group: the reported turning point is that of the group-mean
# curve itself (pcpr()'s own group-mean coefficients, constant included) --
# the same closed-form vertex check as the single-fit case above, not an
# average of each unit's own (individually computed) turning point. Those
# two generally differ (x* = -b1/(2*b2) is a nonlinear function of the
# coefficients, so averaging coefficients first vs. solving for x* first
# don't commute) -- deliberately using the group-mean-curve version here so
# the reported/plotted turning point always sits exactly on the plotted
# curve.
const_mg <- fit_mg$coefficients["const"]
b1_mg <- fit_mg$coefficients["x1^1"]
b2_mg <- fit_mg$coefficients["x1^2"]
expected_x_mg <- unname(-b1_mg / (2 * b2_mg))
tp_mg <- turning_points(fit_mg)
stopifnot(nrow(tp_mg) == 1)
stopifnot(isTRUE(all.equal(tp_mg$x, expected_x_mg)))
expected_y_mg <- unname(const_mg + b1_mg * tp_mg$x + b2_mg * tp_mg$x^2)
stopifnot(isTRUE(all.equal(tp_mg$y, expected_y_mg)))
# ... and it differs from the (no longer reported) average of each unit's
# own turning point, confirming the two really are different quantities:
unit_tp_x <- vapply(fit_mg$unit_fits, function(f) {
  d <- turning_points(f)
  if (nrow(d) == 1) d$x else NA_real_
}, numeric(1))
stopifnot(!isTRUE(all.equal(tp_mg$x, mean(unit_tp_x, na.rm = TRUE))))
cat("[OK] turning_points.pcpr(type='mg') is the group-mean curve's own turning point (matches the plotted curve exactly)\n")

# Panel pooled (pmg): a single common slope, so at most one turning point
# per type; here it happens to fall outside the observed x-range for both
# effects specifications, which is a real (if unexciting) finding -- still
# reported, flagged interior = FALSE, not dropped just for being outside
# the data (checked directly against the unrestricted root).
fit_pmg2 <- pcpr(panel$NOIP / 1000, panel$GNIPC / 1000, id = panel$COUNTRY, time = panel$YEAR,
                  orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")
beta_pmg <- unname(fit_pmg2$coefficients[c("x1^1", "x1^2")])
expected_root_pmg <- -beta_pmg[1] / (2 * beta_pmg[2])
tp_pmg_unrestricted <- poly_turning_points(beta_pmg, c(1, 2), const = 0, x_range = NULL)
stopifnot(isTRUE(all.equal(tp_pmg_unrestricted$x, expected_root_pmg)))
stopifnot(is.na(tp_pmg_unrestricted$interior))  # no x_range given -> not classified
tp_pmg <- turning_points(fit_pmg2)
stopifnot(identical(names(tp_pmg), c("x", "y", "type", "interior")))
stopifnot(nrow(tp_pmg) == 1)  # reported even though it's outside the observed range
stopifnot(isTRUE(!tp_pmg$interior))
cat("[OK] turning_points.pcpr(type='pmg') reports the single common-slope root even when it's outside the observed range (flagged, not dropped)\n")

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

# plot.cpr()'s optional `id` label is purely cosmetic (appended to the plot
# title) -- must not change the returned turning-point data, and must
# default to leaving the title alone.
plot_dev_id <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_id)
tp_cz_id <- plot(fit_cz, id = "Czechia")
grDevices::dev.off()
unlink(plot_dev_id)
stopifnot(isTRUE(all.equal(tp_cz_id, tp_cz)))
cat("[OK] plot.cpr(id = ...) labels the plot without changing the returned data\n")

# plot.cpr()'s axis limits must cover the actual data, not just the fitted
# curve: an observation's residual can put it outside the curve's own
# range (e.g. a noisy point right at the edge of x), and plot()'s first
# call sets the visible region before points() adds the data -- anything
# outside it gets silently clipped unless the axis limits account for it.
set.seed(7)
x_edge <- cumsum(rnorm(40))
y_edge <- 2 + 0.5 * x_edge + 0.05 * x_edge^2 + rnorm(40, sd = 0.5)
y_edge[which.max(x_edge)] <- y_edge[which.max(x_edge)] - 20  # forced large residual right at the edge
fit_edge <- cpr(y_edge, x_edge, orders = 2)

plot_dev_edge <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_edge)
plot(fit_edge)
usr <- graphics::par("usr")  # c(x1, x2, y1, y2) of the actual plotting region
grDevices::dev.off()
unlink(plot_dev_edge)
stopifnot(min(fit_edge$y) >= usr[3] && max(fit_edge$y) <= usr[4])
stopifnot(min(fit_edge$x) >= usr[1] && max(fit_edge$x) <= usr[2])
cat("[OK] plot.cpr()'s axis limits cover the actual data even when a point falls outside the fitted curve's own range\n")

# An extrapolated turning point (outside the observed x-range, like
# fit_pmg2's above) must still actually appear on the plot -- the x-axis
# has to extend far enough to show it, not just report it in the data.
stopifnot(!tp_pmg$interior)  # sanity: this really is the exterior case
plot_dev_extrap <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_extrap)
plot(fit_pmg2)
usr_extrap <- graphics::par("usr")
grDevices::dev.off()
unlink(plot_dev_extrap)
stopifnot(tp_pmg$x >= usr_extrap[1] && tp_pmg$x <= usr_extrap[2])
cat("[OK] plot() extends the axis to actually show a turning point outside the observed data, not just report it\n")

# `extrapolated = FALSE` opts back out of that: the axis should NOT be
# stretched to include an exterior turning point, and turning_points()'s
# own (unaffected) data must still report it either way.
plot_dev_noextrap <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_dev_noextrap)
tp_pmg_noextrap <- plot(fit_pmg2, extrapolated = FALSE)
usr_noextrap <- graphics::par("usr")
grDevices::dev.off()
unlink(plot_dev_noextrap)
stopifnot(!(tp_pmg$x >= usr_noextrap[1] && tp_pmg$x <= usr_noextrap[2]))
stopifnot(isTRUE(all.equal(tp_pmg_noextrap, tp_pmg)))  # returned data unaffected by the drawing choice
cat("[OK] plot(..., extrapolated = FALSE) leaves an exterior turning point off the plot without changing the returned data\n")

cat("\nAll tests passed.\n")
