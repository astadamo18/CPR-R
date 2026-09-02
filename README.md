# CPR

Cointegrating Polynomial Regressions in R, following Wagner and co-authors
(e.g. Wagner; Wagner & Reichold; de Jong & Wagner; Veldhuis & Wagner). This is
an R port of the MATLAB `CTPUTests` toolbox's FM-OLS machinery, structured so
further estimators and a panel version can be added later.

## Status

- `cpr()`: single-equation cointegrating polynomial regression.
  - Estimator: `"FMOLS"` (fully modified OLS) is implemented.
    `"DOLS"`, `"MOLS"`, `"IMOLS"` are reserved and registered as
    placeholders (see `R/estimators.R`) — adding one is a matter of writing
    a `fit_*_cpr(y, x, orders, w, deter, kernel, bandwidth)` function and
    plugging it into `.cpr_estimators`.
  - Bandwidth selection: `"And91"` (Andrews, 1991), `"AM92"` (Andrews &
    Monahan, 1992, VAR(1) pre-whitened), `"NW"` (Newey & West, 1994), or a
    fixed numeric bandwidth.
  - Kernels: truncated (`"tr"`), Bartlett (`"ba"`), Parzen (`"pa"`), Bohman
    (`"bo"`), Daniell (`"da"`), Quadratic Spectral (`"qs"`).
  - Demeaning of `Delta(x_t)` in the long-run variance step is always
    applied (Wagner & Reichold 2023, Remark 5) and is not a user-facing
    option.
- `ct_test()`: KPSS/Shin-type CT cointegration test for a fitted CPR
  (port of `CT_test.m`). Critical values are bundled only for the
  intercept-only, one-regressor, max-power-2 case (`d = 0, m = 1, p = 2`)
  used in the CEE panel example below; more can be added to
  `.ct_critval_table` in `R/ct-test.R` from the original
  `CTcritval/*.mat` files.
- `pcpr()`: panel cointegrating polynomial regression.
  - `type = "mg"` (mean group, default): calls [`cpr()`] once per
    cross-sectional unit -- literally the same estimation as a standalone
    `cpr()` call, not a parallel reimplementation -- and averages the
    unit-specific coefficients (Pesaran & Smith, 1995, mean-group
    estimator and between-unit inference). `object$unit_fits` holds the
    individual `"cpr"` objects.
  - `type = "pmg"` (panel / pooled mean group, common slope across units)
    is registered as a placeholder (see `.pcpr_types` in `R/pcpr.R`) and
    currently raises an informative error.
  - Requires a **balanced panel** (every unit observed at the same set of
    time points): this guarantees the long-run variance for every unit's
    model is estimated over the same number of time points, so units stay
    comparable before being averaged. The auto-selected bandwidth itself
    is still allowed to differ across units under `type = "mg"` -- that is
    the point of allowing full slope heterogeneity; only the sample length
    feeding into it is held fixed.
- `pu_test()`: Phillips-Ouliaris-type PU cointegration test (port of
  `PU_test.m`). Note this tests the *opposite* null from `ct_test()` (PU:
  H0 = no cointegration; CT: H0 = cointegration), operates directly on the
  raw `y`/`x` series rather than a fitted `cpr()` object's residuals, and
  is a genuinely different statistic (a VAR(1)-based long-run-variance
  ratio, not a KPSS-type partial-sum statistic) -- so its results
  shouldn't be expected to mirror `ct_test()`'s automatically. Same
  critical-value bundling approach and caveat as `ct_test()`: only
  `d = 0, m = 1, p = 2` is available so far, in `.pu_critval_table` in
  `R/pu-test.R`.
- Homogeneity tests and turning-point analysis are not implemented yet.

## Usage

```r
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "estimators.R", "cpr.R", "pcpr.R", "ct-test.R", "pu-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

fit <- cpr(y, x, orders = 2)   # y ~ const + x + x^2, FM-OLS, And91/Bartlett
summary(fit)

fit_mg <- pcpr(y, x, id = country, time = year, orders = 2)   # mean-group panel version
summary(fit_mg)
```

See `examples/example_cpr.R` for a fuller `cpr()` walkthrough (trend
terms, stationary regressors, alternate bandwidths).

`examples/example_panel_cee.R` reproduces the original MATLAB
`FM_OLS_panel.m` / `CT_test.m` country-by-country output (NOIP ~ GNIPC +
GNIPC^2 for 13 Central/Eastern European economies, from
`deJongWagner2022`) using `inst/extdata/cee_panel.csv` (COUNTRY, YEAR,
NOIP, GNIPC extracted from that example's `panel.xlsx`, sheet "CEE").
The R output matches the original MATLAB output to 3 decimal places,
including the CT test decisions at the 10%/5%/1% levels.

`examples/example_pcpr_mg.R` reruns the same CEE panel through
`pcpr(type = "mg")` instead of a manual per-country loop, and reports the
group-mean coefficient across all 13 countries alongside the per-country
estimates that were averaged.

`examples/example_pu_test.R` runs `pu_test()` alongside `ct_test()` for
all 13 CEE countries. In this data, CT never rejects for any country
(consistent with cointegration) while PU rejects only for Romania
(consistent with *no* cointegration for the rest) -- the standard
inconclusive-but-not-contradictory joint outcome in a short (T = 28)
panel, since the two tests have opposite nulls.

## Tests

```
Rscript tests/test-cpr.R
```

Base-R sanity tests (no external package dependencies): coefficient
recovery on simulated data, all valid kernel/bandwidth combinations,
informative errors for invalid combinations and unimplemented estimators,
closed-form checks on the low-level building blocks, a regression test
pinning the CEE panel output to the original MATLAB results, checks
that `pcpr(type = "mg")`'s per-unit fits are bit-for-bit identical to
standalone `cpr()` calls, that it rejects unbalanced panels, that
`type = "pmg"` fails clearly as not yet implemented, and that `pu_test()`'s
bundled critical values match the original `PUcritval/*.mat` table.
