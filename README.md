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
- `pcpr()` (panel version), the PU test, homogeneity tests, and
  turning-point analysis are not implemented yet.

## Usage

```r
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "estimators.R", "cpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

fit <- cpr(y, x, orders = 2)   # y ~ const + x + x^2, FM-OLS, And91/Bartlett
summary(fit)
```

See `examples/example_cpr.R` for a fuller walkthrough (trend terms,
stationary regressors, alternate bandwidths).

`examples/example_panel_cee.R` reproduces the original MATLAB
`FM_OLS_panel.m` / `CT_test.m` country-by-country output (NOIP ~ GNIPC +
GNIPC^2 for 13 Central/Eastern European economies, from
`deJongWagner2022`) using `inst/extdata/cee_panel.csv` (COUNTRY, YEAR,
NOIP, GNIPC extracted from that example's `panel.xlsx`, sheet "CEE").
The R output matches the original MATLAB output to 3 decimal places,
including the CT test decisions at the 10%/5%/1% levels.

## Tests

```
Rscript tests/test-cpr.R
```

Base-R sanity tests (no external package dependencies): coefficient
recovery on simulated data, all valid kernel/bandwidth combinations,
informative errors for invalid combinations and unimplemented estimators,
closed-form checks on the low-level building blocks, and a regression
test pinning the CEE panel output to the original MATLAB results.
