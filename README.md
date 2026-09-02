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
- `pcpr()` (panel version) and the PU/CT cointegration tests, homogeneity
  tests, and turning-point analysis are not implemented yet.

## Usage

```r
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "estimators.R", "cpr.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

fit <- cpr(y, x, orders = 2)   # y ~ const + x + x^2, FM-OLS, And91/Bartlett
summary(fit)
```

See `examples/example_cpr.R` for a fuller walkthrough (trend terms,
stationary regressors, alternate bandwidths).

## Tests

```
Rscript tests/test-cpr.R
```

Base-R sanity tests (no external package dependencies): coefficient
recovery on simulated data, all valid kernel/bandwidth combinations,
informative errors for invalid combinations and unimplemented estimators,
and closed-form checks on the low-level building blocks.
