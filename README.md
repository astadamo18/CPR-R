# CPR

Cointegrating Polynomial Regressions in R, following Wagner and co-authors
(e.g. Wagner; Wagner & Reichold; de Jong & Wagner; Veldhuis & Wagner). This is
an R port of the MATLAB `CTPUTests` toolbox's FM-OLS machinery, structured so
further estimators and a panel version can be added later.

## Status

- `cpr()`: single-equation cointegrating polynomial regression.
  - Estimator: `"FMOLS"` (fully modified OLS, `R/fmols.R`) and `"DOLS"`
    (dynamic OLS, `R/dols.R`, port of `MonitoringCPR_MatlabCode/MonitoringCPR/
    DOLS_CPR.m` and its `GenLeadLag.m` helper) are implemented. `"MOLS"`
    remains a placeholder (see `R/estimators.R`) — no standalone
    single-equation MATLAB source for it was found in the original
    toolbox, only a panel-embedded bias correction inside `PanelEKC_*.m`,
    already ported as `pcpr(type = "pmg")`'s `beta_Mod`. `"IMOLS"` also
    remains a placeholder, but real source exists to port later
    (`IM-SCMPR-ExemplaryCode/im_scmpr.m` / `MonitoringCPR_MatlabCode/
    MonitoringCPR/IMOLS_NL.m`, Vogelsang & Wagner's Integrated Modified
    OLS). Adding a new estimator is a matter of writing a
    `fit_*_cpr(y, x, orders, w, deter, kernel, bandwidth, n_lag, n_lead)`
    function and plugging it into `.cpr_estimators`.
    DOLS augments the polynomial regression with leads and lags of
    `Delta(x)` (`n_lag`/`n_lead`, both default `0`) and gets consistency
    from that augmentation via plain OLS, with HAC standard errors
    afterward -- a genuinely different mechanism from FM-OLS's
    Schur-complement correction, not a variant of it. It does not support
    `w`, and (unlike FM-OLS) does not drop the first observation when
    `n_lag = n_lead = 0`.
  - Bandwidth selection: `"And91"` (Andrews, 1991), `"AM92"` (Andrews &
    Monahan, 1992, VAR(1) pre-whitened), `"NW"` (Newey & West, 1994), or a
    fixed numeric bandwidth.
  - Kernels: truncated (`"tr"`), Bartlett (`"ba"`), Parzen (`"pa"`), Bohman
    (`"bo"`), Daniell (`"da"`), Quadratic Spectral (`"qs"`).
  - Demeaning of `Delta(x_t)` in the long-run variance step is always
    applied (Wagner & Reichold 2023, Remark 5) and is not a user-facing
    option.
  - Accepts an `lm()`-like formula/`data` interface alongside the original
    vectors: `cpr(y ~ x1 + x2, data = df, orders = 2)`, with `w`/`deter`
    also acceptable as one-sided formulas (`w = ~ w1 + w2`) against the
    same `data`. Purely additive -- the original `cpr(y, x, orders, ...)`
    vector form still works unchanged, and is what the formula form
    resolves to internally. See the file-level comment in
    `R/formula-data.R` for why a formula's right-hand side must name
    columns of `data` verbatim (no `log(x1)`-style transformed terms).
    The returned fit also carries the resolved `y`/`x` (`$y`/`$x`), which
    `ct_test()`/`pu_test()` reuse when dispatched on the fit directly.
- `ct_test()`: KPSS/Shin-type CT cointegration test for a fitted CPR
  (port of `CT_test.m`). An S3 generic: call it directly on a fitted
  `cpr` object -- `ct_test(fit)` -- and `uplus`/`omega`/`m`/`p`/`d` are all
  read straight off the fit (`d` is inferred from `deter`: none/const-only/
  const+trend map to `d = -1/0/1`; pass `d` explicitly to override, e.g.
  for a custom `deter` the inference can't classify), with the matching
  critical value table picked automatically -- the full grid of 48 tables
  bundled (`d` in `{-1, 0, 1}`, `m` and `p` in `{1, 2, 3, 4}`, extracted
  from every `CTcritval/*.mat` file in the original toolbox, not just the
  one the CEE example needs) covers every combination the original
  toolbox itself tabulates. Only supports `estimator = "FMOLS"` fits: the
  bundled critical values were Monte Carlo simulated specifically for
  FM-OLS residuals (`CT_test.m`'s own docstring: "uplus...FM-OLS
  residuals"), and nothing establishes that a different estimator's
  residuals (e.g. DOLS's) share that null distribution -- `ct_test()` on a
  non-FMOLS fit errors rather than silently reusing the table. (`pu_test()`,
  below, has no such restriction -- it never touches a fit's
  estimator-specific residuals at all.) `ct_test(uplus, omega, d, m, p)`
  still works too, for the raw ingredients.

  The result has a `print()` method showing the test statistic, critical
  values (10%/5%/1%), the reject/do-not-reject decision at each, and
  H0/H1. No p-value or significance stars for now: the null distribution
  is only known via the 9 tabulated percentiles per table, which isn't
  enough points to interpolate a p-value with any real precision.
- `pcpr()`: panel cointegrating polynomial regression.
  - `type = "mg"` (mean group, default): calls [`cpr()`] once per
    cross-sectional unit -- literally the same estimation as a standalone
    `cpr()` call, not a parallel reimplementation -- and averages the
    unit-specific coefficients (Pesaran & Smith, 1995, mean-group
    estimator and between-unit inference). `object$unit_fits` holds the
    individual `"cpr"` objects.
  - `type = "pmg"` (pooled panel, de Jong & Wagner 2016; port of
    `deJongWagner2022/PanelEKC_indiv_eff_only.m` and `PanelEKC_two_eff.m`):
    a genuinely different estimator from `"mg"` -- a single common slope
    shared by all units (removed via a within/LSDV transformation,
    `effects = "oneway"` for individual fixed effects only or `"twoway"`
    for individual + time fixed effects), with per-unit long-run variances
    pooled (averaged) into one bias correction applied uniformly. Produces
    three point estimates (`beta_lsdv`, `beta_Mod`, `beta_FM`; the
    coefficient table reports `beta_FM`) and three VCV matrices, in
    `object$unit_fits`. See the file-level comment in `R/pooled-panel.R`
    for the full derivation. Restricted to a single integrated regressor
    with `orders` exactly `2` or `3` (the theoretical bias-correction
    matrices are only tabulated for those cases in the original source),
    and does not support `w`. Unlike `cpr()`/`"mg"`, it does not drop the
    first time observation.
  - Requires a **balanced panel** (every unit observed at the same set of
    time points): this guarantees the long-run variance for every unit's
    model is estimated over the same number of time points, so units stay
    comparable before being averaged or pooled. The auto-selected
    bandwidth itself is still allowed to differ across units -- that is
    the point of allowing (`"mg"`) or accounting for (`"pmg"`, per-unit
    Omega_i) heterogeneity; only the sample length feeding into it is held
    fixed.
  - Accepts the same formula/`data` interface as `cpr()`. `id`/`time`
    take a *bare* column name too, `lm()`-like -- no quotes needed:
    `pcpr(y ~ x1, data = df, id = country, time = year, orders = 2)`
    (a quoted `id = "country"`, or a raw vector `id = df$country`, both
    still work). `w`/`deter` accept column-name strings/vectors or
    one-sided formulas against `data` the same way `cpr()`'s do.
- `pu_test()`: Phillips-Ouliaris-type PU cointegration test (port of
  `PU_test.m`), rewritten the same way `ct_test()` was: an S3 generic with
  the full 48-table `(d, m, p)` critical-value grid bundled (from every
  `PUcritval/*.mat` file in the original toolbox) and a matching
  `print()` method (statistic, critical values, decisions, H0/H1; no
  p-value/stars, same reasoning as `ct_test()`'s). `pu_test(fit)`
  dispatches straight off a fitted `cpr` object -- reusing its stored
  `$y`/`$x` and `kernel`/`bandwidth`, and inferring `d` the same way
  `ct_test.cpr()` does -- for a single-regressor fit with sequential
  powers `1:p`; anything else (multiple regressors, non-sequential
  powers) needs `pu_test.default(y, x, d, m, orders, kernel, bandwidth)`
  called directly. Note PU tests the *opposite* null from `ct_test()`
  (PU: H0 = no cointegration; CT: H0 = cointegration) and is a genuinely
  different statistic (its own internal plain-OLS regression plus a local
  VAR(1)-with-deterministics on the stacked `[y, x]` system, not a
  KPSS-type partial-sum statistic built from FM-OLS residuals) -- so its
  results shouldn't be expected to mirror `ct_test()`'s automatically, and
  it works off any `cpr` fit regardless of `estimator` (it never touches
  the fit's residuals).
- `turning_points()` / `plot()`: EKC-style turning-point analysis (port of
  the analysis underlying `deJongWagner2022`'s income/emissions curve), for
  a single-regressor fit's quadratic/cubic (or higher) relationship in `x`.
  A turning point is where the fitted curve's slope in `x` is zero; its
  *location* only depends on the slope coefficients, but the curve's
  *level* (and so the plotted/labeled turning-point value) also needs the
  constant, which is always included even though it never moves the
  turning point's x-position -- see the file-level comment in
  `R/turning-points.R`.
  - `turning_points(fit)`: an S3 generic. For a `cpr` fit, returns a data
    frame of `x`/`y`/`type` (`"maximum"`/`"minimum"`/`"inflection"`),
    restricted by default to turning points inside the observed range of
    `x` (interior turning points only -- pass `x_range = NULL` to keep
    extrapolated roots too). Zero rows for a purely linear fit. Only
    supports a single integrated regressor.
  - For a `pcpr` fit: `type = "mg"` computes each unit's own turning
    point first (mean-group philosophy: average a nonlinear function of
    the per-unit estimates, the same way the coefficients themselves are
    averaged), then averages by type across the units that have one,
    reporting the count (`n_units`) and using the panel's own group-mean
    curve (constant included) to compute the labeled `y`. `type = "pmg"`
    has a single common slope, so at most one turning point per type; the
    pooled model has no single estimated constant (fixed effects absorb
    it), so its curve uses the average, across units, of each one's own
    implied fixed effect instead -- reconstructed from each unit's own
    *raw* (not demeaned) data, so it is a real, data-scale level rather
    than an artifact of the within-transformed estimation; verified
    against an independent dummy-variable regression (see the `@details`
    on `turning_points.pcpr()` in `R/turning-points.R` and the
    corresponding test in `tests/test-cpr.R`).
  - `plot(fit)`: draws the fitted curve (`cpr`: with the observed data
    scatter; `pcpr(type = "mg")`: the group-mean curve, with each unit's
    own curve shown faintly for context; `pcpr(type = "pmg")`: the single
    pooled curve) with turning point(s) marked and labeled, and invisibly
    returns the same data `turning_points()` would.
- Homogeneity tests are not implemented yet.

## Installation

```r
# install.packages("remotes")  # if you don't have it yet
remotes::install_github("astadamo18/CPR-R", ref = "claude/new-session-9n2jmo")
library(CPR)
```

(`devtools::install_github()` works identically -- it wraps `remotes` for
this.) The `ref` is required because everything currently lives on that
branch, not `main`. There is no `man/` directory yet (the roxygen-style
comments in `R/*.R` were never run through `roxygen2`), so `library(CPR)`
prints a harmless "No man pages found" note and `?cpr` won't return
anything -- read the source or the examples below instead.

Alternatively, without installing, source the files directly from a local
clone (also how the test suite and example scripts run):

```r
source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R",
  "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R",
  "turning-points.R", "plot.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))
```

## Usage

```r
fit <- cpr(y, x, orders = 2)   # y ~ const + x + x^2, FM-OLS, And91/Bartlett
summary(fit)

fit <- cpr(y ~ x, data = df, orders = 2)   # equivalent, lm()-like formula/data form
summary(fit)

fit_dols <- cpr(y, x, orders = 2, estimator = "DOLS", n_lag = 1, n_lead = 1)
summary(fit_dols)

fit_mg <- pcpr(y, x, id = country, time = year, orders = 2)   # mean-group panel version
summary(fit_mg)

fit_mg <- pcpr(y ~ x, data = df, id = "country", time = "year", orders = 2)   # equivalent

fit_pmg <- pcpr(y, x, id = country, time = year, orders = 2, type = "pmg")   # pooled, common slope
summary(fit_pmg)

ct_test(fit)   # dispatches on the fitted cpr object either way
pu_test(fit)

turning_points(fit)   # EKC-style turning point(s): x, y, type ("maximum"/"minimum")
plot(fit)              # ... plus the fitted curve, labeled at the turning point(s)

turning_points(fit_mg)   # panel: averaged (by type) across units' own turning points
plot(fit_mg)              # ... group-mean curve, each unit's own curve shown faintly
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

`examples/example_ct_test_print.R` shows `ct_test()`'s print output
(statistic, critical values, decisions, H0/H1) on Czechia for `d = 0`,
`d = 1` (with a trend), and `p = 3` (cubic) -- each automatically pulling a genuinely different
critical value table from the bundled grid.

`examples/example_pu_test.R` runs `pu_test(fit)` alongside `ct_test(fit)`
for all 13 CEE countries. In this data, CT never rejects for any country
(consistent with cointegration) while PU rejects only for Romania
(consistent with *no* cointegration for the rest) -- the standard
inconclusive-but-not-contradictory joint outcome in a short (T = 28)
panel, since the two tests have opposite nulls.

`examples/example_formula_data.R` walks through the `lm()`-like
formula/`data` interface for both `cpr()` and `pcpr()` (including `w` as
a one-sided formula and `pcpr()`'s `id`/`time` as *bare* column names),
confirming it gives identical fits to the vector interface, then runs
`ct_test()`/`pu_test()` off the resulting fit.

`examples/example_turning_points.R` computes and plots EKC-style turning
points for Czechia alone, the mean-group panel, and the pooled panel --
including the pooled case's honest empty result (its common-slope curve's
vertex falls outside every country's observed GNIPC range in this data, so
there is no interior turning point to report). Writes PNGs into `examples/`
(not tracked by git; see `.gitignore`) since it's meant to run headlessly.

`examples/example_pcpr_pmg.R` fits the pooled panel model (`type = "pmg"`)
on the same CEE panel, both `oneway` and `twoway`, and compares its
(single, common) slope against `pcpr(type = "mg")`'s group-mean slope --
directionally consistent (negative linear term, small positive quadratic
term) in this data, as expected since the two estimators target the same
underlying relationship under different homogeneity assumptions.

`examples/example_dols.R` compares `cpr(estimator = "DOLS")` against
`"FMOLS"` on Czechia: with `n_lag = n_lead = 0` DOLS is plain OLS on the
full (untruncated) sample; with leads/lags added the point estimates move
but stay in the same neighborhood as FM-OLS (const ~13-16, GNIPC ~-1.2 to
-1.4, GNIPC^2 ~0.013-0.017, all significant either way).

## Tests

```
Rscript tests/test-cpr.R
```

Base-R sanity tests (no external package dependencies): coefficient
recovery on simulated data, all valid kernel/bandwidth combinations,
informative errors for invalid combinations and unimplemented estimators,
that DOLS with `n_lag = n_lead = 0` matches plain OLS on the full
(untruncated) sample and that its truncation with leads/lags is exactly
right, that DOLS rejects `w`, closed-form checks on `gen_lead_lag()` and
the other low-level building blocks, a regression test pinning the CEE
panel output to the original MATLAB results, checks that
`pcpr(type = "mg")`'s per-unit fits are bit-for-bit identical to
standalone `cpr()` calls and that it rejects unbalanced panels,
that `pcpr(type = "pmg")` runs for oneway/twoway effects and q = 2/3 and
rejects its unsupported cases (q outside `{2, 3}`, `w`) with clear errors,
a qualitative consistency check between `pmg` (N=1) and standalone `cpr()`,
that all 48 bundled CT and PU critical-value tables load and are monotone,
that `ct_test()`/`pu_test()` dispatch correctly on a fitted `cpr` object
with their `print()` methods showing the expected sections, that
`ct_test()` rejects a non-`"FMOLS"` fit with a clear error (its critical
values are only valid for FM-OLS residuals) while `pu_test()` runs fine on
the same `"DOLS"` fit (it never touches estimator-specific residuals), and
that `cpr()`/`pcpr()`'s formula/`data`
interface (including `w`/`deter` as formulas and `pcpr()`'s `id`/`time`
as column-name strings) gives results identical to the vector interface
and errors clearly when `data` is missing or a named column isn't found,
that `turning_points()` matches the closed-form vertex of a quadratic fit
and handles the linear (no turning point) and multi-regressor (error)
edge cases, that the panel `"mg"` average is exactly the mean (by type) of
the per-unit turning points recomputed independently, that the panel
`"pmg"` case correctly restricts to the observed range, and that the
`plot()` methods run without error and return the same data
`turning_points()` does.

### A cross-platform bug this port found and fixed

`lr_weights()`/`lr_var()` had a latent bug inherited from a MATLAB-to-R
translation subtlety: MATLAB's colon operator returns an *empty* range for
reversed bounds (e.g. `76:28` is `[]`), so an auto-selected bandwidth
exceeding `T-1` is harmless there -- out-of-range lags just contribute
nothing. R's `:` instead returns a *descending* sequence (`76:28` is
`76, 75, ..., 28`), which crashed downstream indexing in `lr_var()`. This
surfaces whenever Andrews (1991) picks a bandwidth larger than `T-1`,
which happens easily for a short, persistent series (it did for Hungary,
Belarus, and Ukraine with T=28 while building `pcpr(type = "pmg")`).
Fixed by capping the weight vector's `upper` index at `T-1` in
`lr_weights()` for every kernel (mathematically the correct fix, not just
a defensive patch: lags beyond `T-1` have no data to form an
autocovariance from). Covered by the `pmg` tests above.
