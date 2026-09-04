# Example: turning-point analysis and its plot() methods.
#
# For a quadratic (or cubic) cointegrating polynomial relationship
# (const + b1*x + b2*x^2 [+ b3*x^3]), a "turning point" is where the fitted
# curve's slope in x is zero -- the EKC-style vertex/inflection of the
# income/emissions relationship this toolbox was originally built around
# (de Jong & Wagner, 2016). turning_points() computes it; plot() draws the
# fitted curve with it labeled.
#
# Run from the package root with: Rscript examples/example_turning_points.R
# (writes PNG files into examples/, since this is meant to run headlessly)

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R",
  "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R",
  "turning-points.R", "plot.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
panel$noip1000 <- panel$NOIP / 1000
panel$gnipc1000 <- panel$GNIPC / 1000

## ---- Single case: Czechia ----
cz <- panel[panel$COUNTRY == "Czechia", ]
cz <- cz[order(cz$YEAR), ]
fit_cz <- cpr(noip1000 ~ gnipc1000, data = cz, orders = 2, kernel = "ba", bandwidth = "And91")

cat("=== turning_points(fit_cz) ===\n")
print(turning_points(fit_cz))

## The curve plotted is const + b1*x + b2*x^2 -- the constant is always
## added in (here 13.327, from fit_cz$coefficients["const"]) even though
## the turning point's x-location (-b1/(2*b2)) never depends on it; leaving
## it out would draw a curve at the wrong *level* while still marking the
## turning point at the right x.
grDevices::png("examples/turning_points_single.png", width = 800, height = 600)
plot(fit_cz)
grDevices::dev.off()
cat("Wrote examples/turning_points_single.png\n\n")

## ---- Panel case: mean group, averaged turning points ----
fit_mg <- pcpr(noip1000 ~ gnipc1000, data = panel, id = COUNTRY, time = YEAR,
               orders = 2, kernel = "ba", bandwidth = "And91", type = "mg")

cat("=== turning_points(fit_mg) ===\n")
print(turning_points(fit_mg))

## Each of the 13 countries gets its own turning point from its own
## (own-constant-included) curve; those falling inside that country's own
## observed x-range are averaged by type (here: all "minimum") to get the
## panel's single reported turning point. Its y-value is read off the
## group-mean curve -- pcpr(type = "mg")'s own averaged coefficients, which
## already average the constant across units the same way it averages the
## slopes, so the plotted level is the "proper" (group-mean) one, not zero
## and not an arbitrary single country's.
grDevices::png("examples/turning_points_mg.png", width = 800, height = 600)
plot(fit_mg)
grDevices::dev.off()
cat("Wrote examples/turning_points_mg.png\n\n")

## ---- Panel case: pooled (common slope) ----
fit_pmg <- pcpr(noip1000 ~ gnipc1000, data = panel, id = COUNTRY, time = YEAR,
                orders = 2, kernel = "ba", bandwidth = "And91", type = "pmg")

cat("=== turning_points(fit_pmg) ===\n")
print(turning_points(fit_pmg))
## Empty here: the pooled model's single common-slope curve has its vertex
## outside the observed GNIPC range for every country, so there is no
## *interior* turning point to report -- a real finding (this data does not
## support a common EKC-style turning point under full slope pooling), not
## a bug. The pooled model also has no single estimated constant (fixed
## effects absorb it); the curve/label instead use the average, across
## countries, of each one's own implied fixed effect -- see the file-level
## comment in R/turning-points.R.
grDevices::png("examples/turning_points_pmg.png", width = 800, height = 600)
plot(fit_pmg)
grDevices::dev.off()
cat("Wrote examples/turning_points_pmg.png\n")
