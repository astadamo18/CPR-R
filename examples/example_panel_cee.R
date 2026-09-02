# Reproduces the country-by-country FM-OLS panel example from
# deJongWagner2022 (FM_OLS_panel.m applied to the "CEE" sheet of
# panel.xlsx): NOIP ~ const + GNIPC + GNIPC^2, per country, FM-OLS with
# a Bartlett kernel and Andrews (1991) bandwidth, plus the CT
# cointegration test at the 10%/5%/1% levels.
#
# Run from the package root with: Rscript examples/example_panel_cee.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "cpr.R", "ct-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
countries <- unique(panel$COUNTRY)

cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s %8s %-12s %-12s %-12s\n",
            "COUNTRY", "const", "p-value", "GNIPC", "p-value", "GNIPC^2", "p-value",
            "CT", "10%", "5%", "1%"))

for (cname in countries) {
  sub <- panel[panel$COUNTRY == cname, ]
  sub <- sub[order(sub$YEAR), ]

  y <- sub$NOIP / 1000   # scaled for numerical stability, matching FM_OLS_panel.m
  x <- sub$GNIPC / 1000

  fit <- cpr(y, x, orders = 2, kernel = "ba", bandwidth = "And91")

  ct <- ct_test(fit$fit$residuals, fit$fit$Omega_udotv1, d = 0, m = 1, p = 2)
  decision <- ifelse(ct$reject, "rejection", "no rejection")

  ct_tbl <- fit$coef_table
  cat(sprintf("%-10s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %-12s %-12s %-12s\n",
              cname,
              ct_tbl["const", "Estimate"], ct_tbl["const", "Pr(>|z|)"],
              ct_tbl["x1^1", "Estimate"], ct_tbl["x1^1", "Pr(>|z|)"],
              ct_tbl["x1^2", "Estimate"], ct_tbl["x1^2", "Pr(>|z|)"],
              ct$statistic, decision[1], decision[2], decision[3]))
}
