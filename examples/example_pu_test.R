# Example: PU cointegration test, alongside CT, on the CEE panel.
#
# CT (ct_test()) and PU (pu_test()) test *opposite* nulls:
#   - CT: H0 = cointegration (reject => evidence AGAINST cointegration)
#   - PU: H0 = no cointegration (reject => evidence FOR cointegration)
# so agreement between them is not guaranteed, especially in short panels
# (T = 28 here). The standard confirmatory reading (Wagner et al.) is:
# strong evidence for cointegration needs CT to *not* reject *and* PU to
# reject; when both fail to reject (as for most countries below), the
# joint result is inconclusive rather than contradictory.
#
# pu_test(), like ct_test(), is now an S3 generic with a print() method and
# the full (d, m, p) critical-value grid bundled -- pu_test(fit) dispatches
# straight off a fitted cpr object (see also examples/example_ct_test_print.R
# for the print output on its own).
#
# Run from the package root with: Rscript examples/example_pu_test.R

source_order <- c(
  "lr-weights.R", "lr-var.R", "bandwidth.R", "prewhiten.R", "poly-terms.R",
  "fmols.R", "dols.R", "estimators.R", "formula-data.R", "cpr.R", "pooled-panel.R", "pcpr.R", "ct-test.R", "pu-test.R", "methods.R"
)
invisible(lapply(file.path("R", source_order), source))

panel <- read.csv("inst/extdata/cee_panel.csv", stringsAsFactors = FALSE)
countries <- unique(panel$COUNTRY)

cat(sprintf("%-10s %8s %-12s %-12s %-12s | %8s %-12s %-12s %-12s\n",
            "COUNTRY", "CT", "10%", "5%", "1%", "PU", "10%", "5%", "1%"))

for (cname in countries) {
  sub <- panel[panel$COUNTRY == cname, ]
  sub <- sub[order(sub$YEAR), ]
  y <- sub$NOIP / 1000
  x <- sub$GNIPC / 1000

  fit <- cpr(y, x, orders = 2, kernel = "ba", bandwidth = "And91")
  ct <- ct_test(fit)
  ct_dec <- ifelse(ct$reject, "rejection", "no rejection")

  pu <- pu_test(fit)
  pu_dec <- ifelse(pu$reject, "rejection", "no rejection")

  cat(sprintf("%-10s %8.3f %-12s %-12s %-12s | %8.3f %-12s %-12s %-12s\n",
              cname, ct$statistic, ct_dec[1], ct_dec[2], ct_dec[3],
              pu$statistic, pu_dec[1], pu_dec[2], pu_dec[3]))
}
