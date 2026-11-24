# Reference fixtures for StatisticalModels.EmpiricalNull (Test/StatisticalModels/EmpiricalNullTests.cs).
# R 4.6.1, locfdr 1.1.8 (CRAN; Efron, Turnbull and Narasimhan). Run from this folder:
#   Rscript make_empirical_null_fixtures.R
# Writes empirical_null_inputs.tsv, empirical_null_results.tsv and PROVENANCE_empirical_null.txt. R is never run by the
# tests.
#
# Each scenario is one set of z-values. The expected empirical null is locfdr's maximum-likelihood one:
#   locfdr(z, nulltype = 1, plot = 0), the "mlest" row of fp0 (delta, sigma, p0).
# Its first pass is also written, as locfdr's default window gives it: median(z) +- b * (Q3 - Q1) / (2 qnorm(0.75)),
# b = 4.3 * exp(-0.26 * log10(N)). The second pass refits on first delta +- b * first sigma.
# locfdr's own estimate stops its damped Newton iteration when the natural parameters move by less than 1e-5, from the
# start (0, 1), and on some score sets it stops with "ML estimation failed" although the maximum exists. So each scenario
# also carries the exact maximum: the same two windows, each fitted by locfdr:::locmle from the window scores' mean and
# standard deviation with ep = 1e-14 and Jmle = 1000 (mle_*). locfdr_status is locfdr's own outcome, estimated or failed.
suppressPackageStartupMessages(library(locfdr))
set.seed(20261124)

mix <- function(n, parts) {
  # parts: list of c(share, mean, sd); sd 0 means a point mass (a clamped score)
  counts <- round(n * sapply(parts, `[`, 1)); counts[1] <- n - sum(counts[-1])
  unlist(Map(function(k, p) if (p[3] == 0) rep(p[2], k) else rnorm(k, p[2], p[3]), counts, parts))
}
scenarios <- list(
  list(name = "null_200",            z = function() rnorm(200)),
  list(name = "null_2000",           z = function() rnorm(2000)),
  list(name = "null_10000",          z = function() rnorm(10000)),
  list(name = "small_50",            z = function() rnorm(50)),
  list(name = "shifted_narrow_1000", z = function() rnorm(1000, 0.4, 0.9)),
  list(name = "shifted_wide_3000",   z = function() rnorm(3000, -0.2, 1.3)),
  list(name = "one_sided_10pct",     z = function() mix(2000, list(c(0.9, 0.1, 1), c(0.1, 3, 1)))),
  list(name = "one_sided_30pct",     z = function() mix(2500, list(c(0.7, -0.3, 1), c(0.3, 2.5, 1)))),
  list(name = "two_sided_20pct",     z = function() mix(1500, list(c(0.8, 0, 1), c(0.1, -3, 1), c(0.1, 3, 1)))),
  list(name = "clamped_tails",       z = function() mix(500, list(c(0.94, 0, 1), c(0.03, -8, 0), c(0.03, 8, 0)))),
  list(name = "clamped_one_side_30", z = function() mix(2000, list(c(0.7, 0, 1), c(0.3, 8, 0)))),
  list(name = "far_modes_25",        z = function() mix(2000, list(c(0.5, 0, 1), c(0.25, 6, 0.5), c(0.25, -6, 0.5)))),
  list(name = "ties_800",            z = function() round(rnorm(800, 0.2, 1), 2)),
  list(name = "clamped_both_25",     z = function() mix(2000, list(c(0.5, 0, 1), c(0.25, 8, 0), c(0.25, -8, 0))))
)

g <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
inputs <- data.frame(); results <- data.frame()
for (s in scenarios) {
  z <- s$z()
  n <- length(z)
  b <- 4.3 * exp(-0.26 * log(n, 10))
  med <- median(z); sc <- diff(quantile(z)[c(2, 4)]) / (2 * qnorm(0.75))
  fit <- tryCatch(suppressWarnings(locfdr(z, nulltype = 1, plot = 0)), error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    if (!grepl("ML estimation failed|CM and ML Estimation failed", fit)) stop(s$name, ": unexpected locfdr error: ", fit)
    status <- "failed"; est <- c(NA, NA, NA)
  } else {
    status <- "estimated"; est <- fit$fp0["mlest", c("delta", "sigma", "p0")]
  }
  # The exact maximum on locfdr's windows. locmle returns its named estimate vector (del0, sig0, p0, ...) when called
  # without sw or Cov.in.
  exact <- function(centre, half) {
    w <- z[z >= centre - half & z <= centre + half]
    locfdr:::locmle(z, xlim = c(centre, half), d = mean(w), s = sd(w), ep = 1e-14, Jmle = 1000)
  }
  first <- exact(med, b * sc)
  second <- exact(first[1], b * first[2])
  inputs <- rbind(inputs, data.frame(scenario = s$name, z = g(z)))
  results <- rbind(results, data.frame(scenario = s$name, n = n, locfdr_status = status,
    locfdr_delta0 = g(est[1]), locfdr_sigma0 = g(est[2]), locfdr_p0 = g(est[3]),
    first_center = g(med), first_half_width = g(b * sc),
    mle_first_delta0 = g(first[1]), mle_first_sigma0 = g(first[2]),
    mle_delta0 = g(second[1]), mle_sigma0 = g(second[2]), mle_p0 = g(second[3])))
}
write.table(inputs, "empirical_null_inputs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results, "empirical_null_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("locfdr", as.character(packageVersion("locfdr"))),
  "seed 20261124",
  "script make_empirical_null_fixtures.R",
  "call locfdr(z, nulltype = 1, plot = 0), fp0[\"mlest\", c(\"delta\", \"sigma\", \"p0\")] (locfdr_*)",
  "exact locfdr:::locmle(z, xlim = c(centre, half), d = mean(window), s = sd(window), ep = 1e-14, Jmle = 1000) on centre median, half b * IQR / (2 qnorm(0.75)), then on the first estimate +- b * its sigma (mle_*)"
), "PROVENANCE_empirical_null.txt")
