# Reference fixtures for StatisticalModels.EmpiricalNull (Test/StatisticalModels/EmpiricalNullTests.cs).
# R 4.6.1, locfdr 1.1.8 (CRAN; Efron, Turnbull and Narasimhan). Run from this folder:
#   Rscript make_empirical_null_fixtures.R
# Writes empirical_null_results.tsv and PROVENANCE_empirical_null.txt. R is never run by the tests.
#
# Each scenario's z-values are built from its recipe, which the test builds the same way, so no z-value is stored:
#   parts separated by ";", each "n,mean,sd": mean + sd * qnorm((i - 0.5) / n) for i = 1..n, or n copies of mean when
#   sd is 0 (a clamped score); an optional "|k" repeats every value k times, which makes ties.
# The expected empirical null is locfdr's maximum-likelihood one: locfdr(z, nulltype = 1, plot = 0), the "mlest" row of
# fp0 (delta, sigma, p0), written as locfdr_*. Its first window is locfdr's default: median(z) +- b * (Q3 - Q1) /
# (2 qnorm(0.75)), b = 4.3 * exp(-0.26 * log10(N)); the second refits on the first delta +- b * the first sigma.
# locfdr stops its damped Newton iteration when the natural parameters move by less than 1e-5, from the start (0, 1),
# and on some score sets it stops with "ML estimation failed" although the maximum exists. So each scenario also carries
# the exact maximum: the same two windows, each fitted by locfdr:::locmle from the window scores' mean and standard
# deviation with ep = 1e-14 and Jmle = 1000 (mle_*). locfdr_status is locfdr's own outcome, estimated or failed.
suppressPackageStartupMessages(library(locfdr))

scenarios <- c(
  null_200            = "200,0,1",
  null_2000           = "2000,0,1",
  null_10000          = "10000,0,1",
  small_50            = "50,0,1",
  shifted_narrow_1000 = "1000,0.4,0.9",
  shifted_wide_1500   = "1500,-0.2,1.3",
  one_sided_10pct     = "1800,0.1,1;200,3,1",
  one_sided_30pct     = "1750,-0.3,1;750,2.5,1",
  two_sided_20pct     = "1200,0,1;150,-3,1;150,3,1",
  clamped_tails       = "470,0,1;15,-8,0;15,8,0",
  clamped_one_side_30 = "1400,0,1;600,8,0",
  far_modes_25        = "1000,0,1;500,6,0.5;500,-6,0.5",
  ties_800            = "200,0.2,1|4",
  clamped_both_25     = "1000,0,1;500,8,0;500,-8,0"
)

build <- function(recipe) {
  r <- strsplit(recipe, "|", fixed = TRUE)[[1]]
  each <- if (length(r) > 1) as.integer(r[2]) else 1L
  z <- unlist(lapply(strsplit(r[1], ";", fixed = TRUE)[[1]], function(p) {
    v <- as.numeric(strsplit(p, ",", fixed = TRUE)[[1]])
    n <- v[1]
    if (v[3] == 0) rep(v[2], n) else v[2] + v[3] * qnorm((seq_len(n) - 0.5) / n)
  }))
  rep(z, each = each)
}

g <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
results <- data.frame()
for (name in names(scenarios)) {
  z <- build(scenarios[[name]])
  n <- length(z)
  b <- 4.3 * exp(-0.26 * log(n, 10))
  med <- median(z); sc <- diff(quantile(z)[c(2, 4)]) / (2 * qnorm(0.75))
  fit <- tryCatch(suppressWarnings(locfdr(z, nulltype = 1, plot = 0)), error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    if (!grepl("ML estimation failed|CM and ML Estimation failed", fit)) stop(name, ": unexpected locfdr error: ", fit)
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
  results <- rbind(results, data.frame(scenario = name, recipe = scenarios[[name]], n = n, locfdr_status = status,
    locfdr_delta0 = g(est[1]), locfdr_sigma0 = g(est[2]), locfdr_p0 = g(est[3]),
    first_center = g(med), first_half_width = g(b * sc),
    mle_first_delta0 = g(first[1]), mle_first_sigma0 = g(first[2]),
    mle_delta0 = g(second[1]), mle_sigma0 = g(second[2]), mle_p0 = g(second[3])))
}
write.table(results, "empirical_null_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("locfdr", as.character(packageVersion("locfdr"))),
  "inputs built from each scenario's recipe (no random draw): parts n,mean,sd as mean + sd * qnorm((i - 0.5) / n), sd 0 a point mass, |k repeats each value k times",
  "script make_empirical_null_fixtures.R",
  "call locfdr(z, nulltype = 1, plot = 0), fp0[\"mlest\", c(\"delta\", \"sigma\", \"p0\")] (locfdr_*)",
  "exact locfdr:::locmle(z, xlim = c(centre, half), d = mean(window), s = sd(window), ep = 1e-14, Jmle = 1000) on centre median, half b * IQR / (2 qnorm(0.75)), then on the first estimate +- b * its sigma (mle_*)"
), "PROVENANCE_empirical_null.txt")
