# Reference values for NormalizationCheck's empirical null (Test/Quantification/Differential/NormalizationCheckTests.cs).
# R 4.6.1, locfdr 1.1.8 (Efron, Turnbull and Narasimhan). Run from this folder: Rscript make_normalization_check_fixtures.R
# Writes normalization_check_null.tsv and PROVENANCE_normalization_check.txt. R is never run by the tests.
#
# No z-value is stored. Each scenario's are built from its recipe, as the test builds them: parts "n,mean,sd" separated by
# ";", each mean + sd * qnorm((i - 0.5) / n) for i = 1..n, or n copies of mean when sd is 0; "|k" repeats each value k times.
# locfdr_*: locfdr(z, nulltype = 1, plot = 0)'s fp0["mlest", ] (NA where it stops with "ML estimation failed"). It stops
# its Newton iteration about 1e-5 short, and on some sets fails though the maximum exists, so mle_* is the exact maximum on
# the same two windows (median +- b * IQR / (2 qnorm(0.75)), b = 4.3 * exp(-0.26 * log10(N)); then delta +- b * sigma),
# from locfdr:::locmle run from the window's mean and sd with ep = 1e-14.
suppressPackageStartupMessages(library(locfdr))
scenarios <- c(
  null_200 = "200,0,1", null_2000 = "2000,0,1", null_10000 = "10000,0,1", small_50 = "50,0,1",
  shifted_narrow_1000 = "1000,0.4,0.9", shifted_wide_1500 = "1500,-0.2,1.3",
  one_sided_10pct = "1800,0.1,1;200,3,1", one_sided_30pct = "1750,-0.3,1;750,2.5,1",
  two_sided_20pct = "1200,0,1;150,-3,1;150,3,1", clamped_tails = "470,0,1;15,-8,0;15,8,0",
  clamped_one_side_30 = "1400,0,1;600,8,0", far_modes_25 = "1000,0,1;500,6,0.5;500,-6,0.5",
  ties_800 = "200,0.2,1|4", clamped_both_25 = "1000,0,1;500,8,0;500,-8,0")
build <- function(recipe) {
  r <- strsplit(recipe, "|", fixed = TRUE)[[1]]
  z <- unlist(lapply(strsplit(r[1], ";", fixed = TRUE)[[1]], function(p) {
    v <- as.numeric(strsplit(p, ",", fixed = TRUE)[[1]])
    if (v[3] == 0) rep(v[2], v[1]) else v[2] + v[3] * qnorm((seq_len(v[1]) - 0.5) / v[1])
  }))
  rep(z, each = if (length(r) > 1) as.integer(r[2]) else 1L)
}
g <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))
rows <- lapply(names(scenarios), function(name) {
  z <- build(scenarios[[name]])
  b <- 4.3 * exp(-0.26 * log(length(z), 10))
  fit <- tryCatch(suppressWarnings(locfdr(z, nulltype = 1, plot = 0))$fp0["mlest", c("delta", "sigma")],
                  error = function(e) if (grepl("ML estimation failed", conditionMessage(e))) c(NA, NA) else stop(e))
  exact <- function(centre, half) {
    w <- z[z >= centre - half & z <= centre + half]
    locfdr:::locmle(z, xlim = c(centre, half), d = mean(w), s = sd(w), ep = 1e-14, Jmle = 1000)
  }
  first <- exact(median(z), b * diff(quantile(z)[c(2, 4)]) / (2 * qnorm(0.75)))
  second <- exact(first[1], b * first[2])
  data.frame(scenario = name, recipe = scenarios[[name]], locfdr_delta0 = g(fit[1]), locfdr_sigma0 = g(fit[2]),
             mle_delta0 = g(second[1]), mle_sigma0 = g(second[2]))
})
write.table(do.call(rbind, rows), "normalization_check_null.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("R", R.version.string),
             paste("locfdr", as.character(packageVersion("locfdr"))), "script make_normalization_check_fixtures.R"),
           "PROVENANCE_normalization_check.txt")
