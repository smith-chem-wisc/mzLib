# Reference outputs for RandomEffectsMeta.Pool with RandomEffectsInference.HartungKnapp, from metafor.
#
# Run ONCE, by hand, and the outputs checked in beside this script. R is never a build, runtime or test
# dependency of mzLib: the tests read the frozen TSVs.
#
#   Rscript make_hartung_knapp_fixtures.R <out_dir>
#
# What is compared: metafor rma(yi, sei = se, method = "DL", test = "knha"), the DerSimonian-Laird estimate and tau2
# with the Hartung-Knapp standard error, t(k - 1) interval and two-sided p, at levels 0.95 and 0.90 (the set test's
# interval). The cases are k = 2, 3 and 5 (the numbers of datasets a set hypothesis pools), one with tau2 = 0 where the
# Hartung-Knapp interval is narrower than DerSimonian-Laird's, one with every estimate equal (standard error 0), and two
# larger ones.

suppressPackageStartupMessages(library(metafor))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
fmt <- function(df) {                                  # full precision; NA written as NaN for .NET
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20261114)
cases <- list(
  two_datasets          = list(yi = c(-0.40, 0.10), sei = c(0.20, 0.30)),
  three_datasets        = list(yi = c(0.62, 0.35, 0.91), sei = c(0.21, 0.18, 0.33)),
  five_datasets         = list(yi = c(0.10, 0.60, -0.20, 0.45, 0.05), sei = c(0.10, 0.15, 0.12, 0.20, 0.08)),
  agreeing_three        = list(yi = c(0.300, 0.302, 0.299), sei = c(0.10, 0.20, 0.15)),
  equal_estimates       = list(yi = c(0.25, 0.25, 0.25), sei = c(0.10, 0.20, 0.15)),
  heterogeneous_six     = list(yi = c(0.10, 0.60, -0.20, 0.45, 0.05, 0.80), sei = c(0.10, 0.15, 0.12, 0.20, 0.08, 0.25)),
  random_eight          = list(yi = rnorm(8, 0.2, 0.3), sei = runif(8, 0.05, 0.3))
)
inputs <- list(); results <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  for (i in seq_along(cs$yi)) inputs[[length(inputs) + 1]] <- data.frame(case = nm, yi = cs$yi[i], sei = cs$sei[i])
  for (level in c(95, 90)) {
    r <- rma(yi = cs$yi, sei = cs$sei, method = "DL", test = "knha", level = level)
    results[[length(results) + 1]] <- data.frame(case = nm, level = level / 100, k = r$k, estimate = as.numeric(r$beta),
      se = r$se, df = r$dfs, tau2 = r$tau2, ci_lb = r$ci.lb, ci_ub = r$ci.ub, pval = r$pval)
  }
}
w(fmt(do.call(rbind, inputs)), "hk_inputs.tsv")
w(fmt(do.call(rbind, results)), "hk_results.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("metafor", as.character(packageVersion("metafor"))),
  "seed 20261114",
  "script make_hartung_knapp_fixtures.R",
  "call rma(yi, sei, method = \"DL\", test = \"knha\", level = 95 and 90)"
), file.path(out, "PROVENANCE_hartung_knapp.txt"))
