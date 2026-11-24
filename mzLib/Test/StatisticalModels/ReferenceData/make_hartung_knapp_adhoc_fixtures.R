# Reference outputs for RandomEffectsMeta.Pool with RandomEffectsInference.HartungKnappAdHoc, from metafor.
#
# Run ONCE, by hand, and the outputs checked in beside this script. R is never a build, runtime or test
# dependency of mzLib: the tests read the frozen TSVs.
#
#   Rscript make_hartung_knapp_adhoc_fixtures.R <out_dir>
#
# What is compared: metafor rma(yi, sei = se, method = "DL", test = "adhoc"), the DerSimonian-Laird estimate and tau2
# with the Hartung-Knapp standard error under Knapp and Hartung's (2003) ad hoc correction (never below the
# DerSimonian-Laird standard error), the t(k - 1) interval and two-sided p, at levels 0.95 and 0.90. The cases are those
# of make_hartung_knapp_fixtures.R (same values and seed), where the correction binds in some and not in others, and
# two more:
#   near_identical_two: two studies that agree to 0.001 log2, each with se 0.3. Plain Hartung-Knapp gives se 0.0005;
#                       the correction gives the DerSimonian-Laird se, sqrt(0.09 / 2).
#   disagreeing_two:    two studies far apart. With k = 2 and tau2 > 0 the two standard errors are equal, so the
#                       correction changes nothing. (three_datasets, five_datasets, heterogeneous_six and random_eight
#                       are the cases where Hartung-Knapp is the larger; agreeing_three and equal_estimates are where
#                       the correction binds.)
#   heterogeneous_binding: five studies, two precise and three not, with tau2 > 0 (0.069) where the correction still
#                       binds: Hartung-Knapp's se is 0.086 and DerSimonian-Laird's 0.152. In every other binding case
#                       tau2 = 0, where the DerSimonian-Laird se equals the fixed-effect one (0.038 here), so only this
#                       case tells the two floors apart (mzLib #1616 review). Found by a seeded search over rounded
#                       values, kept as written; it uses no random draw, so the cases above are unchanged.
# Each row also carries metafor's test = "knha" se and the DerSimonian-Laird se (test = "z"), so a reader can see which
# one the correction took.

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
  random_eight          = list(yi = rnorm(8, 0.2, 0.3), sei = runif(8, 0.05, 0.3)),
  near_identical_two    = list(yi = c(1.000, 1.001), sei = c(0.30, 0.30)),
  disagreeing_two       = list(yi = c(-0.50, 0.90), sei = c(0.10, 0.12)),
  heterogeneous_binding = list(yi = c(0.34, 0.40, 0.43, 0.24, 0.66), sei = c(0.37, 0.32, 0.39, 0.05, 0.06))
)
inputs <- list(); results <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  for (i in seq_along(cs$yi)) inputs[[length(inputs) + 1]] <- data.frame(case = nm, yi = cs$yi[i], sei = cs$sei[i])
  for (level in c(95, 90)) {
    r <- rma(yi = cs$yi, sei = cs$sei, method = "DL", test = "adhoc", level = level)
    knha <- rma(yi = cs$yi, sei = cs$sei, method = "DL", test = "knha", level = level)
    z <- rma(yi = cs$yi, sei = cs$sei, method = "DL", test = "z", level = level)
    results[[length(results) + 1]] <- data.frame(case = nm, level = level / 100, k = r$k, estimate = as.numeric(r$beta),
      se = r$se, df = r$dfs, tau2 = r$tau2, ci_lb = r$ci.lb, ci_ub = r$ci.ub, pval = r$pval,
      se_knha = knha$se, se_dl = z$se)
  }
}
w(fmt(do.call(rbind, inputs)), "hk_adhoc_inputs.tsv")
w(fmt(do.call(rbind, results)), "hk_adhoc_results.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("metafor", as.character(packageVersion("metafor"))),
  "seed 20261114",
  "script make_hartung_knapp_adhoc_fixtures.R",
  "call rma(yi, sei, method = \"DL\", test = \"adhoc\", level = 95 and 90); se_knha from test = \"knha\", se_dl from test = \"z\""
), file.path(out, "PROVENANCE_hartung_knapp_adhoc.txt"))
