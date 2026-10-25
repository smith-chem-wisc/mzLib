# Reference outputs for mzLib's SpearmanCorrelation.PartialCorrelate, from ppcor::pcor.test(method = "spearman").
#
# Run once, and the outputs checked in beside this script. R is never a runtime or test dependency of mzLib:
# the tests read the frozen TSVs.
#
#   Rscript make_partial_spearman_fixtures.R <out_dir>
#
# Cases: one to three covariates, small and larger n, with and without ties, a strong partial correlation and one
# that a covariate explains away. ppcor ranks every variable with ties averaged, takes the partial correlation from
# the inverse of the rank correlation matrix, and tests it with t on n - 2 - k df, which is what mzLib does.

suppressPackageStartupMessages(library(ppcor))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20261012)
inputs <- list(); results <- list()
add_case <- function(name, x, y, z) {
  z <- as.matrix(z)
  pt <- pcor.test(x, y, z, method = "spearman")
  zz <- matrix(NA_real_, length(x), 3)
  zz[, seq_len(ncol(z))] <- z
  inputs[[length(inputs) + 1]] <<- data.frame(case = name, x = x, y = y, z1 = zz[, 1], z2 = zz[, 2], z3 = zz[, 3])
  results[[length(results) + 1]] <<- data.frame(case = name, k = ncol(z), n = pt$n, rho = pt$estimate, statistic = pt$statistic,
                                                p = pt$p.value)
}
for (n in c(6, 10, 30)) for (k in 1:3) {
  if (n - 2 - k < 2) next
  z <- matrix(rnorm(n * k), n, k)
  x <- drop(z %*% rnorm(k)) + rnorm(n)
  y <- 0.5 * x + drop(z %*% rnorm(k)) + rnorm(n)
  add_case(sprintf("n%d_k%d", n, k), x, y, z)
}
# Ties in every variable.
z <- matrix(sample(1:4, 40, replace = TRUE), 20, 2)
add_case("ties_n20_k2", sample(1:6, 20, replace = TRUE), sample(1:5, 20, replace = TRUE), z)
# A covariate that drives both: a clear raw correlation, little partial correlation.
z <- rnorm(25)
add_case("explained_n25_k1", z + rnorm(25, 0, 0.3), z + rnorm(25, 0, 0.3), z)
# A strong partial correlation, given a run-level covariate.
z <- rnorm(40); e <- rnorm(40)
add_case("strong_n40_k1", 2 * z + e, -z + e + rnorm(40, 0, 0.2), z)
w(fmt(do.call(rbind, inputs)), "partial_spearman_inputs.tsv")
w(fmt(do.call(rbind, results)), "partial_spearman_results.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("ppcor", as.character(packageVersion("ppcor"))),
  "seed 20261012",
  "script make_partial_spearman_fixtures.R",
  "run locally (Rscript make_partial_spearman_fixtures.R <this folder>), not in CI"
), file.path(out, "PROVENANCE_partial_spearman.txt"))
