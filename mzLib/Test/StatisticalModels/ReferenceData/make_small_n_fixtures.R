# Reference outputs for the intensity trend on FEW features, where limma's default spline size depends
# on how many features there are (fitFDist: 1 + (n >= 3) + (n >= 6) + (n >= 30) basis functions,
# capped at the number of distinct covariate values; below 2 it fits no trend).
#
# Run ONCE, in CI on a throwaway branch, and the output checked in beside this script, as for
# make_reference_fixtures.R. It reads that script's frozen inputs, so it needs no seed:
#
#   Rscript make_small_n_fixtures.R <this_dir> <out_dir>
#
# Each case is the first n features of limma_responses.tsv, fitted with lmFit and moderated with
# eBayes(trend = TRUE, robust = FALSE, legacy = TRUE).

suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
inp <- if (length(args) >= 1) args[[1]] else "."
out <- if (length(args) >= 2) args[[2]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
                                   row.names = FALSE, na = "NaN")
fmt <- function(df) {                                  # full precision; NA written as NaN for .NET
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

Y <- as.matrix(read.delim(file.path(inp, "limma_responses.tsv"), na.strings = "NaN"))
design <- as.matrix(read.delim(file.path(inp, "limma_design.tsv")))

features <- list(); priors <- list()
for (n in c(2L, 3L, 5L, 6L, 20L, 29L, 30L)) {
  fit <- lmFit(Y[seq_len(n), , drop = FALSE], design)
  eb <- eBayes(fit, trend = TRUE, robust = FALSE, legacy = TRUE)
  features[[length(features) + 1]] <- data.frame(
    n = n, feature = seq_len(n) - 1L,
    s2_prior = if (length(eb$s2.prior) == 1) rep(eb$s2.prior, n) else eb$s2.prior,
    s2_post = eb$s2.post, df_total = eb$df.total, t_age = eb$t[, 2], p_age = eb$p.value[, 2])
  priors[[length(priors) + 1]] <- data.frame(n = n, df_prior = eb$df.prior)
}
w(fmt(do.call(rbind, features)), "limma_ebayes_trend_small_n.tsv")
w(fmt(do.call(rbind, priors)), "limma_ebayes_trend_small_n_prior.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("limma", as.character(packageVersion("limma"))),
  "script make_small_n_fixtures.R (inputs: limma_responses.tsv, limma_design.tsv)"
), file.path(out, "PROVENANCE_small_n.txt"))
