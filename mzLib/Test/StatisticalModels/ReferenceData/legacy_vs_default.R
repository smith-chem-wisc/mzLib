# limma's CURRENT DEFAULT eBayes on the synthetic fixture: the benchmark for mzLib's MarginalLikelihood
# variance prior (MarginalLikelihoodPriorTests). A benchmark, not a parity target: legacy = FALSE uses
# fitFDistUnequalDF1 whenever residual df differ between features, and mzLib does not implement it.
#
# Run ONCE, in CI, and the outputs checked in beside this script. R is never a runtime or test dependency
# of mzLib: the tests read the frozen TSVs.
#
#   Rscript legacy_vs_default.R <out_dir>        (run from the repository root)
#
# Writes limma_default_notrend.tsv, limma_default_trend.tsv and limma_default_prior.tsv. This is the
# synthetic case of the one-off comparison first run on branch ci/statistics-legacy-gap (commit c601fd18),
# which also compared real deposits; those inputs are not part of mzLib, so that part is left out here.

suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
ref <- "mzLib/Test/StatisticalModels/ReferenceData"

read_matrix <- function(path) {
  m <- as.matrix(read.delim(path, check.names = FALSE, na.strings = c("NaN", "NA")))
  storage.mode(m) <- "double"
  m
}

# Keep features whose observed design is full rank with at least one residual df, as the mzLib fixture does:
# limma drops aliased coefficients and keeps the feature, mzLib reports it RankDeficient and leaves it out.
fittable <- function(Y, design) {
  apply(Y, 1, function(y) {
    obs <- !is.na(y)
    sum(obs) > ncol(design) && qr(design[obs, , drop = FALSE])$rank == ncol(design)
  })
}

fmt <- function(x) ifelse(is.na(x), "NaN", sprintf("%.17g", x))
Y <- read_matrix(file.path(ref, "limma_responses.tsv"))
design <- read_matrix(file.path(ref, "limma_design.tsv"))
keep <- fittable(Y, design)
Y <- Y[keep, , drop = FALSE]
fit <- lmFit(Y, design)
j <- match("age_decades", colnames(design))
prior <- list()
for (trend in c(FALSE, TRUE)) {
  new <- eBayes(fit, trend = trend, legacy = FALSE)
  prior[[length(prior) + 1]] <- data.frame(trend = trend, df_prior_default = new$df.prior)
  per <- data.frame(feature = seq_len(nrow(Y)), df_residual = fit$df.residual, amean = fit$Amean,
                    s2_prior_default = fmt(rep_len(new$s2.prior, nrow(Y))), p_default = fmt(new$p.value[, j]))
  write.table(per, file.path(out, sprintf("limma_default_%s.tsv", if (trend) "trend" else "notrend")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
write.table(do.call(rbind, prior), file.path(out, "limma_default_prior.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
# Record these in PROVENANCE_limma_default.txt with the generation time and the CI run.
cat(R.version.string, "\n", "limma", as.character(packageVersion("limma")), "\n")
