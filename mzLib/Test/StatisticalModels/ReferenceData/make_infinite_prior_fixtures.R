# Reference outputs for the variance prior when its degrees of freedom are infinite and there is no trend: the
# feature variances spread no more than sampling explains, so every feature takes the prior variance. limma's
# fitFDist then sets the prior variance to the mean of the (zero-floored) variances, s20 <- mean(x); limma before
# 3.32.0 (2017) returned the limiting value of the finite-df formula, exp(mean(e)), usually a little larger.
#
#   Rscript make_infinite_prior_fixtures.R <out_dir>
#
# The data: 40 features on 6 samples (3 per group), one true residual variance for all, a group effect on the first
# 10. The first seed from 20261110 whose eBayes(legacy = TRUE, trend = FALSE) prior df is infinite is used.

suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}
w <- function(x, name) write.table(fmt(x), file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)

features <- 40; group <- rep(c(0, 1), each = 3)
design <- cbind(intercept = 1, group = group)
make <- function(seed) {
  set.seed(seed)
  effect <- c(rep(1, 10), rep(0, features - 10))
  Y <- matrix(20 + rnorm(features * 6, 0, 0.2), features, 6) + outer(effect, group)
  colnames(Y) <- paste0("s", 1:6)
  Y
}
seed <- 20261110
repeat {
  Y <- make(seed)
  eb <- eBayes(lmFit(Y, design), trend = FALSE, robust = FALSE, legacy = TRUE)
  if (is.infinite(eb$df.prior)) break
  seed <- seed + 1
}
fit <- lmFit(Y, design)
s2 <- fit$sigma^2; d <- fit$df.residual
e <- log(s2) - digamma(d / 2) + log(d / 2)

w(as.data.frame(Y), "limma_infinite_prior_responses.tsv")
w(as.data.frame(design), "limma_infinite_prior_design.tsv")
w(data.frame(s2_prior = rep(eb$s2.prior, features), s2_post = eb$s2.post, df_total = eb$df.total,
             t_group = eb$t[, 2], p_group = eb$p.value[, 2]), "limma_infinite_prior_ebayes.tsv")
w(data.frame(df_prior = eb$df.prior, var_prior = eb$s2.prior, mean_variance = mean(s2),
             pre2017_exp_mean_e = exp(mean(e))), "limma_infinite_prior_prior.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("limma", as.character(packageVersion("limma"))),
  paste("seed", seed, "(the first from 20261110 whose legacy, untrended prior df is infinite)"),
  "script make_infinite_prior_fixtures.R",
  "eBayes(lmFit(Y, design), trend = FALSE, robust = FALSE, legacy = TRUE); fitFDist's infinite-df branch: s20 <- mean(x)",
  sprintf("prior variance %.17g = mean of the variances %.17g; limma before 2017 (exp(mean(e))) would give %.17g",
          eb$s2.prior, mean(s2), exp(mean(e)))
), file.path(out, "PROVENANCE_infinite_prior.txt"))
