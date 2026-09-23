# Reference outputs for mzLib's Statistics project, from the published implementations.
#
# Run ONCE, in CI (a manual-dispatch job on a throwaway branch), and the outputs checked in beside
# this script. R is never a runtime or test dependency of mzLib: the tests read the frozen TSVs.
#
#   Rscript make_reference_fixtures.R <out_dir>
#
# What is compared, and why these settings:
#   * limma lmFit + eBayes, trend = FALSE and trend = TRUE, robust = FALSE.
#   * legacy = TRUE. From limma 3.61 eBayes switches to a different prior estimator
#     (fitFDistUnequalDF1) whenever residual df differ between features, which missing values cause.
#     mzLib implements the method-of-moments estimator of Smyth (2004), which is the legacy one.
#   * metafor rma(method = "DL") for DerSimonian-Laird pooling.
#   * p.adjust(method = "BH").
# Every feature's observed design is kept full-rank, because limma drops aliased coefficients and
# keeps the feature, while mzLib reports it RankDeficient and leaves it out of the prior.

suppressPackageStartupMessages({ library(limma); library(metafor) })
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
                                   row.names = FALSE, na = "NaN")
fmt <- function(df) {                                  # full precision; NA written as NaN for .NET
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20260923)

# ---- a design with a continuous covariate and a categorical one ------------------------------------
G <- 400; n <- 12
age <- seq(20, 80, length.out = n)
design <- cbind(Intercept = 1, age_decades = (age - 50) / 10, sex = rep(0:1, n / 2))

A <- rnorm(G, 22, 2)                                   # average log2 intensity
s2 <- 0.05 * 20 / rchisq(G, 20) * exp(-(A - 22) / 3)   # variance falls with intensity: a real trend
beta <- ifelse(runif(G) < 0.1, rnorm(G, 0, 0.5), 0)
Y <- matrix(A, G, n) + outer(beta, design[, "age_decades"]) + matrix(rnorm(G * n), G, n) * sqrt(s2)
miss <- matrix(runif(G * n) < 0.08, G, n)
for (g in seq_len(G)) {                                # keep each feature's observed design full-rank
  obs <- !miss[g, ]
  if (sum(obs) <= ncol(design) || qr(design[obs, , drop = FALSE])$rank < ncol(design)) miss[g, ] <- FALSE
}
Y[miss] <- NA

w(fmt(as.data.frame(Y)), "limma_responses.tsv")
w(fmt(as.data.frame(design)), "limma_design.tsv")

fit <- lmFit(Y, design)
w(fmt(data.frame(
  coef_intercept = fit$coefficients[, 1], coef_age = fit$coefficients[, 2], coef_sex = fit$coefficients[, 3],
  unscaled_age = fit$stdev.unscaled[, 2], sigma = fit$sigma, df_residual = fit$df.residual, amean = fit$Amean
)), "limma_fit.tsv")

for (trend in c(FALSE, TRUE)) {
  eb <- eBayes(fit, trend = trend, robust = FALSE, legacy = TRUE)
  tag <- if (trend) "trend" else "notrend"
  w(fmt(data.frame(
    s2_prior = if (length(eb$s2.prior) == 1) rep(eb$s2.prior, G) else eb$s2.prior,
    s2_post = eb$s2.post, df_total = eb$df.total,
    t_age = eb$t[, 2], p_age = eb$p.value[, 2], bh_age = p.adjust(eb$p.value[, 2], "BH")
  )), sprintf("limma_ebayes_%s.tsv", tag))
  w(fmt(data.frame(df_prior = eb$df.prior)), sprintf("limma_ebayes_%s_prior.tsv", tag))
}

# ---- DerSimonian-Laird, several heterogeneity regimes ------------------------------------------------
meta <- list()
cases <- list(
  homogeneous   = list(yi = c(0.30, 0.28, 0.35, 0.31, 0.29), sei = c(0.10, 0.12, 0.09, 0.11, 0.10)),
  heterogeneous = list(yi = c(0.10, 0.60, -0.20, 0.45, 0.05, 0.80), sei = c(0.10, 0.15, 0.12, 0.20, 0.08, 0.25)),
  two_studies   = list(yi = c(-0.40, 0.10), sei = c(0.20, 0.30)),
  random_eight  = list(yi = rnorm(8, 0.2, 0.3), sei = runif(8, 0.05, 0.3))
)
rows <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  r <- rma(yi = cs$yi, sei = cs$sei, method = "DL")
  for (i in seq_along(cs$yi)) rows[[length(rows) + 1]] <- data.frame(case = nm, yi = cs$yi[i], sei = cs$sei[i])
  meta[[nm]] <- data.frame(case = nm, estimate = as.numeric(r$beta), se = r$se, tau2 = r$tau2,
                           q = r$QE, i2 = r$I2 / 100, ci_lb = r$ci.lb, ci_ub = r$ci.ub, pval = r$pval)
}
w(fmt(do.call(rbind, rows)), "dl_inputs.tsv")
w(fmt(do.call(rbind, meta)), "dl_results.tsv")

# ---- provenance ---------------------------------------------------------------------------------------
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("limma", as.character(packageVersion("limma"))),
  paste("metafor", as.character(packageVersion("metafor"))),
  "seed 20260923",
  "script make_reference_fixtures.R"
), file.path(out, "PROVENANCE.txt"))
