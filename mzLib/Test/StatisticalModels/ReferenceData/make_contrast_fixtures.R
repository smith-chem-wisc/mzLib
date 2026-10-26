# Reference values for moderated CONTRASTS with confidence intervals (STAT1 milestone M2):
# EmpiricalBayes.ModerateContrasts against limma.
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   Rscript make_contrast_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# Two data sets share one design (intercept, group B, group C, age per decade):
#   complete  no missing values. Reference = limma itself: lmFit, contrasts.fit, eBayes(legacy = TRUE),
#             topTable(confint = 0.95).
#   na        about 10% missing, omitted per feature. limma's contrasts.fit APPROXIMATES a contrast's
#             unscaled SD when missing values make a feature's coefficients correlate differently from the
#             full design, so it is not the reference here. Reference = the exact per-feature value:
#             lm on the feature's observed rows, unscaled SD = sqrt(c' (Xo'Xo)^-1 c), prior from
#             limma::squeezeVar(legacy = TRUE) (what eBayes(legacy = TRUE) uses), df.total capped at the
#             pooled residual df, CI = estimate +/- qt((1 + 0.95) / 2, df.total) * s.post * unscaled SD.
# Each with and without the intensity trend (eBayes(trend = TRUE) / squeezeVar(covariate = Amean)).

suppressPackageStartupMessages(library(limma))
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", ifelse(is.infinite(v), ifelse(v > 0, "Inf", "-Inf"), sprintf("%.17g", v)))
put <- function(df, name) write.table(df, file.path(here, name), sep = "\t", quote = FALSE, row.names = FALSE)

seed <- 20261008
set.seed(seed)
G <- 300
group <- factor(rep(c("A", "B", "C"), each = 4), levels = c("A", "B", "C"))
age <- c(23, 35, 47, 61, 29, 41, 55, 68, 26, 38, 52, 70)
age_decade <- (age - 50) / 10
X <- model.matrix(~ group + age_decade)
colnames(X) <- c("intercept", "groupB", "groupC", "age_decade")
n <- nrow(X); p <- ncol(X)

cm <- cbind(C_vs_B = c(0, -1, 1, 0), age = c(0, 0, 0, 1), BC_vs_A = c(0, 0.5, 0.5, 0))
rownames(cm) <- colnames(X)

base <- rnorm(G, 20, 2)
beta <- cbind(base, rnorm(G, 0, 0.6), rnorm(G, 0, 0.6), rnorm(G, 0, 0.3))
sdev <- sqrt(0.2 * 4 / rchisq(G, 4)) * exp(-(base - 20) / 8)   # variances spread, trending with intensity
Y <- beta %*% t(X) + matrix(rnorm(G * n), G, n) * sdev
rownames(Y) <- sprintf("f%03d", seq_len(G))

Yna <- Y
for (g in seq_len(G)) {
  repeat {
    miss <- runif(n) < 0.10
    obs <- which(!miss)
    if (length(obs) >= p + 2 && qr(X[obs, , drop = FALSE])$rank == p) break
  }
  Yna[g, miss] <- NA
}

put(data.frame(X, check.names = FALSE), "limma_contrast_design.tsv")
put(data.frame(contrast = colnames(cm), t(cm), check.names = FALSE), "limma_contrast_weights.tsv")
put(data.frame(feature = rownames(Y), apply(Y, 2, fmt), check.names = FALSE), "limma_contrast_responses_complete.tsv")
put(data.frame(feature = rownames(Yna), apply(Yna, 2, fmt), check.names = FALSE), "limma_contrast_responses_na.tsv")

rows_out <- function(est, su, s2post, s2prior, dfp, dft, label) {
  out <- NULL
  for (k in seq_len(ncol(cm))) {
    se <- sqrt(s2post) * su[, k]
    t <- est[, k] / se
    pv <- 2 * pt(-abs(t), dft)
    q <- qt(0.975, dft)   # as topTable: t on df.total, which the pooled-df cap keeps finite even when df.prior is Inf
    out <- rbind(out, data.frame(
      contrast = colnames(cm)[k], feature = rownames(Y),
      estimate = fmt(est[, k]), se = fmt(se), t = fmt(t), df_total = fmt(dft), p = fmt(pv),
      adj_p = fmt(p.adjust(pv, "BH")), ci_low = fmt(est[, k] - q * se), ci_high = fmt(est[, k] + q * se),
      s2_post = fmt(s2post), s2_prior = fmt(s2prior), df_prior = fmt(rep(dfp, G)), check.names = FALSE))
  }
  put(out, sprintf("limma_contrast_%s.tsv", label))
  cat(label, ": df.prior ", dfp, "\n", sep = "")
}

# complete: limma itself
for (tr in c(FALSE, TRUE)) {
  fit <- lmFit(Y, X)
  fit2 <- contrasts.fit(fit, cm)
  eb <- eBayes(fit2, trend = tr, legacy = TRUE)
  label <- sprintf("complete_%s", if (tr) "trend" else "none")
  rows_out(fit2$coefficients, fit2$stdev.unscaled, eb$s2.post, rep_len(eb$s2.prior, G), eb$df.prior, eb$df.total, label)
  # topTable's own CI must agree with the formula written above
  for (k in seq_len(ncol(cm))) {
    tt <- topTable(eb, coef = k, number = Inf, sort.by = "none", confint = 0.95)
    se <- sqrt(eb$s2.post) * fit2$stdev.unscaled[, k]
    q <- qt(0.975, eb$df.total)
    stopifnot(max(abs(tt$CI.L - (fit2$coefficients[, k] - q * se))) < 1e-12)
    stopifnot(max(abs(tt$P.Value - eb$p.value[, k])) < 1e-15)
  }
}

# na: exact per-feature values, prior as eBayes(legacy = TRUE) fits it
est <- su <- matrix(NA_real_, G, ncol(cm)); s2 <- dfr <- amean <- numeric(G)
for (g in seq_len(G)) {
  obs <- which(!is.na(Yna[g, ]))
  Xo <- X[obs, , drop = FALSE]; yo <- Yna[g, obs]
  lf <- lm.fit(Xo, yo)
  dfr[g] <- length(obs) - p
  s2[g] <- sum(lf$residuals^2) / dfr[g]
  V <- chol2inv(qr.R(qr(Xo)))
  est[g, ] <- drop(t(cm) %*% lf$coefficients)
  su[g, ] <- sqrt(diag(t(cm) %*% V %*% cm))
  amean[g] <- mean(yo)
}
rownames(est) <- rownames(su) <- rownames(Y)
for (tr in c(FALSE, TRUE)) {
  sv <- squeezeVar(s2, dfr, covariate = if (tr) amean else NULL, legacy = TRUE)
  dft <- pmin(dfr + sv$df.prior, sum(dfr))
  label <- sprintf("na_%s", if (tr) "trend" else "none")
  rows_out(est, su, sv$var.post, rep_len(sv$var.prior, G), sv$df.prior, dft, label)
}

writeLines(c(
  sprintf("generated_utc %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R %s", R.version.string),
  sprintf("limma %s", as.character(packageVersion("limma"))),
  sprintf("seed %d", seed),
  "script make_contrast_fixtures.R",
  "files limma_contrast_{design,weights,responses_complete,responses_na}.tsv and limma_contrast_{complete,na}_{none,trend}.tsv",
  "complete = limma lmFit + contrasts.fit + eBayes(legacy = TRUE) + topTable(confint = 0.95); na = exact per-feature unscaled SD with squeezeVar(legacy = TRUE)"
), file.path(here, "PROVENANCE_contrast.txt"))
