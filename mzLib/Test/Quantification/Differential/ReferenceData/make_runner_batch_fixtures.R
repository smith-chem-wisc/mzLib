# Reference values for the STAT1 runner's means when a protein was never measured in a batch (STAT1 M5b ruling 18;
# DifferentialAnalysisReferenceTests): a small label-free dataset with batches, generated here with a fixed seed, and R's
# answer for every row the runner should fit.
#
#   Rscript make_runner_batch_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# Data: one stratum, groups A (reference) and B, batches b1 (reference), b2, b3, two samples per group in each batch (12
# samples, so every batch holds the same number and emmeans's equal weights over batches are the stratum's proportions).
# log2 values = 20 + peptide offset + group effect + batch effect + a per-sample effect of the protein + noise, about 5% of
# values missing. Proteins:
#   B01-B08  3-4 peptides, every batch
#   B09      never measured in b3          B10  never measured in b1 (the reference batch)
#   B11      measured in b2 only           B12  measured in b3 only in group B (b3 is still a batch it was measured in)
#
# What R does for each protein, as the runner should (rulings 1, 3, 6, 16, 18): its usable peptides (valued in at least 2
# samples); y ~ peptide + group + batch + (1 | sample), REML (lme4 + lmerTest, bobyqa to 1e-12), on its own values with the
# batches it has (batch left out when it has one; the sample term only when lme4's checks allow it, else lm); moderation over
# the proteins with squeezeVar(legacy = TRUE), variance df = the group term's Satterthwaite DenDF; means by emmeans, equal
# weight over the peptides and over the batches the protein has.
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(limma); library(emmeans) })
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
emm_options(lmer.df = "asymptotic")
tight <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                     check.conv.singular = "ignore", check.conv.grad = "ignore")
set.seed(20261111)

batches <- c("b1", "b2", "b3")
samples <- expand.grid(rep = 1:2, group = c("A", "B"), batch = batches, stringsAsFactors = FALSE)
samples$sample <- sprintf("%s_%s%d", samples$batch, samples$group, samples$rep)
batchEffect <- c(b1 = 0, b2 = 0.4, b3 = -0.6)
values <- NULL
protein <- function(name, npep, effect, keep = function(s) TRUE, drop = 0.05) {
  sdSample <- 0.15 * exp(rnorm(1, 0, 0.4)); sdNoise <- 0.1 * exp(rnorm(1, 0, 0.4))
  u <- setNames(rnorm(nrow(samples), 0, sdSample), samples$sample)
  for (k in seq_len(npep)) {
    off <- rnorm(1, 0, 1)
    for (i in seq_len(nrow(samples))) {
      s <- samples[i, ]
      if (!keep(s) || runif(1) < drop) next
      y <- 20 + off + (if (s$group == "B") effect else 0) + batchEffect[[s$batch]] + u[[s$sample]] + rnorm(1, 0, sdNoise)
      values <<- rbind(values, data.frame(sample = s$sample, group = s$group, batch = s$batch,
                                          peptide = sprintf("%s_PEP%d", name, k), groups = name, y = y))
    }
  }
}
for (p in sprintf("B%02d", 1:8)) protein(p, sample(3:4, 1), rnorm(1, 0, 0.8))
protein("B09", 3, 1, keep = function(s) s$batch != "b3")
protein("B10", 3, -0.7, keep = function(s) s$batch != "b1")
protein("B11", 3, 0.5, keep = function(s) s$batch == "b2", drop = 0)
protein("B12", 3, 0.9, keep = function(s) s$batch != "b3" || s$group == "B")

rows <- NULL; fits <- list()
for (p in sprintf("B%02d", 1:12)) {
  x <- values[values$groups == p, ]
  usable <- unique(x$peptide)[sapply(unique(x$peptide), function(pp) length(unique(x$sample[x$peptide == pp])) >= 2)]
  x <- x[x$peptide %in% usable, ]
  x$peptide <- factor(x$peptide, levels = sort(usable, method = "radix"))
  x$group <- factor(x$group, levels = c("A", "B"))
  present <- batches[batches %in% x$batch]
  x$batch <- factor(x$batch, levels = present)
  fixed <- if (length(present) > 1) "y ~ peptide + group + batch" else "y ~ peptide + group"
  k <- length(unique(x$sample))
  kind <- if (k >= 2 && k < nrow(x)) "lmer" else "lm"
  m <- if (kind == "lmer") lmer(as.formula(paste(fixed, "+ (1 | sample)")), data = x, REML = TRUE, control = tight) else
    lm(as.formula(fixed), data = x)
  em <- summary(emmeans(m, ~ group))
  mean_of <- setNames(em$emmean, as.character(em$group))
  b <- if (kind == "lmer") fixef(m) else coef(m)
  L <- setNames(rep(0, length(b)), names(b)); L["groupB"] <- 1
  if (kind == "lmer") {
    ct <- contest(m, L, joint = FALSE); est <- ct$Estimate; se <- ct$`Std. Error`; cdf <- ct$df
    vdf <- anova(m)["group", "DenDF"]
  } else {
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% vcov(m) %*% L)); cdf <- df.residual(m); vdf <- cdf
  }
  fits[[p]] <- list(kind = kind, s2 = sigma(m)^2, vdf = vdf, est = est, se = se, cdf = cdf,
                    mean_num = mean_of[["B"]], mean_den = mean_of[["A"]],
                    seen = paste(present, collapse = ","))
}
for (family in c("lmer", "lm")) {
  members <- names(fits)[sapply(fits, function(f) f$kind == family)]
  if (length(members) == 0) next
  s2 <- sapply(fits[members], `[[`, "s2"); vdf <- sapply(fits[members], `[[`, "vdf")
  if (length(members) >= 3) {
    sq <- squeezeVar(s2, vdf, legacy = TRUE); post <- sq$var.post; d0 <- sq$df.prior
  } else {
    post <- s2; d0 <- 0
  }
  cap <- sum(vdf)
  for (i in seq_along(members)) {
    f <- fits[[members[i]]]
    se <- sqrt(f$se^2 / f$s2 * post[i]); df <- min(f$cdf + d0, cap)
    rows <- rbind(rows, data.frame(protein = members[i], contrast = "c1", model = f$kind, batches_seen = f$seen,
      log2 = f$est, se = se, df = df, p = 2 * pt(-abs(f$est / se), df), mean_num = f$mean_num, mean_den = f$mean_den))
  }
}
rows <- rows[order(rows$protein), ]
for (col in c("log2", "se", "df", "p", "mean_num", "mean_den")) rows[[col]] <- fmt(rows[[col]])
values$y <- fmt(values$y)
write.table(values[, c("sample", "group", "batch", "peptide", "groups", "y")], file.path(here, "runner_batch_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rows, file.path(here, "runner_batch_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4")), "lmerTest", as.character(packageVersion("lmerTest")),
        "limma", as.character(packageVersion("limma")), "emmeans", as.character(packageVersion("emmeans"))),
  "script make_runner_batch_fixtures.R (seed 20261111)",
  paste("rows", nrow(rows), "fitted; models", paste(names(table(rows$model)), table(rows$model), collapse = " "))
), file.path(here, "PROVENANCE_runner_batch.txt"))
cat("written", nrow(rows), "rows\n")
print(rows[, c("protein", "model", "batches_seen", "log2", "mean_num", "mean_den")])
