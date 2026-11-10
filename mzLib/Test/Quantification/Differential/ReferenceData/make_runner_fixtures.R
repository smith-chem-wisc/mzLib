# Reference values for the STAT1 runner end to end (DifferentialAnalysisReferenceTests): a small label-free dataset with
# known structure, generated here with a fixed seed, and R's answer for every row the runner should fit.
#
#   Rscript make_runner_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# Data: two strata (tissue brain, liver), groups A (reference), B, C with 4 samples each; log2 values = 20 + peptide offset
# + group effect + a per-sample effect of the protein + noise, their sds drawn per protein around 0.15 and 0.1 (log-normal,
# sd 0.4, so the variances differ as real proteins' do and the prior's df is finite), about 8% of values missing. Proteins:
#   P01-P30  3-5 peptides; P01 and P02 also share one peptide, which is unique to neither and so never used
#   P31      no value in group C              P32  no value in the reference group A
#   P33      one peptide seen in group B only P34  two peptides, each sample holds one value of one of them
#   S01-S05  one peptide (the moderated t-test)
#   DECOY, ENTRAP (decoy and entrapment groups) and UNDESC (a group the protein table does not describe): no rows
#
# What R does per stratum, as STAT1 M5b's rulings say (design/stat1/M5B_CONTRACTS_2026-11-04.md), with robust weighting
# and normalization off (the runner's settings in the test):
#   usable peptide: unique to the protein and valued in at least 2 of the stratum's samples (ruling 6);
#   2+ usable peptides: y ~ peptide + group + (1 | sample), REML (lme4 + lmerTest, bobyqa to 1e-12), on the protein's own values with the
#     groups it has; the sample term only when lme4's checks allow it (2 <= levels < values), else lm (ruling 16);
#   1 usable peptide: lm(y ~ group) on its per-sample values (ruling 6);
#   a comparison is answered when both groups have at least 2 samples with a value (rulings 1, 2);
#   moderation per family (peptide model; t-test), squeezeVar(legacy = TRUE); variance df = the joint Satterthwaite
#     F df of the group terms (anova DenDF) or the residual df; contrast df = its own Satterthwaite (or residual) df + d0,
#     capped at the sum of the family's variance df; fewer than 3 fitted in a family: no prior (ruling 5);
#   means: emmeans, equal weight over peptides (ruling 3); Benjamini-Hochberg per (stratum, comparison) over both
#     families' rows (ruling 4); global shift per (stratum, comparison) on every peptide of a non-decoy, non-entrapment
#     group valued in every sample of both sides (GR-27, ruling 9).
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(limma); library(emmeans) })
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
emm_options(lmer.df = "asymptotic")
# lme4's default optimizer can stop short of the optimum (liver P10 by 8e-9 in the criterion); bobyqa to 1e-12, as the other
# mixed-model fixtures, so the reference is the optimum itself.
tight <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                     check.conv.singular = "ignore", check.conv.grad = "ignore")
set.seed(20261110)

groups <- c("A", "B", "C")
values <- NULL
add <- function(tissue, protein, peptide, groupsCell, sample, y)
  values <<- rbind(values, data.frame(tissue = tissue, sample = sample, peptide = peptide, groups = groupsCell, y = y))
for (tissue in c("brain", "liver")) {
  samples <- data.frame(sample = paste0(tissue, "_", rep(groups, each = 4), 1:4), group = rep(groups, each = 4))
  protein <- function(name, npep, effect, keep = function(s, k) TRUE, drop = 0.08) {
    sdSample <- 0.15 * exp(rnorm(1, 0, 0.4)); sdNoise <- 0.1 * exp(rnorm(1, 0, 0.4))
    u <- setNames(rnorm(nrow(samples), 0, sdSample), samples$sample)
    for (k in seq_len(npep)) {
      off <- rnorm(1, 0, 1)
      for (i in seq_len(nrow(samples))) {
        s <- samples$sample[i]; g <- samples$group[i]
        if (!keep(s, k) || runif(1) < drop) next
        add(tissue, name, sprintf("%s_PEP%d", name, k), name, s, 20 + off + effect[[g]] + u[[s]] + rnorm(1, 0, sdNoise))
      }
    }
  }
  for (p in sprintf("P%02d", 1:30)) protein(p, sample(3:5, 1), list(A = 0, B = rnorm(1, 0, 0.8), C = rnorm(1, 0, 0.8)))
  for (s in samples$sample) add(tissue, "P01", "SHARED_PEP", "P01;P02", s, 21 + rnorm(1, 0, 0.1))
  protein("P31", 4, list(A = 0, B = 1, C = 2), keep = function(s, k) !grepl("_C", s))
  protein("P32", 4, list(A = 0, B = 1, C = 2), keep = function(s, k) !grepl("_A", s))
  protein("P33", 3, list(A = 0, B = 0.7, C = -0.5))
  for (s in samples$sample[samples$group == "B"]) add(tissue, "P33", "P33_ONLYB", "P33", s, 22.5 + rnorm(1, 0, 0.1))
  protein("P34", 2, list(A = 0, B = 1.2, C = 0.4), keep = function(s, k) (match(s, samples$sample) %% 2) == (k - 1), drop = 0)
  for (p in sprintf("S%02d", 1:5)) protein(p, 1, list(A = 0, B = rnorm(1, 0, 0.8), C = rnorm(1, 0, 0.8)), drop = 0.05)
  protein("DECOY", 3, list(A = 0, B = 3, C = 3))
  protein("ENTRAP", 3, list(A = 0, B = 3, C = 3))
  protein("UNDESC", 3, list(A = 0, B = 0.4, C = 0.4))
}
values$group <- substr(sub("^[a-z]+_", "", values$sample), 1, 1)

contrasts <- list(c1 = c("B", "A"), c2 = c("C", "A"), c3 = c("C", "B"))
rows <- NULL; shifts <- NULL
for (tissue in c("brain", "liver")) {
  v <- values[values$tissue == tissue, ]
  stratum <- paste0("tissue=", tissue)
  # Global shift (GR-27): eligible peptides, every sample of both sides valued (no sample here is a lost run).
  eligible <- v[!(v$groups %in% c("DECOY", "ENTRAP")), ]
  for (cid in names(contrasts)) {
    cc <- contrasts[[cid]]
    ns <- unique(v$sample[v$group == cc[1]]); ds <- unique(v$sample[v$group == cc[2]])
    diffs <- c()
    for (pep in unique(eligible$peptide)) {
      x <- eligible[eligible$peptide == pep, ]
      if (!all(ns %in% x$sample) || !all(ds %in% x$sample)) next
      diffs <- c(diffs, mean(x$y[x$sample %in% ns]) - mean(x$y[x$sample %in% ds]))
    }
    shifts <- rbind(shifts, data.frame(stratum = stratum, contrast = cid, value = fmt(if (length(diffs)) median(diffs) else NA),
                                       peptides = length(diffs)))
  }
  fits <- list()
  for (p in c(sprintf("P%02d", 1:34), sprintf("S%02d", 1:5))) {
    x <- v[v$groups == p, ]
    counts <- table(x$peptide)
    usable <- names(counts)[sapply(names(counts), function(pp) length(unique(x$sample[x$peptide == pp])) >= 2)]
    x <- x[x$peptide %in% usable, ]
    if (length(usable) == 0) next
    present <- groups[groups %in% x$group]
    x$group <- factor(x$group, levels = present)
    single <- length(usable) == 1
    if (single) {
      m <- lm(y ~ group, data = x); kind <- "lm"
    } else {
      x$peptide <- factor(x$peptide, levels = sort(usable, method = "radix"))
      k <- length(unique(x$sample))
      kind <- if (k >= 2 && k < nrow(x)) "lmer" else "lm"
      m <- if (kind == "lmer") lmer(y ~ peptide + group + (1 | sample), data = x, REML = TRUE, control = tight) else lm(y ~ peptide + group, data = x)
    }
    b <- if (kind == "lmer") fixef(m) else coef(m)
    vdf <- if (kind == "lmer") anova(m)["group", "DenDF"] else df.residual(m)
    em <- summary(emmeans(m, ~ group))
    mean_of <- setNames(em$emmean, as.character(em$group))
    out <- list()
    for (cid in names(contrasts)) {
      cc <- contrasts[[cid]]
      n1 <- length(unique(x$sample[x$group == cc[1]])); n2 <- length(unique(x$sample[x$group == cc[2]]))
      if (n1 < 2 || n2 < 2) next
      L <- setNames(rep(0, length(b)), names(b))
      for (j in 1:2) { nm <- paste0("group", cc[j]); if (nm %in% names(L)) L[nm] <- L[nm] + c(1, -1)[j] }
      if (kind == "lmer") {
        ct <- contest(m, L, joint = FALSE); est <- ct$Estimate; se <- ct$`Std. Error`; df <- ct$df
      } else {
        est <- sum(L * b); se <- sqrt(drop(t(L) %*% vcov(m) %*% L)); df <- df.residual(m)
      }
      out[[cid]] <- data.frame(contrast = cid, est = est, se = se, cdf = df, mean_num = mean_of[[cc[1]]], mean_den = mean_of[[cc[2]]])
    }
    if (length(out) == 0) next
    fits[[p]] <- list(family = if (single) "t" else "mixed", s2 = sigma(m)^2, vdf = vdf, rows = do.call(rbind, out),
                      kind = kind)
  }
  for (family in c("mixed", "t")) {
    members <- names(fits)[sapply(fits, function(f) f$family == family)]
    s2 <- sapply(fits[members], `[[`, "s2"); vdf <- sapply(fits[members], `[[`, "vdf")
    if (length(members) >= 3) {
      sq <- squeezeVar(s2, vdf, legacy = TRUE); post <- sq$var.post; d0 <- sq$df.prior
    } else {
      post <- s2; d0 <- 0
    }
    cap <- sum(vdf)
    for (i in seq_along(members)) {
      f <- fits[[members[i]]]; r <- f$rows
      se <- sqrt(r$se^2 / f$s2 * post[i]); df <- pmin(r$cdf + d0, cap)
      rows <- rbind(rows, data.frame(stratum = stratum, protein = members[i], contrast = r$contrast,
        status = if (family == "t") "fitted_single_peptide" else "fitted", model = f$kind,
        log2 = r$est, se = se, df = df, p = 2 * pt(-abs(r$est / se), df), mean_num = r$mean_num, mean_den = r$mean_den))
    }
  }
}
rows$p_adj <- ave(rows$p, paste(rows$stratum, rows$contrast), FUN = function(p) p.adjust(p, "BH"))
for (col in c("log2", "se", "df", "p", "p_adj", "mean_num", "mean_den")) rows[[col]] <- fmt(rows[[col]])
write.table(values[, c("tissue", "sample", "group", "peptide", "groups", "y")], file.path(here, "runner_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
values$y <- fmt(values$y)
write.table(values[, c("tissue", "sample", "group", "peptide", "groups", "y")], file.path(here, "runner_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rows, file.path(here, "runner_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(shifts, file.path(here, "runner_shifts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4")), "lmerTest", as.character(packageVersion("lmerTest")),
        "limma", as.character(packageVersion("limma")), "emmeans", as.character(packageVersion("emmeans"))),
  "script make_runner_fixtures.R (seed 20261110)",
  paste("rows", nrow(rows), "fitted; lm fits", paste(unique(rows$protein[rows$model == "lm" & rows$status == "fitted"]), collapse = ","))
), file.path(here, "PROVENANCE_runner.txt"))
cat("written", nrow(rows), "rows\n")
