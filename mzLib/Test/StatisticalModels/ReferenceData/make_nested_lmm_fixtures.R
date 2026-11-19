# Reference values for NestedMixedModel (STAT1 milestone M3a): a per-protein linear mixed model on
# peptide-level log2 intensities, with Satterthwaite degrees of freedom and MSstatsTMT-style variance
# moderation (QuantProject design/STATS-FRAMEWORK.md GR-6, GR-7, GR-21).
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   Rscript make_nested_lmm_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# Per protein:  y ~ common fixed effects + peptide (treatment-coded) + random intercepts, REML (lme4::lmer),
#   random part  (1 | individual) + (1 | individual:sample)  in scenario "nested"
#                (1 | sample)                                in scenario "sample_only"
#   df           lmerTest::contest (Satterthwaite) per contrast; the JOINT F-test of the condition terms gives
#                each protein's variance df (DenDF), as MSstatsTMT uses the Group F-test's DenDF.
#   moderation   limma::squeezeVar(sigma^2, DenDF, legacy = TRUE); SE = sqrt(c' unsc c * s2.post);
#                df = Satterthwaite df of the contrast + df.prior, capped at sum(DenDF); t-based 95% CI.
#                (MSstatsTMT 2.20.0 MSstatsTestSingleProteinTMT, extended to several condition terms.)

suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(limma) })
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", ifelse(is.infinite(v), ifelse(v > 0, "Inf", "-Inf"), sprintf("%.17g", v)))
put <- function(df, name) write.table(df, file.path(here, name), sep = "\t", quote = FALSE, row.names = FALSE)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)

seed <- 20261008
set.seed(seed)

scenario <- function(name) {
  if (name == "nested") {
    samples <- data.frame(sample = sprintf("s%02d", 1:12), individual = sprintf("m%d", rep(1:6, each = 2)),
                          condition = rep(c("A", "B"), 6), genotype = rep(c("wt", "ko"), each = 6))
    common <- cbind(intercept = 1, conditionB = as.numeric(samples$condition == "B"),
                    genotypeko = as.numeric(samples$genotype == "ko"))
    conditionRows <- list(c(0, 1, 0), c(0, 0, 1))                 # the joint test: both condition terms
    contrasts <- rbind(conditionB = c(0, 1, 0), genotypeko = c(0, 0, 1), sum = c(0, 1, 1))
    random <- "(1 | individual) + (1 | individual:sample)"
  } else {
    samples <- data.frame(sample = sprintf("s%02d", 1:12), individual = NA,
                          group = rep(c("A", "B", "C"), each = 4))
    common <- cbind(intercept = 1, groupB = as.numeric(samples$group == "B"), groupC = as.numeric(samples$group == "C"))
    conditionRows <- list(c(0, 1, 0), c(0, 0, 1))
    contrasts <- rbind(groupB = c(0, 1, 0), groupC = c(0, 0, 1), C_vs_B = c(0, -1, 1))
    random <- "(1 | sample)"
  }
  colnames(contrasts) <- colnames(common)
  list(name = name, samples = samples, common = common, conditionRows = conditionRows, contrasts = contrasts, random = random)
}

simulate <- function(sc, P) {
  n <- nrow(sc$samples)
  rows <- NULL
  for (p in seq_len(P)) {
    k <- sample(2:6, 1)
    if (p == 1) k <- 1                                              # one peptide: not fittable
    beta <- c(rnorm(1, 22, 1.5), rnorm(ncol(sc$common) - 1, 0, 0.5))
    pep <- c(0, rnorm(k - 1, 0, 1))
    tauInd <- c(0, 0.15, 0.4)[(p %% 3) + 1]                         # some proteins with no individual variance
    tauSam <- c(0.3, 0, 0.2)[(p %% 3) + 1]
    sigma <- sqrt(0.08 * 4 / rchisq(1, 4))
    ind <- if (all(is.na(sc$samples$individual))) rep(0, n) else rnorm(6, 0, tauInd)[as.integer(factor(sc$samples$individual))]
    sam <- rnorm(n, 0, tauSam)
    for (j in seq_len(k)) {
      y <- drop(sc$common %*% beta) + pep[j] + ind + sam + rnorm(n, 0, sigma)
      y[runif(n) < 0.08] <- NA
      rows <- rbind(rows, data.frame(protein = sprintf("p%03d", p), peptide = sprintf("pep%d", j),
                                     sample = sc$samples$sample, y = y))
    }
  }
  rows
}

run <- function(sc, P) {
  long <- simulate(sc, P)
  put(data.frame(sc$samples, sc$common, check.names = FALSE), sprintf("lmm_%s_samples.tsv", sc$name))
  put(data.frame(protein = long$protein, peptide = long$peptide, sample = long$sample, y = fmt(long$y)),
      sprintf("lmm_%s_values.tsv", sc$name))
  put(data.frame(contrast = rownames(sc$contrasts), sc$contrasts, check.names = FALSE), sprintf("lmm_%s_contrasts.tsv", sc$name))

  fits <- list()
  for (prot in unique(long$protein)) {
    d <- merge(long[long$protein == prot & !is.na(long$y), ], sc$samples, by = "sample")
    d$peptide <- factor(d$peptide, levels = sort(unique(d$peptide), method = "radix"))
    for (cn in colnames(sc$common)[-1]) d[[cn]] <- sc$common[match(d$sample, sc$samples$sample), cn]
    terms <- paste(colnames(sc$common)[-1], collapse = " + ")
    f <- as.formula(sprintf("y ~ %s%s + %s", terms, if (nlevels(d$peptide) > 1) " + peptide" else "", sc$random))
    fit <- tryCatch(lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl), error = function(e) NULL)
    if (is.null(fit) || nlevels(d$peptide) < 2) { fits[[prot]] <- NULL; next }
    p <- ncol(sc$common)
    vc <- as.data.frame(VarCorr(fit))
    getv <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v else NA }
    sig2 <- sigma(fit)^2
    unsc <- as.matrix(vcov(fit))[1:p, 1:p] / sig2
    L <- matrix(0, length(sc$conditionRows), length(fixef(fit)))
    for (r in seq_along(sc$conditionRows)) L[r, 1:p] <- sc$conditionRows[[r]]
    denDf <- contest(fit, L, joint = TRUE)$DenDF
    cdf <- sapply(seq_len(nrow(sc$contrasts)), function(i) {
      l <- numeric(length(fixef(fit))); l[1:p] <- sc$contrasts[i, ]; contest(fit, l, joint = FALSE)$df })
    fits[[prot]] <- list(beta = fixef(fit)[1:p], unsc = unsc, sig2 = sig2, ind = getv("individual"),
                         sam = if (sc$name == "nested") getv("individual:sample") else getv("sample"),
                         reml = REMLcrit(fit), denDf = denDf, cdf = cdf)
  }

  prots <- names(fits)
  s2 <- sapply(prots, function(q) fits[[q]]$sig2); dd <- sapply(prots, function(q) fits[[q]]$denDf)
  sv <- squeezeVar(s2, dd, legacy = TRUE)
  totalDf <- sum(dd)
  out <- NULL
  for (i in seq_along(prots)) {
    q <- prots[i]; F <- fits[[q]]
    post <- sv$var.post[i]
    for (k in seq_len(nrow(sc$contrasts))) {
      cc <- sc$contrasts[k, ]
      est <- sum(cc * F$beta)
      se <- sqrt(drop(t(cc) %*% F$unsc %*% cc) * post)
      dfp <- min(F$cdf[k] + sv$df.prior, totalDf)
      t <- est / se
      qt975 <- qt(0.975, dfp)
      out <- rbind(out, data.frame(
        protein = q, contrast = rownames(sc$contrasts)[k], estimate = fmt(est),
        se_unmoderated = fmt(sqrt(drop(t(cc) %*% F$unsc %*% cc) * F$sig2)), df_satterthwaite = fmt(F$cdf[k]),
        se = fmt(se), df = fmt(dfp), t = fmt(t), p = fmt(2 * pt(-abs(t), dfp)),
        ci_low = fmt(est - qt975 * se), ci_high = fmt(est + qt975 * se),
        sigma2 = fmt(F$sig2), var_individual = fmt(F$ind), var_sample = fmt(F$sam), reml = fmt(F$reml),
        variance_df = fmt(F$denDf), s2_post = fmt(post), df_prior = fmt(sv$df.prior), s2_prior = fmt(sv$var.prior)))
    }
  }
  # BH per contrast over the fitted proteins
  out$adj_p <- NA
  for (cn in unique(out$contrast)) {
    w <- out$contrast == cn
    out$adj_p[w] <- fmt(p.adjust(as.numeric(out$p[w]), "BH"))
  }
  put(out, sprintf("lmm_%s_results.tsv", sc$name))
  unscRows <- do.call(rbind, lapply(prots, function(q) data.frame(protein = q, t(fmt(as.vector(fits[[q]]$unsc))))))
  put(unscRows, sprintf("lmm_%s_unscaled.tsv", sc$name))
  cat(sc$name, ": ", length(prots), " of ", P, " proteins fitted; df.prior ", sv$df.prior, "\n", sep = "")
}

run(scenario("nested"), 40)
run(scenario("sample_only"), 40)

writeLines(c(
  sprintf("generated_utc %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R %s", R.version.string),
  sprintf("lme4 %s; lmerTest %s; limma %s", packageVersion("lme4"), packageVersion("lmerTest"), packageVersion("limma")),
  sprintf("seed %d", seed),
  "script make_nested_lmm_fixtures.R",
  "lmer REML, bobyqa rhoend 1e-12; Satterthwaite via lmerTest::contest; squeezeVar(legacy = TRUE) on sigma^2 with the joint condition F-test DenDF"
), file.path(here, "PROVENANCE_nested_lmm.txt"))
