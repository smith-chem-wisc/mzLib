# Reference values for NestedMixedModel's ROBUST option (STAT1 milestone M3b; QuantProject GR-22):
# msqrob2 1.20.0's robust reweighting (.robust_fitting), run to convergence, then Satterthwaite df and the
# MSstatsTMT-style moderation of GR-21 on the final weighted fit.
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   Rscript make_robust_lmm_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# The loop follows msqrob2's rules:
#   weights <- MASS::psi.huber(resid(model) / mad(resid(model), 0)); refit; stop when |pwrss_old - pwrss| / pwrss_old <= 1e-6
# with the iteration cap raised from msqrob2's default of 1 to 100 (GR-22: iterate to convergence), and each refit done
# from scratch (see robustFit for why).
# Reweighting that does not settle within the cap falls back to msqrob2's default, ONE round (GR-22 as amended
# 2026-11-03): the weights from the unweighted fit's residuals, one weighted fit, then lmerTest and moderation on it.
# Data: the M3a scenarios, with about 6% of values pushed 2-4 log2 units away (outliers); and two fixtures that never
# settle: "fallback" (2-peptide proteins on 4 samples, STAT1's protein model) and "set" (one SETTEST1 set model).

suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(limma); library(MASS); library(msqrob2) })
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", ifelse(is.infinite(v), ifelse(v > 0, "Inf", "-Inf"), sprintf("%.17g", v)))
put <- function(df, name) write.table(df, file.path(here, name), sep = "\t", quote = FALSE, row.names = FALSE)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)
maxitRob <- 100; tolRob <- 1e-6

seed <- 20261009
set.seed(seed)

scenario <- function(name) {
  if (name == "nested") {
    samples <- data.frame(sample = sprintf("s%02d", 1:12), individual = sprintf("m%d", rep(1:6, each = 2)),
                          condition = rep(c("A", "B"), 6), genotype = rep(c("wt", "ko"), each = 6))
    common <- cbind(intercept = 1, conditionB = as.numeric(samples$condition == "B"),
                    genotypeko = as.numeric(samples$genotype == "ko"))
    contrasts <- rbind(conditionB = c(0, 1, 0), genotypeko = c(0, 0, 1), sum = c(0, 1, 1))
    random <- "(1 | individual) + (1 | individual:sample)"
  } else if (name == "fallback") {
    # STAT1's protein model on the smallest design it fits: 2 + 2 samples, one random intercept per sample.
    samples <- data.frame(sample = sprintf("s%02d", 1:4), individual = NA, condition = rep(c("A", "B"), each = 2))
    common <- cbind(intercept = 1, conditionB = as.numeric(samples$condition == "B"))
    contrasts <- rbind(conditionB = c(0, 1))
    colnames(contrasts) <- colnames(common)
    return(list(name = name, samples = samples, common = common, conditionRows = list(c(0, 1)),
                contrasts = contrasts, random = "(1 | sample)"))
  } else {
    samples <- data.frame(sample = sprintf("s%02d", 1:12), individual = NA,
                          group = rep(c("A", "B", "C"), each = 4))
    common <- cbind(intercept = 1, groupB = as.numeric(samples$group == "B"), groupC = as.numeric(samples$group == "C"))
    contrasts <- rbind(groupB = c(0, 1, 0), groupC = c(0, 0, 1), C_vs_B = c(0, -1, 1))
    random <- "(1 | sample)"
  }
  colnames(contrasts) <- colnames(common)
  list(name = name, samples = samples, common = common, conditionRows = list(c(0, 1, 0), c(0, 0, 1)),
       contrasts = contrasts, random = random)
}

simulate <- function(sc, P, peptides = NULL) {
  n <- nrow(sc$samples)
  rows <- NULL
  for (p in seq_len(P)) {
    k <- if (is.null(peptides)) sample(3:6, 1) else peptides
    beta <- c(rnorm(1, 22, 1.5), rnorm(ncol(sc$common) - 1, 0, 0.5))
    pep <- c(0, rnorm(k - 1, 0, 1))
    tauInd <- c(0.15, 0.4, 0.25)[(p %% 3) + 1]
    tauSam <- c(0.3, 0.1, 0.2)[(p %% 3) + 1]
    sigma <- sqrt(0.08 * 4 / rchisq(1, 4))
    ind <- if (all(is.na(sc$samples$individual))) rep(0, n) else rnorm(6, 0, tauInd)[as.integer(factor(sc$samples$individual))]
    sam <- rnorm(n, 0, tauSam)
    for (j in seq_len(k)) {
      y <- drop(sc$common %*% beta) + pep[j] + ind + sam + rnorm(n, 0, sigma)
      out <- runif(n) < 0.06
      y[out] <- y[out] + sample(c(-1, 1), sum(out), TRUE) * runif(sum(out), 2, 4)
      y[runif(n) < 0.05] <- NA
      rows <- rbind(rows, data.frame(protein = sprintf("p%03d", p), peptide = sprintf("pep%d", j),
                                     sample = sc$samples$sample, y = y))
    }
  }
  rows
}

# msqrob2's loop, with one change: each round refits FROM SCRATCH (lmer with the weights in the call) instead of
# lme4::refit, which warm-starts from the previous theta. On this data the warm start stayed stuck at theta = 0 on two
# proteins (sample_only p016: REML 57.664840794 vs 57.664346682 fitted afresh; p023: 51.287 vs 51.092), so its later
# weights were built on a sub-optimal fit. The weight rule and the stopping rule are msqrob2's unchanged. A fresh fit
# with weights in the call is also what lmerTest needs: it rebuilds its deviance from the call, and weights injected
# into the frame (as msqrob2 does) are invisible to it (p001: df 171 on 44 observations).
#
# Reweighting that does not settle within maxitRob rounds falls back to ONE round (GR-22 as amended 2026-11-03), which is
# what msqrob2 does by default (maxitRob = 1): weights from the UNWEIGHTED fit's residuals, then one weighted fit, made
# here from scratch with those weights in the call. That those weights are msqrob2's own is checked against
# msqrob2:::.robust_fitting(model, maxitRob = 1) on the same unweighted lme4 fit (oneRoundCheck).
robustFit <- function(f, d) {
  model <- lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl)
  sseOld <- model@devcomp$cmp["pwrss"]
  it <- 0; converged <- FALSE; first <- NULL
  while (it < maxitRob) {
    it <- it + 1
    res <- resid(model)
    d$w <- MASS::psi.huber(res / mad(res, 0))
    if (it == 1) first <- d$w
    model <- lmerTest::lmer(f, data = d, weights = w, REML = TRUE, control = ctrl)
    sse <- model@devcomp$cmp["pwrss"]
    if (abs(sseOld - sse) / sseOld <= tolRob) { converged <- TRUE; break }
    sseOld <- sse
  }
  if (!converged) {
    d$w <- first
    model <- lmerTest::lmer(f, data = d, weights = w, REML = TRUE, control = ctrl)
    oneRoundCheck(f, d, first)
  }
  list(model = model, weights = d$w, iterations = it, converged = converged)
}
# msqrob2 1.20.0 with maxitRob = 1, as .noridge_msqrobLmer calls it: an lme4 fit whose frame weights are set to 1, then
# .robust_fitting(model, 1, pwrss, tol). Its weights must be the fallback's exactly. Its own refit warm-starts (lme4::refit),
# so its fixed effects are reported beside the fresh fit's, not compared.
oneRoundCheck <- function(f, d, weights) {
  m0 <- lme4::lmer(f, data = d, REML = TRUE, control = ctrl)
  m0@frame$`(weights)` <- rep(1, nrow(m0@frame))
  ms <- msqrob2:::.robust_fitting(m0, 1, m0@devcomp$cmp["pwrss"], tolRob)
  diff <- max(abs(ms@frame$`(weights)` - weights))
  stopifnot(diff < 1e-12)
  d$w <- weights
  fresh <- lme4::lmer(f, data = d, weights = w, REML = TRUE, control = ctrl)
  oneRound[[length(oneRound) + 1]] <<- c(weights = diff, fixef = max(abs(fixef(ms) - fixef(fresh))))
}
oneRound <- list()
run <- function(sc, P, peptides = NULL, keep = NULL) {
  long <- simulate(sc, P, peptides)
  put(data.frame(sc$samples, sc$common, check.names = FALSE), sprintf("rlmm_%s_samples.tsv", sc$name))
  put(data.frame(contrast = rownames(sc$contrasts), sc$contrasts, check.names = FALSE), sprintf("rlmm_%s_contrasts.tsv", sc$name))

  fits <- list(); weightRows <- NULL
  for (prot in unique(long$protein)) {
    d <- merge(long[long$protein == prot & !is.na(long$y), ], sc$samples, by = "sample")
    d <- d[order(d$peptide, d$sample, method = "radix"), ]
    d$peptide <- factor(d$peptide, levels = sort(unique(d$peptide), method = "radix"))
    for (cn in colnames(sc$common)[-1]) d[[cn]] <- sc$common[match(d$sample, sc$samples$sample), cn]
    terms <- paste(colnames(sc$common)[-1], collapse = " + ")
    f <- as.formula(sprintf("y ~ %s + peptide + %s", terms, sc$random))
    rb <- tryCatch(robustFit(f, d), error = function(e) { message(prot, ": ", conditionMessage(e)); NULL })
    if (is.null(rb)) next
    fit <- rb$model
    p <- ncol(sc$common)
    vc <- as.data.frame(VarCorr(fit))
    getv <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v else NA }
    sig2 <- sigma(fit)^2
    L <- matrix(0, length(sc$conditionRows), length(fixef(fit)))
    for (r in seq_along(sc$conditionRows)) L[r, 1:p] <- sc$conditionRows[[r]]
    cdf <- sapply(seq_len(nrow(sc$contrasts)), function(i) {
      l <- numeric(length(fixef(fit))); l[1:p] <- sc$contrasts[i, ]; contest(fit, l, joint = FALSE)$df })
    fits[[prot]] <- list(beta = fixef(fit)[1:p], unsc = as.matrix(vcov(fit))[1:p, 1:p] / sig2, sig2 = sig2,
                         ind = getv("individual"), sam = if (sc$name == "nested") getv("individual:sample") else getv("sample"),
                         reml = REMLcrit(fit), denDf = contest(fit, L, joint = TRUE)$DenDF, cdf = cdf,
                         iterations = rb$iterations, converged = rb$converged)
    weightRows <- rbind(weightRows, data.frame(protein = prot, peptide = as.character(d$peptide), sample = d$sample,
                                               weight = fmt(rb$weights)))
  }
  if (!is.null(keep)) {
    kept <- keep(fits)
    fits <- fits[kept]; long <- long[long$protein %in% kept, ]; weightRows <- weightRows[weightRows$protein %in% kept, ]
  }
  put(data.frame(protein = long$protein, peptide = long$peptide, sample = long$sample, y = fmt(long$y)),
      sprintf("rlmm_%s_values.tsv", sc$name))
  put(weightRows, sprintf("rlmm_%s_weights.tsv", sc$name))

  prots <- names(fits)
  s2 <- sapply(prots, function(q) fits[[q]]$sig2); dd <- sapply(prots, function(q) fits[[q]]$denDf)
  sv <- squeezeVar(s2, dd, legacy = TRUE)
  totalDf <- sum(dd)
  out <- NULL
  for (i in seq_along(prots)) {
    q <- prots[i]; F <- fits[[q]]; post <- sv$var.post[i]
    for (k in seq_len(nrow(sc$contrasts))) {
      cc <- sc$contrasts[k, ]
      est <- sum(cc * F$beta)
      se <- sqrt(drop(t(cc) %*% F$unsc %*% cc) * post)
      dfp <- min(F$cdf[k] + sv$df.prior, totalDf)
      t <- est / se
      qt975 <- qt(0.975, dfp)
      out <- rbind(out, data.frame(
        protein = q, contrast = rownames(sc$contrasts)[k], estimate = fmt(est), df_satterthwaite = fmt(F$cdf[k]),
        se = fmt(se), df = fmt(dfp), t = fmt(t), p = fmt(2 * pt(-abs(t), dfp)),
        ci_low = fmt(est - qt975 * se), ci_high = fmt(est + qt975 * se),
        sigma2 = fmt(F$sig2), var_individual = fmt(F$ind), var_sample = fmt(F$sam), reml = fmt(F$reml),
        variance_df = fmt(F$denDf), s2_post = fmt(post), df_prior = fmt(sv$df.prior),
        iterations = F$iterations, converged = F$converged))
    }
  }
  out$adj_p <- NA
  for (cn in unique(out$contrast)) { w <- out$contrast == cn; out$adj_p[w] <- fmt(p.adjust(as.numeric(out$p[w]), "BH")) }
  put(out, sprintf("rlmm_%s_results.tsv", sc$name))
  put(do.call(rbind, lapply(prots, function(q) data.frame(protein = q, t(fmt(as.vector(fits[[q]]$unsc)))))),
      sprintf("rlmm_%s_unscaled.tsv", sc$name))
  its <- sapply(prots, function(q) fits[[q]]$iterations)
  cat(sc$name, ": ", length(prots), " of ", P, " fitted; iterations ", min(its), "-", max(its),
      "; converged ", sum(sapply(prots, function(q) fits[[q]]$converged)), " (the rest fell back to one round)",
      "; df.prior ", sv$df.prior, "\n", sep = "")
}

run(scenario("nested"), 30)
run(scenario("sample_only"), 30)
# 2-peptide proteins on 4 samples, the case where reweighting most often never settles (QuantProject's STAT1 measurement,
# 2026-11-03). Kept: the first 4 of 600 that never settle and the first 6 that settle within 40 rounds, so that the
# moderation mixes both kinds and no kept protein sits near the cap, where R and C# could fall either side of it.
run(scenario("fallback"), 600, peptides = 2, keep = function(fits) {
  its <- sapply(fits, function(x) x$iterations); conv <- sapply(fits, function(x) x$converged)
  cat("fallback: ", sum(!conv), " of ", length(fits), " never settle\n", sep = "")
  stopifnot(sum(!conv) >= 2)
  sort(c(head(names(fits)[!conv], 4), head(names(fits)[conv & its <= 40], 6)), method = "radix")
})

# "set": the set model of SETTEST1 S2's AnIncompleteMember_StaysInTheEstimate_AndBlocksSupported (mzLib
# Test/Quantification/Differential/SetTestModelTests.cs). Three members of two peptides each over o1, o2 (old) and y1, y2
# (young), the third member missing at o2; each second peptide 0.4 higher with a fixed jitter. Fixed: per member an
# intercept, old, and its second peptide; random: (1 | sample) + (1 | sample:member), which NestedMixedModel fits as
# individual = sample, sample = member. Reweighting never settles (R cycles for all 100 rounds), so it falls back.
setRun <- function() {
  jitter <- c(0.03, -0.05, 0.02, -0.01, 0.04, -0.03)
  four <- c("o1", "o2", "y1", "y2")
  shifted <- function(s) c(20 + s, 21 + s, 20.5, 20.7)
  vals <- list(G1 = shifted(0.05), G2 = shifted(-0.05), G3 = c(20, NA, 20.5, 20.7))
  rows <- NULL
  for (k in seq_along(vals)) for (pp in 1:2) for (i in 1:4) {
    v <- if (pp == 1) vals[[k]][i] else vals[[k]][i] + 0.4 + jitter[(i - 1 + k - 1) %% 6 + 1]
    if (!is.na(v)) rows <- rbind(rows, data.frame(member = names(vals)[k],
                                                  peptide = paste0(names(vals)[k], "PEP", c("A", "B")[pp]), sample = four[i], y = v))
  }
  g <- match(rows$member, names(vals))
  X <- matrix(0, nrow(rows), 9)
  for (r in seq_len(nrow(rows))) {
    X[r, 2 * g[r] - 1] <- 1
    X[r, 2 * g[r]] <- as.numeric(startsWith(rows$sample[r], "o"))
    if (endsWith(rows$peptide[r], "B")) X[r, 6 + g[r]] <- 1
  }
  common <- paste0(rep(names(vals), each = 2), c(":intercept", ":old"))
  colnames(X) <- c(paste0("common_", common), paste0("nuisance_", names(vals), ":PEPB"))
  d <- rows; d$sm <- paste(d$sample, d$member, sep = ":"); d$X <- X
  rb <- robustFit(y ~ 0 + X + (1 | sample) + (1 | sm), d)
  stopifnot(!rb$converged, rb$iterations == maxitRob)
  fit <- rb$model
  contrasts <- rbind(old_mean = c(0, 1, 0, 1, 0, 1) / 3, G1_old = c(0, 1, 0, 0, 0, 0), G3_old = c(0, 0, 0, 0, 0, 1))
  colnames(contrasts) <- common
  vc <- as.data.frame(VarCorr(fit))
  sig2 <- sigma(fit)^2
  out <- do.call(rbind, lapply(rownames(contrasts), function(cn) {
    l <- c(contrasts[cn, ], 0, 0, 0)
    ct <- contest(fit, l, joint = FALSE)
    data.frame(contrast = cn, estimate = fmt(ct$Estimate), se = fmt(ct$`Std. Error`), df_satterthwaite = fmt(ct$df),
               sigma2 = fmt(sig2), var_sample = fmt(vc$vcov[vc$grp == "sample"]), var_sample_member = fmt(vc$vcov[vc$grp == "sm"]),
               reml = fmt(REMLcrit(fit)), iterations = rb$iterations, converged = rb$converged)
  }))
  put(data.frame(rows[, c("member", "peptide", "sample")], y = fmt(rows$y), apply(X, 2, fmt), check.names = FALSE),
      "rlmm_set_values.tsv")
  put(data.frame(contrast = rownames(contrasts), contrasts, check.names = FALSE), "rlmm_set_contrasts.tsv")
  put(data.frame(rows[, c("member", "peptide", "sample")], weight = fmt(rb$weights)), "rlmm_set_weights.tsv")
  put(out, "rlmm_set_results.tsv")
  cat("set: ", nrow(rows), " values; ", rb$iterations, " rounds, not settled; one round down-weights ", sum(rb$weights < 1),
      ", min weight ", min(rb$weights), "; old_mean ", out$estimate[1], " df ", out$df_satterthwaite[1], "\n", sep = "")
}
setRun()

checks <- do.call(rbind, oneRound)
cat("one-round fallbacks checked against msqrob2:::.robust_fitting(maxitRob = 1): ", nrow(checks),
    "; largest weight difference ", max(checks[, "weights"]), "; largest fixef difference of msqrob2's warm-started refit vs the fresh fit ",
    max(checks[, "fixef"]), "\n", sep = "")

writeLines(c(
  sprintf("generated_utc %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R %s", R.version.string),
  sprintf("lme4 %s; lmerTest %s; limma %s; MASS %s; msqrob2 %s", packageVersion("lme4"), packageVersion("lmerTest"),
          packageVersion("limma"), packageVersion("MASS"), packageVersion("msqrob2")),
  sprintf("seed %d", seed),
  "script make_robust_lmm_fixtures.R",
  "robust loop = msqrob2 1.20.0 .robust_fitting rules (psi.huber k = 1.345 on resid / mad(res, 0); pwrss tol 1e-6), cap 100, each round refitted from scratch with weights in the call; then lmerTest + squeezeVar(legacy) as GR-21",
  "not settled within the cap (rlmm_fallback_*, rlmm_set_*) = one round, as msqrob2's default maxitRob = 1 (GR-22 as amended 2026-11-03): weights from the unweighted fit's residuals, identical to msqrob2:::.robust_fitting(maxitRob = 1)'s, then a fresh weighted fit with those weights in the call",
  sprintf("one-round fallbacks: %d, each checked against msqrob2:::.robust_fitting; largest weight difference %.3g; msqrob2's warm-started refit vs the fresh fit, largest fixef difference %.3g",
          nrow(checks), max(checks[, "weights"]), max(checks[, "fixef"]))
), file.path(here, "PROVENANCE_robust_lmm.txt"))
