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
# Data: the M3a scenarios, with about 6% of values pushed 2-4 log2 units away (outliers).

suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(limma); library(MASS) })
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

simulate <- function(sc, P) {
  n <- nrow(sc$samples)
  rows <- NULL
  for (p in seq_len(P)) {
    k <- sample(3:6, 1)
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
robustFit <- function(f, d) {
  model <- lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl)
  sseOld <- model@devcomp$cmp["pwrss"]
  it <- 0; converged <- FALSE
  while (it < maxitRob) {
    it <- it + 1
    res <- resid(model)
    d$w <- MASS::psi.huber(res / mad(res, 0))
    model <- lmerTest::lmer(f, data = d, weights = w, REML = TRUE, control = ctrl)
    sse <- model@devcomp$cmp["pwrss"]
    if (abs(sseOld - sse) / sseOld <= tolRob) { converged <- TRUE; break }
    sseOld <- sse
  }
  list(model = model, weights = d$w, iterations = it, converged = converged)
}
run <- function(sc, P) {
  long <- simulate(sc, P)
  put(data.frame(sc$samples, sc$common, check.names = FALSE), sprintf("rlmm_%s_samples.tsv", sc$name))
  put(data.frame(protein = long$protein, peptide = long$peptide, sample = long$sample, y = fmt(long$y)),
      sprintf("rlmm_%s_values.tsv", sc$name))
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
      "; converged ", sum(sapply(prots, function(q) fits[[q]]$converged)), "; df.prior ", sv$df.prior, "\n", sep = "")
}

run(scenario("nested"), 30)
run(scenario("sample_only"), 30)

writeLines(c(
  sprintf("generated_utc %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R %s", R.version.string),
  sprintf("lme4 %s; lmerTest %s; limma %s; MASS %s", packageVersion("lme4"), packageVersion("lmerTest"),
          packageVersion("limma"), packageVersion("MASS")),
  sprintf("seed %d", seed),
  "script make_robust_lmm_fixtures.R",
  "robust loop = msqrob2 1.20.0 .robust_fitting rules (psi.huber k = 1.345 on resid / mad(res, 0); pwrss tol 1e-6), cap 100, each round refitted from scratch with weights in the call; then lmerTest + squeezeVar(legacy) as GR-21"
), file.path(here, "PROVENANCE_robust_lmm.txt"))
