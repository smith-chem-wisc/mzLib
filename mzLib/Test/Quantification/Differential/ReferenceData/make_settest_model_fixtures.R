# Reference values for SETTEST1's set model (QuantProject DEF-DIFF-SET-MODEL, DATA-DEFINITIONS v3.12):
#   log2 value ~ peptide + member:(the stratum's fixed design) + (1 | sample) + (1 | sample:member), REML,
# with beta-bar the equal-weight average of the members' contrasts, its se from the fit and its df Satterthwaite's.
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   E:/Programs/R-4.6.1/bin/Rscript.exe make_settest_model_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# The fixed part is built as an explicit matrix, so its coding is the C# one: per member, an intercept and the design's
# other columns (member:(design)); per member, one indicator for each of its peptides but the first in ordinal order.
#
# Robust weighting (S2b; GR-8, GR-22): the same data refitted with msqrob2 1.20.0's reweighting run to convergence, with
# the loop of StatisticalModels/ReferenceData/make_robust_lmm_fixtures.R (#1499) unchanged:
#   weights <- MASS::psi.huber(resid(model) / mad(resid(model), 0)); refit; stop when |pwrss_old - pwrss| / pwrss_old <= 1e-6
# with the round cap raised from msqrob2's 1 to 100, and each round refitted from scratch with the weights IN THE lmer CALL
# (lmerTest rebuilds its deviance from the call; weights put in the frame afterwards are invisible to it).
# Reweighting that does not settle within the cap falls back to msqrob2's default, ONE round (GR-22 as amended
# 2026-11-03; #1499's script does the same): weights from the unweighted fit's residuals, one weighted fit. The
# "incomplete" set below never settles.

suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(MASS); library(msqrob2) })
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
put <- function(df, name) write.table(df, file.path(here, name), sep = "\t", quote = FALSE, row.names = FALSE)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)
maxitRob <- 100; tolRob <- 1e-6

# msqrob2's loop, as make_robust_lmm_fixtures.R (#1499) runs it: each round refits FROM SCRATCH with the weights in the call
# (lme4::refit warm-starts from the previous theta and can stay stuck at theta = 0). Not settled within the cap: one round,
# refitted with the first round's weights in the call, whose weights must equal msqrob2:::.robust_fitting(maxitRob = 1)'s.
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
    # msqrob2 1.20.0 with maxitRob = 1, as .noridge_msqrobLmer calls it, on the same unweighted lme4 fit.
    m0 <- lme4::lmer(f, data = d, REML = TRUE, control = ctrl)
    m0@frame$`(weights)` <- rep(1, nrow(m0@frame))
    ms <- msqrob2:::.robust_fitting(m0, 1, m0@devcomp$cmp["pwrss"], tolRob)
    stopifnot(max(abs(ms@frame$`(weights)` - first)) < 1e-12)
    cat("one round: weights equal msqrob2's (largest difference ", max(abs(ms@frame$`(weights)` - first)), ")\n", sep = "")
  }
  list(model = model, weights = d$w, iterations = it, converged = converged)
}

seed <- 20261026
set.seed(seed)

scenario <- function(name) {
  if (name == "balanced") {
    samples <- data.frame(sample = sprintf("s%02d", 1:8), group = rep(c("young", "old"), each = 4), treatment = "none")
    peptides <- c(3, 2, 4, 3); missing <- 0; context <- FALSE
  } else if (name == "missing") {
    samples <- data.frame(sample = sprintf("s%02d", 1:10), group = rep(c("young", "old"), 5), treatment = "none")
    peptides <- c(3, 4, 2); missing <- 0.12; context <- FALSE
  } else {
    samples <- data.frame(sample = sprintf("s%02d", 1:12), group = rep(c("young", "old"), 6),
                          treatment = rep(c("a", "b"), each = 6))
    peptides <- c(1, 2, 3, 4, 5); missing <- 0.05; context <- TRUE
  }
  list(name = name, samples = samples, peptides = peptides, missing = missing, context = context)
}

simulate <- function(sc) {
  n <- nrow(sc$samples); old <- as.numeric(sc$samples$group == "old"); trt <- as.numeric(sc$samples$treatment == "b")
  shared <- rnorm(n, 0, 0.35)                                       # (1 | sample): every member moves together
  rows <- NULL
  for (m in seq_along(sc$peptides)) {
    mu <- rnorm(1, 22, 1.2); effect <- rnorm(1, 0.6, 0.4); te <- if (sc$context) rnorm(1, 0, 0.5) else 0
    own <- rnorm(n, 0, 0.25)                                        # (1 | sample:member)
    for (k in seq_len(sc$peptides[m])) {
      pep <- if (k == 1) 0 else rnorm(1, 0, 0.8)
      y <- mu + effect * old + te * trt + pep + shared + own + rnorm(n, 0, 0.2)
      if (sc$missing > 0) y[runif(n) < sc$missing] <- NA
      rows <- rbind(rows, data.frame(member = sprintf("G%d", m), peptide = sprintf("G%dpep%d", m, k),
                                     sample = sc$samples$sample, y = y))
    }
  }
  if (sc$name == "missing") rows$y[rows$member == "G2" & rows$sample == "s03"] <- NA   # one member gone from one sample
  rows[!is.na(rows$y), ]
}

run <- function(sc) {
  d <- simulate(sc)
  d <- merge(d, sc$samples, by = "sample")
  d <- d[order(d$member, d$peptide, d$sample, method = "radix"), ]
  members <- sort(unique(d$member), method = "radix")
  X <- NULL; names <- NULL
  for (m in members) {
    is <- as.numeric(d$member == m)
    X <- cbind(X, is, is * (d$group == "old")); names <- c(names, paste0(m, ":intercept"), paste0(m, ":old"))
    if (sc$context) { X <- cbind(X, is * (d$treatment == "b")); names <- c(names, paste0(m, ":treatmentb")) }
  }
  for (m in members) {
    peps <- sort(unique(d$peptide[d$member == m]), method = "radix")
    for (p in peps[-1]) { X <- cbind(X, as.numeric(d$peptide == p)); names <- c(names, p) }
  }
  colnames(X) <- names
  d$sm <- paste(d$sample, d$member, sep = ":")
  fit <- lmerTest::lmer(y ~ 0 + X + (1 | sample) + (1 | sm), data = d, REML = TRUE, control = ctrl)
  L <- numeric(ncol(X)); L[grep(":old$", names)] <- 1 / length(members)
  ct <- contest(fit, L, joint = FALSE)
  vc <- as.data.frame(VarCorr(fit))
  put(data.frame(sample = sc$samples$sample, group = sc$samples$group, treatment = sc$samples$treatment),
      sprintf("settest_model_%s_samples.tsv", sc$name))
  put(data.frame(member = d$member, peptide = d$peptide, sample = d$sample, y = fmt(d$y)),
      sprintf("settest_model_%s_values.tsv", sc$name))
  put(data.frame(estimate = fmt(ct$Estimate), se = fmt(ct$`Std. Error`), df = fmt(ct$df),
                 sigma2 = fmt(sigma(fit)^2), var_sample = fmt(vc$vcov[vc$grp == "sample"]),
                 var_sample_member = fmt(vc$vcov[vc$grp == "sm"]), reml = fmt(REMLcrit(fit)),
                 members = length(members), peptides = length(unique(d$peptide)), observations = nrow(d)),
      sprintf("settest_model_%s_result.tsv", sc$name))
  cat(sc$name, ": estimate ", ct$Estimate, " se ", ct$`Std. Error`, " df ", ct$df, "\n", sep = "")

  # The same model, robust.
  rb <- robustFit(y ~ 0 + X + (1 | sample) + (1 | sm), d)
  stopifnot(rb$converged)
  rfit <- rb$model
  rct <- contest(rfit, L, joint = FALSE)
  # Guard: lmerTest saw the weights. A fresh fit here, with the final weights in the call, gives the same df.
  d$w <- rb$weights
  check <- lmerTest::lmer(y ~ 0 + X + (1 | sample) + (1 | sm), data = d, weights = w, REML = TRUE, control = ctrl)
  stopifnot(isTRUE(all.equal(contest(check, L, joint = FALSE)$df, rct$df, tolerance = 1e-8)))
  rvc <- as.data.frame(VarCorr(rfit))
  put(data.frame(member = d$member, peptide = d$peptide, sample = d$sample, weight = fmt(rb$weights)),
      sprintf("settest_model_%s_robust_weights.tsv", sc$name))
  put(data.frame(estimate = fmt(rct$Estimate), se = fmt(rct$`Std. Error`), df = fmt(rct$df),
                 sigma2 = fmt(sigma(rfit)^2), var_sample = fmt(rvc$vcov[rvc$grp == "sample"]),
                 var_sample_member = fmt(rvc$vcov[rvc$grp == "sm"]), reml = fmt(REMLcrit(rfit)),
                 iterations = rb$iterations, converged = rb$converged, down_weighted = sum(rb$weights < 1),
                 members = length(members), peptides = length(unique(d$peptide)), observations = nrow(d)),
      sprintf("settest_model_%s_robust_result.tsv", sc$name))
  cat(sc$name, " robust: estimate ", rct$Estimate, " se ", rct$`Std. Error`, " df ", rct$df, "; ", rb$iterations,
      " rounds; ", sum(rb$weights < 1), " of ", nrow(d), " values down-weighted, min weight ", min(rb$weights), "\n", sep = "")
}

for (s in c("balanced", "missing", "context")) run(scenario(s))

# "incomplete": the set of SetTestModelTests' IncompleteG3, built as its Members() helper builds it (no random numbers):
# over o1, o2 (old) and y1, y2 (young), G1 = 20.05, 21.05, 20.5, 20.7; G2 = 19.95, 20.95, 20.5, 20.7; G3 = 20, -, 20.5,
# 20.7; each member's second peptide = the first + 0.4 + jitter[(i + k) mod 6] (sample i, member k, both from 0). Robust
# only: its reweighting never settles (R cycles for all 100 rounds), so its fit is the one round (GR-22 amended 2026-11-03).
incomplete <- function() {
  jitter <- c(0.03, -0.05, 0.02, -0.01, 0.04, -0.03)
  four <- c("o1", "o2", "y1", "y2")
  shifted <- function(s) c(20 + s, 21 + s, 20.5, 20.7)
  vals <- list(G1 = shifted(0.05), G2 = shifted(-0.05), G3 = c(20, NA, 20.5, 20.7))
  d <- NULL
  for (k in seq_along(vals)) for (pp in c("A", "B")) for (i in 1:4) {
    v <- if (pp == "A") vals[[k]][i] else vals[[k]][i] + 0.4 + jitter[(i - 1 + k - 1) %% 6 + 1]
    if (!is.na(v)) d <- rbind(d, data.frame(member = names(vals)[k], peptide = paste0(names(vals)[k], "PEP", pp),
                                            sample = four[i], y = v, group = if (i <= 2) "old" else "young"))
  }
  members <- names(vals)
  X <- NULL; names <- NULL
  for (m in members) {
    is <- as.numeric(d$member == m)
    X <- cbind(X, is, is * (d$group == "old")); names <- c(names, paste0(m, ":intercept"), paste0(m, ":old"))
  }
  for (m in members) { X <- cbind(X, as.numeric(d$peptide == paste0(m, "PEPB"))); names <- c(names, paste0(m, "PEPB")) }
  colnames(X) <- names
  d$sm <- paste(d$sample, d$member, sep = ":")
  L <- numeric(ncol(X)); L[grep(":old$", names)] <- 1 / length(members)
  rb <- robustFit(y ~ 0 + X + (1 | sample) + (1 | sm), d)
  stopifnot(!rb$converged, rb$iterations == maxitRob)
  rfit <- rb$model
  rct <- contest(rfit, L, joint = FALSE)
  rvc <- as.data.frame(VarCorr(rfit))
  put(data.frame(member = d$member, peptide = d$peptide, sample = d$sample, y = fmt(d$y)), "settest_model_incomplete_values.tsv")
  put(data.frame(member = d$member, peptide = d$peptide, sample = d$sample, weight = fmt(rb$weights)),
      "settest_model_incomplete_robust_weights.tsv")
  put(data.frame(estimate = fmt(rct$Estimate), se = fmt(rct$`Std. Error`), df = fmt(rct$df),
                 sigma2 = fmt(sigma(rfit)^2), var_sample = fmt(rvc$vcov[rvc$grp == "sample"]),
                 var_sample_member = fmt(rvc$vcov[rvc$grp == "sm"]), reml = fmt(REMLcrit(rfit)),
                 iterations = rb$iterations, converged = rb$converged, down_weighted = sum(rb$weights < 1),
                 members = length(members), peptides = length(unique(d$peptide)), observations = nrow(d)),
      "settest_model_incomplete_robust_result.tsv")
  cat("incomplete robust: not settled in ", rb$iterations, " rounds; one round: estimate ", rct$Estimate, " se ",
      rct$`Std. Error`, " df ", rct$df, "; ", sum(rb$weights < 1), " of ", nrow(d), " values down-weighted, min weight ",
      min(rb$weights), "\n", sep = "")
}
incomplete()

writeLines(c(
  sprintf("generated_utc %s UTC", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
  sprintf("R %s", R.version.string),
  sprintf("lme4 %s; lmerTest %s; MASS %s; msqrob2 %s", packageVersion("lme4"), packageVersion("lmerTest"),
          packageVersion("MASS"), packageVersion("msqrob2")),
  sprintf("seed %d", seed),
  "script make_settest_model_fixtures.R",
  "lmer REML on an explicit fixed matrix (member:(intercept, old[, treatment b]) + peptides within member), (1|sample) + (1|sample:member), bobyqa rhoend 1e-12; beta-bar and its Satterthwaite df via lmerTest::contest",
  "robust (*_robust_*) = msqrob2 1.20.0 .robust_fitting rules (psi.huber k = 1.345 on resid / mad(res, 0); pwrss tol 1e-6), cap 100, each round refitted from scratch with weights in the call, as make_robust_lmm_fixtures.R (#1499); then lmerTest::contest on the final weighted fit",
  "not settled within the cap (settest_model_incomplete_*, hand-made values, no random numbers) = one round, as msqrob2's default maxitRob = 1 (GR-22 as amended 2026-11-03): weights from the unweighted fit's residuals, checked equal to msqrob2:::.robust_fitting(maxitRob = 1)'s, then a fresh weighted fit with those weights in the call"),
  file.path(here, "PROVENANCE_settest_model.txt"))
