# Reference outputs for mzLib's LogisticRegression, MixedModel, SpearmanCorrelation and
# PValueCombination, from R's glm, nlme::lme, cor.test and base distribution functions.
#
# Run ONCE, in CI (a manual-dispatch job on a throwaway branch), and the outputs checked in beside
# this script. R is never a runtime or test dependency of mzLib: the tests read the frozen TSVs.
#
#   Rscript make_glm_lmm_fixtures.R <out_dir>
#
# Settings, and why:
#   * glm(binomial) with epsilon = 1e-14, so R and mzLib both stop at the maximum-likelihood estimate
#     rather than at R's default tolerance. Features where R's fitted probabilities reach 1e-10 of 0 or 1
#     are written as separated (NaN), which is what mzLib reports.
#   * nlme::lme, random intercept per group, REML and ML, with tight optimizer tolerances. nlme's
#     "containment" denominator df are recorded; mzLib implements the same rule.
#   * cor.test(method = "spearman"): exact = TRUE for n <= 9 without ties (enumerated exactly in both),
#     exact = FALSE otherwise (the t approximation in both).
#   * Fisher and Stouffer from pchisq / qnorm / pnorm directly.

suppressPackageStartupMessages(library(nlme))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
                                   row.names = FALSE, na = "NaN")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20260924)

# ---- logistic regression ----------------------------------------------------------------------------
G <- 80; n <- 40
design <- cbind(Intercept = 1, x = rnorm(n), grp = rep(0:1, n / 2))
Y <- matrix(NA_real_, G, n)
for (g in seq_len(G)) {
  beta <- c(rnorm(1, -0.3, 0.6), rnorm(1, 0, 0.8), rnorm(1, 0, 0.8))
  Y[g, ] <- rbinom(n, 1, plogis(drop(design %*% beta)))
}
Y[matrix(runif(G * n) < 0.1, G, n)] <- NA
res <- vector("list", G)
for (g in seq_len(G)) {
  obs <- !is.na(Y[g, ])
  X <- design[obs, , drop = FALSE]; y <- Y[g, obs]
  row <- data.frame(status = "fitted", b0 = NaN, b1 = NaN, b2 = NaN, se0 = NaN, se1 = NaN, se2 = NaN,
                    p1 = NaN, deviance = NaN)
  if (sum(obs) <= ncol(X) || qr(X)$rank < ncol(X)) { row$status <- "not_fitted"; res[[g]] <- row; next }
  fit <- suppressWarnings(glm(y ~ 0 + X, family = binomial,
                              control = glm.control(epsilon = 1e-14, maxit = 100)))
  mu <- fitted(fit)
  if (all(y == y[1]) || any(mu < 1e-10 | mu > 1 - 1e-10) || !fit$converged) {
    row$status <- "separated"; res[[g]] <- row; next
  }
  s <- summary(fit)$coefficients
  row[, c("b0", "b1", "b2")] <- s[, 1]; row[, c("se0", "se1", "se2")] <- s[, 2]
  row$p1 <- s[2, 4]; row$deviance <- fit$deviance
  res[[g]] <- row
}
w(fmt(as.data.frame(Y)), "glm_responses.tsv")
w(fmt(as.data.frame(design)), "glm_design.tsv")
w(fmt(do.call(rbind, res)), "glm_fit.tsv")

# ---- random-intercept mixed model --------------------------------------------------------------------
sizes <- c(3, 5, 4, 6, 2, 4, 5)                    # unbalanced groups, as datasets are
group <- rep(paste0("d", seq_along(sizes)), sizes)
wp <- rep(c(0, 1, 0, 1, 1, 0, 1), sizes)             # constant within a group: a group-level covariate
age <- round(runif(length(group), 20, 80)) / 10 - 5  # varies within groups
N <- length(group); F <- 120
Yl <- matrix(NA_real_, F, N)
for (f in seq_len(F)) {
  tau <- c(0, 0.2, 1)[1 + f %% 3]                    # includes features with no group variance
  u <- rnorm(length(sizes), 0, tau)[as.integer(factor(group, levels = unique(group)))]
  Yl[f, ] <- 0.3 + rnorm(1, 0, 0.3) * age + rnorm(1, 0, 0.5) * wp + u + rnorm(N, 0, 0.5)
}
Yl[matrix(runif(F * N) < 0.08, F, N)] <- NA
ctrl <- lmeControl(maxIter = 500, msMaxIter = 500, niterEM = 100, tolerance = 1e-12, msTol = 1e-14,
                   opt = "nlminb", returnObject = TRUE)
lmm <- list()
for (method in c("REML", "ML")) for (f in seq_len(F)) {
  obs <- !is.na(Yl[f, ])
  d <- data.frame(y = Yl[f, obs], age = age[obs], wp = wp[obs], group = group[obs])
  row <- data.frame(feature = f - 1, method = method, status = "fitted",
                    b0 = NaN, b1 = NaN, b2 = NaN, se0 = NaN, se1 = NaN, se2 = NaN,
                    df0 = NaN, df1 = NaN, df2 = NaN, p1 = NaN, p2 = NaN,
                    sigma2 = NaN, tau2 = NaN, loglik = NaN)
  X <- cbind(1, d$age, d$wp)
  if (length(unique(d$group)) < 2 || max(table(d$group)) < 2 || qr(X)$rank < 3 || nrow(d) <= 3) {
    row$status <- "not_fitted"; lmm[[length(lmm) + 1]] <- row; next
  }
  fit <- tryCatch(lme(y ~ age + wp, random = ~ 1 | group, data = d, method = method, control = ctrl),
                  error = function(e) NULL)
  if (is.null(fit)) { row$status <- "r_error"; lmm[[length(lmm) + 1]] <- row; next }
  tt <- summary(fit)$tTable
  row[, c("b0", "b1", "b2")] <- tt[, "Value"]; row[, c("se0", "se1", "se2")] <- tt[, "Std.Error"]
  row[, c("df0", "df1", "df2")] <- tt[, "DF"]; row$p1 <- tt[2, "p-value"]; row$p2 <- tt[3, "p-value"]
  row$sigma2 <- fit$sigma^2
  row$tau2 <- as.numeric(VarCorr(fit)[1, "Variance"])
  row$loglik <- as.numeric(logLik(fit))
  lmm[[length(lmm) + 1]] <- row
}
w(fmt(as.data.frame(Yl)), "lmm_responses.tsv")
w(fmt(data.frame(Intercept = 1, age = age, wp = wp)), "lmm_design.tsv")
w(data.frame(group = group), "lmm_groups.tsv")
w(fmt(do.call(rbind, lmm)), "lmm_fit.tsv")

# ---- Spearman ----------------------------------------------------------------------------------------
sp_in <- list(); sp_out <- list()
add_case <- function(name, x, y, exact) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = exact))
  sp_in[[length(sp_in) + 1]] <<- data.frame(case = name, x = x, y = y)
  sp_out[[length(sp_out) + 1]] <<- data.frame(case = name, rho = unname(ct$estimate), p = ct$p.value,
                                              method = if (exact) "exact" else "asymptotic")
}
for (n in 4:9) for (r in 1:3) add_case(sprintf("exact_n%d_%d", n, r), sample(n), sample(n), TRUE)
add_case("exact_perfect_9", 1:9, 1:9 * 2, TRUE)
add_case("ties_12", round(rnorm(12), 0), round(rnorm(12), 0), FALSE)
add_case("ties_20", sample(1:5, 20, TRUE), rnorm(20), FALSE)
for (r in 1:3) { x <- rnorm(30); add_case(sprintf("asymptotic_n30_%d", r), x, x * 0.4 * r + rnorm(30), FALSE) }
w(fmt(do.call(rbind, sp_in)), "spearman_inputs.tsv")
w(fmt(do.call(rbind, sp_out)), "spearman_results.tsv")

# ---- Fisher and Stouffer -----------------------------------------------------------------------------
cmb_in <- list(); cmb_out <- list()
cases <- list(few = c(0.04, 0.2, 0.5), tiny = c(1e-12, 0.3, 0.8, 0.02), one = c(0.037),
              many = runif(12), mixed = c(0.001, 0.999, 0.5, 0.05, 0.95))
for (nm in names(cases)) {
  p <- cases[[nm]]; wt <- seq_along(p) + 0.5
  cmb_in[[nm]] <- data.frame(case = nm, p = p, weight = wt)
  x <- -2 * sum(log(p))
  z <- sum(qnorm(1 - p)) / sqrt(length(p)); zw <- sum(wt * qnorm(1 - p)) / sqrt(sum(wt^2))
  cmb_out[[nm]] <- data.frame(case = nm, fisher_x = x, fisher_p = pchisq(x, 2 * length(p), lower.tail = FALSE),
                              stouffer_z = z, stouffer_p = pnorm(z, lower.tail = FALSE),
                              weighted_z = zw, weighted_p = pnorm(zw, lower.tail = FALSE))
}
w(fmt(do.call(rbind, cmb_in)), "combination_inputs.tsv")
w(fmt(do.call(rbind, cmb_out)), "combination_results.tsv")

# ---- provenance --------------------------------------------------------------------------------------
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("nlme", as.character(packageVersion("nlme"))),
  "seed 20260924",
  "script make_glm_lmm_fixtures.R"
), file.path(out, "PROVENANCE_glm_lmm.txt"))
