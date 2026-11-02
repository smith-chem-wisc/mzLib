# Reference outputs for mzLib's LogitNormalMixedModel, from lme4 and glmmTMB.
#
# Run once, and the outputs checked in beside this script. R is never a runtime or test dependency of mzLib:
# the tests read the frozen TSVs.
#
#   Rscript make_lnmm_fixtures.R <out_dir>
#
# The model, per feature: the latent is eta_s = x_s' beta + u_group + e_s, with u ~ N(0, tau^2) and e ~ N(0, sigma^2).
# A measured sample reads eta_s. A counted sample has k of n trials succeed, each with chance g, where
# logit g = a + b eta_s + c d_s (the link; d_s the sample's depth).
#
# Three cases, each one R fits as a standard model:
#   * Every sample measured: lme4::lmer, REML and ML. With nothing counted the model is a linear mixed model.
#   * Every sample counted, sigma^2 given (tiny): lme4::glmer with 25-point adaptive Gauss-Hermite quadrature (AGQ).
#     The link enters as an offset a + c d and a rescaling by b: glmer's coefficients are b beta and its tau is b tau.
#     - ML is glmer's own fit. Its standard errors come from the Hessian of glmer's deviance in beta at the estimated
#       tau (numDeriv, Richardson), which is how mzLib defines them.
#     - REML integrates beta out of glmer's AGQ likelihood by a Laplace approximation, as glmmTMB's REML = TRUE does:
#       l(beta_hat) - log|H| / 2 + (p / 2) log 2 pi, H = -d2 l / d beta2 (numDeriv), maximized over tau here because
#       glmer has no REML for count data.
#     - glmer reports the log-likelihood relative to the saturated model; the full one is written here.
#     - At a large tau, AGQ with 25 points is itself inexact. Every fit's log-likelihood is checked against an exact
#       one (integrate(), per group), and a fit where they differ by more than 1e-6 is marked "reference_inexact".
#   * Every sample counted, sigma^2 estimated: glmmTMB with an observation-level random intercept, REML = TRUE and
#     ML. glmmTMB integrates the random effects by a Laplace approximation where mzLib uses quadrature, so the two
#     agree only to the Laplace error. The trials per sample are many, so that error is small. Under REML glmmTMB
#     reports beta as the joint mode of (beta, u), not the maximum of the integrated likelihood, so only its variances
#     and standard errors compare. The standard errors written are glmmTMB's at its estimated variances (REML: the beta
#     block of the inverse random-effect Hessian; ML: the inverse of the beta block of the fixed-effect Hessian), not
#     vcov(), which adds the variances' uncertainty.

suppressPackageStartupMessages({ library(lme4); library(glmmTMB); library(numDeriv); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20261011)

# One layout for every case: 8 animals (groups) with 2 to 4 runs each, a trait constant within an animal, and a
# covariate that varies between runs.
sizes <- c(2, 3, 2, 4, 3, 2, 3, 4)
group <- rep(paste0("a", seq_along(sizes)), sizes)
trait <- rep(c(0, 1, 1, 0, 1, 0, 0, 1), sizes)
x <- round(rnorm(length(group)), 2)
N <- length(group)
design <- data.frame(Intercept = 1, trait = trait, x = x)
w(fmt(design), "lnmm_design.tsv")
w(data.frame(group = group), "lnmm_groups.tsv")
na_row <- function() data.frame(b0 = NaN, b1 = NaN, b2 = NaN, se0 = NaN, se1 = NaN, se2 = NaN, sigma2 = NaN, tau2 = NaN, loglik = NaN)

# ---- every sample measured: lmer ------------------------------------------------------------------------------
F <- 60
Y <- matrix(NA_real_, F, N)
for (f in seq_len(F)) {
  tau <- c(0, 0.3, 1)[1 + f %% 3]
  u <- rnorm(length(sizes), 0, tau)[as.integer(factor(group, levels = unique(group)))]
  Y[f, ] <- -1 + rnorm(1, 0, 0.8) * trait + rnorm(1, 0, 0.3) * x + u + rnorm(N, 0, 0.6)
}
Y[matrix(runif(F * N) < 0.08, F, N)] <- NA
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore")
res <- list()
for (method in c("REML", "ML")) for (f in seq_len(F)) {
  obs <- !is.na(Y[f, ])
  d <- data.frame(y = Y[f, obs], trait = trait[obs], x = x[obs], group = group[obs])
  row <- cbind(data.frame(feature = f - 1, method = method, status = "fitted"), na_row())
  if (length(unique(d$group)) < 2 || max(table(d$group)) < 2) { row$status <- "not_fitted"; res[[length(res) + 1]] <- row; next }
  fit <- lmer(y ~ trait + x + (1 | group), data = d, REML = method == "REML", control = ctrl)
  row[, c("b0", "b1", "b2")] <- fixef(fit)
  row[, c("se0", "se1", "se2")] <- sqrt(diag(as.matrix(vcov(fit))))
  row$sigma2 <- sigma(fit)^2
  row$tau2 <- as.numeric(VarCorr(fit)$group[1])
  row$loglik <- as.numeric(logLik(fit))
  res[[length(res) + 1]] <- row
}
w(fmt(as.data.frame(Y)), "lnmm_measured_responses.tsv")
w(fmt(do.call(rbind, res)), "lnmm_measured_fit.tsv")

# ---- every sample counted, sigma^2 given: glmer, adaptive quadrature ------------------------------------------
# Two links: the identity, and one with an intercept, a slope below 1 and a depth term.
depth <- round(rnorm(N, 0, 0.7), 2)
links <- data.frame(name = c("identity", "link"), a = c(0, -0.4), b = c(1, 0.6), c = c(0, 0.3))
w(fmt(links), "lnmm_counted_links.tsv")
w(fmt(data.frame(depth = depth)), "lnmm_counted_depth.tsv")
gctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12, maxfun = 1e5), tolPwrss = 1e-12,
                      check.conv.singular = "ignore", check.conv.grad = "ignore", check.conv.hess = "ignore")
counted <- function(F, trials, tau_set, sigma, link, label) {
  K <- matrix(NA_integer_, F, N); M <- matrix(NA_integer_, F, N)
  for (f in seq_len(F)) {
    tau <- tau_set[1 + f %% length(tau_set)]
    u <- rnorm(length(sizes), 0, tau)[as.integer(factor(group, levels = unique(group)))]
    eta <- rnorm(1, -0.3, 0.5) + rnorm(1, 0, 0.8) * trait + rnorm(1, 0, 0.3) * x + u + rnorm(N, 0, sigma)
    M[f, ] <- trials()
    K[f, ] <- rbinom(N, M[f, ], plogis(link$a + link$b * eta + link$c * depth))
  }
  drop <- matrix(runif(F * N) < 0.08, F, N)
  K[drop] <- NA; M[drop] <- NA
  w(as.data.frame(K), sprintf("lnmm_%s_successes.tsv", label))
  w(as.data.frame(M), sprintf("lnmm_%s_trials.tsv", label))
  list(K = K, M = M)
}

# The exact log-likelihood at (beta, tau) on glmer's scale, by integrate() over each group's z = u / tau: in pieces
# around the integrand's mode, scaled by its value there, since integrate() over the whole line can miss a narrow peak.
exact_loglik <- function(d, beta, tau) {
  eta <- d$off + beta[1] + beta[2] * d$trait + beta[3] * d$x
  total <- 0
  for (g in unique(d$group)) {
    s <- which(d$group == g)
    lf <- function(z) sum(dbinom(d$k[s], d$n[s], plogis(eta[s] + tau * z), log = TRUE)) + dnorm(z, log = TRUE)
    if (tau == 0) { total <- total + lf(0) - dnorm(0, log = TRUE); next }
    mode <- optimize(lf, c(-40, 40), maximum = TRUE, tol = 1e-12)$maximum
    top <- lf(mode)
    h <- 1e-4
    width <- 1 / sqrt(max(-(lf(mode + h) - 2 * top + lf(mode - h)) / h^2, 1e-12))
    f <- function(z) sapply(z, function(v) exp(lf(v) - top))
    edges <- sort(c(mode + c(-1, 1) * (abs(mode) + 40), mode + c(-50, -5, 0, 5, 50) * min(width, 1)))
    piece <- sapply(seq_len(length(edges) - 1), function(i)
      integrate(f, edges[i], edges[i + 1], rel.tol = 1e-12, abs.tol = 1e-15 * min(width, 1), subdivisions = 2000L)$value)
    total <- total + top + log(sum(piece))
  }
  total
}

res <- list()
for (li in seq_len(nrow(links))) {
  link <- links[li, ]
  cnt <- counted(16, function() sample(1:8, N, replace = TRUE), c(0.4, 0.9, 1.5), 0, link, paste0("counted_", link$name))
  for (f in seq_len(nrow(cnt$K))) {
    obs <- !is.na(cnt$K[f, ])
    d <- data.frame(k = cnt$K[f, obs], n = cnt$M[f, obs], trait = trait[obs], x = x[obs], group = group[obs],
                    off = link$a + link$c * depth[obs])
    saturated <- sum(dbinom(d$k, d$n, d$k / d$n, log = TRUE))
    rows <- lapply(c("ML", "REML"), function(m) cbind(data.frame(link = link$name, feature = f - 1, method = m, status = "fitted"), na_row(), agq_error = NaN))
    form <- cbind(k, n - k) ~ trait + x + (1 | group) + offset(off)
    fit <- tryCatch(glmer(form, family = binomial, data = d, nAGQ = 25, control = gctrl), error = function(e) NULL)
    if (is.null(fit) || any(abs(fixef(fit) / link$b) > 30)) {
      for (r in rows) { r$status <- if (is.null(fit)) "r_error" else "separated"; res[[length(res) + 1]] <- r }
      next
    }
    dev <- glmer(form, family = binomial, data = d, nAGQ = 25, control = gctrl, devFunOnly = TRUE)
    # l(beta) = -dev / 2 + saturated at a given theta, beta maximized from a start; and H = -d2 l / d beta2 there.
    inner <- function(theta, start) {
      o <- nlminb(start, function(b) dev(c(theta, b)), control = list(rel.tol = 1e-15, x.tol = 1e-13, eval.max = 1e4, iter.max = 1e4))
      list(beta = o$par, loglik = -o$objective / 2 + saturated, h = numDeriv::hessian(function(b) dev(c(theta, b)), o$par) / 2)
    }
    p <- 3
    put <- function(row, theta, at, loglik) {
      row[, c("b0", "b1", "b2")] <- at$beta / link$b
      row[, c("se0", "se1", "se2")] <- sqrt(diag(solve(at$h))) / link$b
      row$tau2 <- (theta / link$b)^2
      row$loglik <- loglik
      exact <- exact_loglik(d, at$beta, theta)
      row$agq_error <- at$loglik - exact
      if (abs(exact - at$loglik) > 1e-6) row$status <- "reference_inexact"
      row
    }
    theta <- getME(fit, "theta")
    ml <- inner(theta, fixef(fit))
    res[[length(res) + 1]] <- put(rows[[1]], theta, ml, ml$loglik)
    # REML on mzLib's scale: l - log|b^2 H| / 2 + (p / 2) log 2 pi, maximized over log theta (theta = 0 if better).
    start <- fixef(fit)
    crit <- function(lt) {
      at <- inner(exp(lt), start)
      start <<- at$beta
      at$loglik - 0.5 * determinant(at$h * link$b^2)$modulus + 0.5 * p * log(2 * pi)
    }
    o <- optimize(function(lt) -crit(lt), c(log(link$b) - 6, log(link$b) + 2.5), tol = 1e-10)
    theta_reml <- exp(o$minimum)
    at0 <- inner(exp(-8) * link$b, start)
    crit0 <- at0$loglik - 0.5 * determinant(at0$h * link$b^2)$modulus + 0.5 * p * log(2 * pi)
    if (crit0 >= -o$objective) { theta_reml <- 0; at <- inner(0, start); value <- at$loglik - 0.5 * determinant(at$h * link$b^2)$modulus + 0.5 * p * log(2 * pi) }
    else { at <- inner(theta_reml, start); value <- -o$objective }
    res[[length(res) + 1]] <- put(rows[[2]], theta_reml, at, value)
  }
}
w(fmt(do.call(rbind, res)), "lnmm_counted_fit.tsv")

# ---- every sample counted, sigma^2 estimated: glmmTMB, Laplace ---------------------------------------------------
cnt <- counted(20, function() sample(20:40, N, replace = TRUE), c(0.4, 0.8), 0.4, links[1, ], "scatter")
res <- list()
for (method in c("REML", "ML")) for (f in seq_len(nrow(cnt$K))) {
  obs <- !is.na(cnt$K[f, ])
  d <- data.frame(k = cnt$K[f, obs], n = cnt$M[f, obs], trait = trait[obs], x = x[obs], group = group[obs],
                  obs = factor(seq_len(sum(obs))))
  row <- cbind(data.frame(feature = f - 1, method = method, status = "fitted"), na_row())
  fit <- tryCatch(glmmTMB(cbind(k, n - k) ~ trait + x + (1 | group) + (1 | obs), family = binomial, data = d,
                          REML = method == "REML",
                          control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-12))),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$sdr$pdHess)) { row$status <- "r_error"; res[[length(res) + 1]] <- row; next }
  vc <- VarCorr(fit)$cond
  row[, c("b0", "b1", "b2")] <- fixef(fit)$cond
  # Standard errors at the estimated variances, as mzLib defines them. vcov() adds the variances' own uncertainty (under
  # REML, by propagating it into beta as a random effect), which is not the same quantity.
  if (method == "REML") {
    env <- fit$obj$env
    v <- solve(env$spHess(env$last.par.best, random = TRUE))
    row[, c("se0", "se1", "se2")] <- sqrt(diag(as.matrix(v))[names(env$last.par.best)[env$random] == "beta"])
  } else {
    hf <- solve(fit$sdr$cov.fixed)
    ib <- which(rownames(hf) == "beta")
    row[, c("se0", "se1", "se2")] <- sqrt(diag(solve(hf[ib, ib])))
  }
  row$sigma2 <- as.numeric(vc$obs[1])
  row$tau2 <- as.numeric(vc$group[1])
  row$loglik <- as.numeric(logLik(fit))
  res[[length(res) + 1]] <- row
}
w(fmt(do.call(rbind, res)), "lnmm_scatter_fit.tsv")

# ---- provenance --------------------------------------------------------------------------------------------------
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4"))),
  paste("glmmTMB", as.character(packageVersion("glmmTMB"))),
  paste("TMB", as.character(packageVersion("TMB"))),
  paste("Matrix", as.character(packageVersion("Matrix"))),
  paste("numDeriv", as.character(packageVersion("numDeriv"))),
  "seed 20261011",
  "script make_lnmm_fixtures.R"
), file.path(out, "PROVENANCE_lnmm.txt"))
