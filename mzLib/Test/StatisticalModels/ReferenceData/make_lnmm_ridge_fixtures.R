# Exact reference values for mzLib's LogitNormalMixedModel where counts are one-sided or nearly so (mzLib #1511).
#
# Run once, and the outputs checked in beside this script. R is never a runtime or test dependency of mzLib:
# the tests read the frozen TSVs.
#
#   Rscript make_lnmm_ridge_fixtures.R <out_dir>
#
# The model (as make_lnmm_fixtures.R): per group, the latent of sample s is eta_s = x_s' beta + u + e_s, with
# u ~ N(0, tau^2) and e_s ~ N(0, sigma^2), and a counted sample has k_s of n_s trials succeed with logit = eta_s (the
# identity link). A group's likelihood is int N(u; 0, tau^2) prod_s int N(e; 0, sigma^2) Binom(k_s; n_s, eta_s) de du.
#
# lme4 and glmmTMB cannot serve here: glmer has no observation-level variance with quadrature, glmmTMB's Laplace error
# is large for 1-4 trials, and on a one-sided group at a large tau every Gauss-Hermite rule misses the step (the
# integrand is the prior times a smoothed step). So the integrals are computed directly:
#   * over e, by 200-node Gauss-Hermite on N(0, sigma^2) (the count probability is entire; its nearest poles are
#     pi / (sigma sqrt 2) from the real axis in GH's variable, so 200 nodes leave an error far below 1e-12);
#   * over u, by the trapezoid rule with step h = min(0.1, tau / 100) over the prior's +-10 tau widened to cover every
#     count's step with 40 to spare. The integrand is analytic in a strip about the real axis, so the trapezoid error
#     falls as exp(-2 pi d / h) and is negligible; the range truncation drops prior mass below e^-50.
#   A check: halving h changes no value by more than 1e-12 (reported as max_halving_change in the provenance).
#
# Outputs:
#   lnmm_ridge_groups.tsv  one group at a time, intercept-only: its log-likelihood over a grid of patterns, b0, tau, sigma.
#   lnmm_ridge_fits.tsv    #1511's layouts A and B (8 groups x 3 runs, 2 trials, group trait), sigma^2 = 0.25 given: the
#                          exact ML estimate (beta, tau^2, l), H = -d2 l / d beta2 there at fixed variances (numDeriv),
#                          and the ML fit with the trait coefficient fixed at 0.

suppressPackageStartupMessages({ library(numDeriv) })
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

gh <- function(n) { # Gauss-Hermite for exp(-x^2), Golub-Welsch, weights normalized to N(0, 1/2)
  i <- seq_len(n - 1); J <- matrix(0, n, n); J[cbind(i, i + 1)] <- sqrt(i / 2); J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, lw = log(e$vectors[1, ]^2))
}
E <- gh(200)

# log int N(e; 0, sigma^2) Binom(k; n, plogis(eta + e)) de, for a vector eta
log_count_direct <- function(k, n, eta, sigma) {
  lv <- outer(eta, sqrt(2) * sigma * E$x, "+")
  lv <- lchoose(n, k) + k * plogis(lv, log.p = TRUE) + (n - k) * plogis(-lv, log.p = TRUE) +
    matrix(E$lw, length(eta), length(E$x), byrow = TRUE)
  top <- lv[cbind(seq_len(nrow(lv)), max.col(lv, ties.method = "first"))]
  top + log(rowSums(exp(lv - top)))
}
# The same, tabulated once per (k, n, sigma) and read by a cubic spline: knots 0.005 apart within 60 of 0 and 0.05
# beyond, out to 1700. Its error is far below 1e-10 (checked against log_count_direct; max_spline_error in the provenance).
tables <- new.env()
log_count <- function(k, n, eta, sigma) {
  key <- paste(k, n, sigma)
  if (is.null(tables[[key]])) {
    knots <- sort(unique(c(seq(-1700, -60, by = 0.05), seq(-60, 60, by = 0.005), seq(60, 1700, by = 0.05))))
    tables[[key]] <- splinefun(knots, log_count_direct(k, n, knots, sigma), method = "fmm")
  }
  tables[[key]](eta)
}

# log int N(u; 0, tau^2) prod_s q_s(eta_s + u) du, eta_s = the samples' fixed parts
log_group <- function(k, n, eta, tau, sigma, halve = FALSE, direct = FALSE) {
  h <- min(0.1, tau / 100); if (halve) h <- h / 2
  steps <- -eta
  lo <- min(-10 * tau, min(steps) - 40); hi <- max(10 * tau, max(steps) + 40)
  u <- seq(lo, hi, by = h)
  lf <- dnorm(u, 0, tau, log = TRUE)
  q <- if (direct) log_count_direct else log_count
  for (s in seq_along(k)) lf <- lf + q(k[s], n[s], eta[s] + u, sigma)
  top <- max(lf)
  top + log(sum(exp(lf - top)) * h)
}

# ---- one group at a time ----------------------------------------------------------------------------------------
patterns <- list(c(0, 0, 0), c(2, 2, 2), c(1, 2, 2), c(0, 0, 1), c(0, 0, 0, 0), c(1, 1, 1, 0), c(4, 4, 4), c(0, 0, 0), c(4, 4, 3))
trials <- c(2, 2, 2, 2, 1, 1, 4, 4, 4)
rows <- list(); worst <- 0; worst_spline <- 0
for (pi in seq_along(patterns)) for (b0 in c(-12, -3, 0, 3, 12)) for (tau in c(0.05, 0.5, 2, 8, 15.25, 40, 148))
  for (sigma in c(0.1, 0.5, 1, 2, 3)) {
    k <- patterns[[pi]]; n <- rep(trials[pi], length(k)); eta <- rep(b0, length(k))
    v <- log_group(k, n, eta, tau, sigma)
    if (length(rows) %% 25 == 0) {
      worst <- max(worst, abs(v - log_group(k, n, eta, tau, sigma, halve = TRUE)))
      worst_spline <- max(worst_spline, abs(v - log_group(k, n, eta, tau, sigma, direct = TRUE)))
    }
    rows[[length(rows) + 1]] <- data.frame(successes = paste(k, collapse = ","), trials = trials[pi], b0 = b0, tau2 = tau^2,
                                           sigma2 = sigma^2, loglik = v)
  }
w(fmt(do.call(rbind, rows)), "lnmm_ridge_groups.tsv")

# ---- #1511's layouts --------------------------------------------------------------------------------------------
trait_g <- c(0, 0, 0, 0, 1, 1, 1, 1)
full_g <- c(1, 0, 0, 0, 1, 1, 1, 0)
layout <- function(which) {
  k <- lapply(1:8, function(g) rep(2 * full_g[g], 3))
  if (which == "A") k[[1]][1] <- 1 else k[[2]][1] <- 1
  k
}
loglik <- function(k, beta, tau, sigma) {
  total <- 0
  cache <- list()
  for (g in 1:8) {
    key <- paste(trait_g[g], paste(k[[g]], collapse = ","))
    if (is.null(cache[[key]])) cache[[key]] <- log_group(k[[g]], rep(2, 3), rep(beta[1] + beta[2] * trait_g[g], 3), tau, sigma)
    total <- total + cache[[key]]
  }
  total
}
sigma <- 0.5
fits <- list()
for (which in c("A", "B")) {
  k <- layout(which)
  # full ML: over (b0, b1, log tau), from a coarse profile's best point
  o <- optim(if (which == "A") c(-11.6, 20.4, log(11.9)) else c(-4.6, 13.7, log(12.4)), function(v) -loglik(k, v[1:2], exp(v[3]), sigma), method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 4000))
  o <- optim(o$par, function(v) -loglik(k, v[1:2], exp(v[3]), sigma), method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))
  beta <- o$par[1:2]; tau <- exp(o$par[3])
  H <- -hessian(function(b) loglik(k, b, tau, sigma), beta)
  # ML with the trait coefficient fixed at 0: over (b0, log tau)
  o0 <- optim(c(0, log(15)), function(v) -loglik(k, c(v[1], 0), exp(v[2]), sigma), method = "Nelder-Mead",
              control = list(reltol = 1e-14, maxit = 4000))
  o0 <- optim(o0$par, function(v) -loglik(k, c(v[1], 0), exp(v[2]), sigma), method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))
  fits[[which]] <- data.frame(layout = which, sigma2 = sigma^2, b0 = beta[1], b1 = beta[2], tau2 = tau^2, loglik = -o$value,
                              h00 = H[1, 1], h01 = H[1, 2], h11 = H[2, 2],
                              null_b0 = o0$par[1], null_tau2 = exp(2 * o0$par[2]), null_loglik = -o0$value)
}
w(fmt(do.call(rbind, fits)), "lnmm_ridge_fits.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("numDeriv", as.character(packageVersion("numDeriv"))),
  paste("max_halving_change", sprintf("%.3g", worst)),
  paste("max_spline_error", sprintf("%.3g", worst_spline)),
  "script make_lnmm_ridge_fixtures.R"
), file.path(out, "PROVENANCE_lnmm_ridge.txt"))
