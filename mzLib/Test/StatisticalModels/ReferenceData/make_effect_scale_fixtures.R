# Reference fixtures for StatisticalModels.EffectScale (Test/StatisticalModels/EffectScaleTests.cs).
# R 4.6.1, mr.raps 0.2 (CRAN archive; Zhao et al. 2020, Ann. Stat. 48:1742). Run from this folder:
#   Rscript make_effect_scale_fixtures.R
# Writes effect_scale_inputs.tsv, effect_scale_results.tsv and PROVENANCE_effect_scale.txt. R is never run by the tests.
#
# Each scenario draws true effects theta, a reference set x = theta + e_x (se s_x) and a scaled set
# y = scale * theta + e_y (se sqrt(s_y^2 + tau^2)), the model mr.raps.overdispersed fits: y ~ N(beta x_true, s_y^2 + tau2).
suppressPackageStartupMessages(library(mr.raps))
set.seed(20261104)

scenarios <- list(
  list(name = "typical",               k = 300,  scale = 0.7,  tau = 0.10, theta_sd = 1.00, se_lo = 0.05, se_hi = 0.30),
  list(name = "no_overdispersion",     k = 100,  scale = 1.0,  tau = 0.00, theta_sd = 1.00, se_lo = 0.05, se_hi = 0.30),
  list(name = "strong_overdispersion", k = 60,   scale = 0.85, tau = 0.40, theta_sd = 1.00, se_lo = 0.05, se_hi = 0.30),
  list(name = "weak_effects",          k = 2000, scale = 0.8,  tau = 0.05, theta_sd = 0.25, se_lo = 0.10, se_hi = 0.40),
  list(name = "negative_scale",        k = 80,   scale = -0.6, tau = 0.15, theta_sd = 1.00, se_lo = 0.05, se_hi = 0.30),
  list(name = "few_features",          k = 5,    scale = 0.9,  tau = 0.10, theta_sd = 1.00, se_lo = 0.05, se_hi = 0.30)
)

g <- function(v) sprintf("%.17g", v)
inputs <- data.frame()
results <- data.frame()
for (s in scenarios) {
  theta <- rnorm(s$k, 0, s$theta_sd)
  sx <- runif(s$k, s$se_lo, s$se_hi)
  sy <- runif(s$k, s$se_lo, s$se_hi)
  x <- theta + rnorm(s$k, 0, sx)
  y <- s$scale * theta + rnorm(s$k, 0, sqrt(sy^2 + s$tau^2))
  fit <- mr.raps.overdispersed(x, y, sx, sy, suppress.warning = TRUE)
  ivw <- sum(x * y / sy^2) / sum(x^2 / sy^2)
  bces <- sum(x * y) / (sum(x^2) - sum(sx^2))
  xw <- sum(x / sx^2) / sum(1 / sx^2)
  q <- sum((x - xw)^2 / sx^2)
  i2gx <- max(0, (q - (s$k - 1)) / q)
  inputs <- rbind(inputs, data.frame(scenario = s$name, reference = g(x), reference_se = g(sx), scaled = g(y), scaled_se = g(sy)))
  results <- rbind(results, data.frame(scenario = s$name, features = s$k, true_scale = g(s$scale), true_tau2 = g(s$tau^2),
    scale = g(fit$beta.hat), scale_se = g(fit$beta.se), tau2 = g(fit$tau2.hat), tau2_se = g(fit$tau2.se),
    ivw_scale = g(ivw), bces_scale = g(bces), i2_gx = g(i2gx)))
}
write.table(inputs, "effect_scale_inputs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results, "effect_scale_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("Generated (UTC):", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R:", R.version.string),
  paste("mr.raps:", as.character(packageVersion("mr.raps")), "(CRAN archive mr.raps_0.2.tar.gz; mr.raps.overdispersed, initialization 'simple', niter 20, tol sqrt(eps))"),
  "Seed: set.seed(20261104)",
  "Script: make_effect_scale_fixtures.R",
  "ivw_scale, bces_scale and i2_gx are the closed forms in the script, not from a package."
), "PROVENANCE_effect_scale.txt")
