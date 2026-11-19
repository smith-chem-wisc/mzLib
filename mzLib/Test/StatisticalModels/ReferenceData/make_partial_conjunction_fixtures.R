# Reference outputs for PValueCombination.FisherPartialConjunction and MultipleTesting.BenjaminiBogomolov.
#
# Run ONCE, by hand, and the outputs checked in beside this script. R is never a build, runtime or test
# dependency of mzLib: the tests read the frozen TSVs.
#
#   Rscript make_partial_conjunction_fixtures.R <out_dir>
#
# Partial conjunction (Benjamini and Heller 2008), Fisher's method: for n p-values and 1 <= u <= n, the statistic is
# X = -2 sum(log(p)) over the n - u + 1 largest p's, and the p is pchisq(X, 2 (n - u + 1), lower.tail = FALSE).
# No package implements it, so the closed form is computed here, with base R only.
#
# Benjamini-Bogomolov (2014): families are selected by Benjamini-Hochberg on one p per family
# (p.adjust(family_p, "BH") <= q); each selected family's p's are adjusted by p.adjust(..., "BH") within the family and
# compared at q * R / m. Families not selected are not adjusted (NaN).

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NaN")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "NaN", sprintf("%.17g", df[[j]]))
  df
}

set.seed(20261114)
# ---- partial conjunction -------------------------------------------------------------------------------
pc_cases <- list(
  two     = c(0.012, 0.31),
  three   = c(0.004, 0.03, 0.47),
  five    = c(0.20, 0.001, 0.049, 0.65, 0.0002),
  tiny    = c(1e-30, 1e-12, 3e-8),
  ties    = c(0.05, 0.05, 0.2, 0.05),
  random  = runif(7)
)
pc_in <- list(); pc_out <- list()
for (nm in names(pc_cases)) {
  p <- pc_cases[[nm]]
  for (i in seq_along(p)) pc_in[[length(pc_in) + 1]] <- data.frame(case = nm, p = p[i])
  n <- length(p)
  for (u in seq_len(n)) {
    kept <- sort(p, decreasing = TRUE)[seq_len(n - u + 1)]
    x <- -2 * sum(log(kept))
    pc_out[[length(pc_out) + 1]] <- data.frame(case = nm, u = u, n = n, statistic = x,
      pval = pchisq(x, 2 * (n - u + 1), lower.tail = FALSE))
  }
}
w(fmt(do.call(rbind, pc_in)), "partial_conjunction_inputs.tsv")
w(fmt(do.call(rbind, pc_out)), "partial_conjunction_results.tsv")

# ---- Benjamini-Bogomolov -------------------------------------------------------------------------------
q <- 0.05
family_p <- c(0.0004, 0.30, 0.011, 0.020, 0.75, 0.0001)
within <- list(c(0.001, 0.02, 0.04), c(0.10, 0.50), c(0.003, 0.07, 0.009, 0.20),
               c(0.0199), c(0.40, 0.90, 0.33), c(0.00001, 0.02, 0.6, 0.03, 0.011))
selected <- p.adjust(family_p, "BH") <= q
m <- length(family_p); R <- sum(selected)
bb_in <- list(); bb_out <- list()
for (f in seq_len(m)) {
  adj <- if (selected[f]) p.adjust(within[[f]], "BH") else rep(NA_real_, length(within[[f]]))
  for (i in seq_along(within[[f]])) {
    bb_in[[length(bb_in) + 1]] <- data.frame(family = f, family_p = family_p[f], p = within[[f]][i])
    bb_out[[length(bb_out) + 1]] <- data.frame(family = f, index = i, selected = as.integer(selected[f]),
      adjusted = adj[i], level = q * R / m)
  }
}
w(fmt(do.call(rbind, bb_in)), "bogomolov_inputs.tsv")
w(fmt(do.call(rbind, bb_out)), "bogomolov_results.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  "packages none (base R: pchisq, p.adjust)",
  "seed 20261114",
  "script make_partial_conjunction_fixtures.R"
), file.path(out, "PROVENANCE_partial_conjunction.txt"))
