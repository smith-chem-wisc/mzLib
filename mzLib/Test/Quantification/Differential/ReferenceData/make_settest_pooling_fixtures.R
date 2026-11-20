# Reference values for SETTEST1 S3, pooling and the families (QuantProject DEF-DIFF-SET-POOL, -FAMILY, -MATURATION).
#
# Run ONCE, by hand, and the outputs checked in beside this script. R is never a build, runtime or test dependency of
# mzLib: SetTestPoolingTests reads the frozen TSVs.
#
#   Rscript make_settest_pooling_fixtures.R <out_dir>
#
# One scenario of seven hypothesis parts. Each dataset contrast is given as a set estimate (beta-bar, se, df), as S2's
# set model or a single's STAT1 row would give it; S2's tests check those. Here R computes, independently of mzLib:
#   * each dataset's p and 90% interval, by the rules of DEF-DIFF-SET-VERDICT (pt / qt);
#   * each pooled row: metafor rma(yi, sei = se, method = "DL", test = "adhoc") for the estimate, tau2 and the
#     Hartung-Knapp se under Knapp and Hartung's ad hoc correction, never below the DerSimonian-Laird se (ruled by
#     QuantProject's user 2026-11-19), then the same rules on t(k - 1); an attenuated_down pool's m-bar is the average
#     of the datasets' m-bar with the random-effects weights 1 / (se^2 + tau2);
#   * the aging family's and the maturation family's Benjamini-Hochberg adjusted p (p.adjust);
#   * R (aging primaries supported: adjusted p <= 0.05 and not contradicted; the scenario has no guard), the level
#     0.05 R / m, and each supported pooled part's datasets adjusted by p.adjust within the part;
#   * n_datasets_agree (datasets with p <= 0.05) and the "supported in at least 2" p: Fisher's combination of the
#     k - 1 largest dataset p's, pchisq(-2 sum(log p), 2 (k - 1), lower.tail = FALSE).

suppressPackageStartupMessages(library(metafor))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "."
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
fmt <- function(df) {
  for (j in seq_along(df)) if (is.numeric(df[[j]]))
    df[[j]] <- ifelse(is.na(df[[j]]), "", sprintf("%.17g", df[[j]]))
  df
}

alpha <- 0.05; delta <- 0.5; f <- 0.5
# part, prediction, dataset, kind (aging / maturation), independence, b, se, df, mrna (m-bar as S2 would give it)
d <- data.frame(stringsAsFactors = FALSE,
  part        = c("HA", "HA", "HA", "HB", "HB", "HC", "HD", "HD", "HE", "HE", "HF", "HG", "HG"),
  prediction  = c("up", "up", "up", "down", "down", "not_up", "unchanged", "unchanged", "attenuated_down",
                  "attenuated_down", "up", "up", "up"),
  dataset     = c("D1", "D2", "D3", "D1", "D4", "D2", "D3", "D5", "D1", "D4", "D6", "D2", "D3"),
  kind        = c(rep("aging", 10), "maturation", "aging", "aging"),
  independence = c(rep("independent", 11), "probably reproduces [12]", "independent"),
  b   = c(1.20, 1.10, 0.75, -0.60, -0.20, 0.10, 0.02, -0.03, -0.50, -0.90, 0.90, 0.70, 0.30),
  se  = c(0.15, 0.20, 0.35, 0.20, 0.30, 0.15, 0.08, 0.10, 0.20, 0.30, 0.30, 0.20, 0.20),
  df  = c(9.1, 12.4, 7.0, 10.0, 8.0, 11.0, 14.0, 9.5, 10.0, 8.0, 6.0, 10.0, 10.0),
  mrna = c(rep(NA, 8), -2.0, -1.6, rep(NA, 3)))

test <- function(prediction, b, se, df, mrna) {
  Tc <- function(t) pt(t, df)
  upper <- function(bound) 1 - Tc((b - bound) / se)
  lower <- function(bound) Tc((b - bound) / se)
  q <- qt(1 - alpha, df); low <- b - q * se; high <- b + q * se
  fm <- f * mrna
  switch(prediction,
    up = list(p = upper(0), contradicted = high < 0, low = low, high = high),
    down = list(p = lower(0), contradicted = low > 0, low = low, high = high),
    unchanged = list(p = max(upper(-delta), lower(delta)), contradicted = low > delta || high < -delta, low = low, high = high),
    not_up = list(p = lower(delta), contradicted = low > delta, low = low, high = high),
    attenuated_down = list(p = max(lower(0), upper(fm)), contradicted = low >= 0 || high < fm, low = low, high = high))
}

dt <- lapply(seq_len(nrow(d)), function(i) test(d$prediction[i], d$b[i], d$se[i], d$df[i], d$mrna[i]))
d$p <- sapply(dt, `[[`, "p"); d$ci_low <- sapply(dt, `[[`, "low"); d$ci_high <- sapply(dt, `[[`, "high")
d$contradicted <- sapply(dt, `[[`, "contradicted")
d$poolable <- d$kind == "aging" & d$independence == "independent"

parts <- unique(d$part)
pooled <- list(); primaries <- list()
for (pt_ in parts) {
  rows <- d[d$part == pt_ & d$poolable, ]
  k <- nrow(rows)
  if (k >= 2) {
    r <- rma(yi = rows$b, sei = rows$se, method = "DL", test = "adhoc")
    est <- as.numeric(r$beta)
    mbar <- NA
    if (rows$prediction[1] == "attenuated_down") {
      wt <- 1 / (rows$se^2 + r$tau2); mbar <- sum(wt * rows$mrna) / sum(wt)
    }
    tt <- test(rows$prediction[1], est, r$se, k - 1, mbar)
    kept <- sort(rows$p, decreasing = TRUE)[seq_len(k - 1)]
    pooled[[pt_]] <- data.frame(part = pt_, k = k, estimate = est, se = r$se, df = k - 1, tau2 = r$tau2, mrna_mean = mbar,
      p = tt$p, ci_low = tt$low, ci_high = tt$high, contradicted = as.integer(tt$contradicted),
      n_agree = sum(rows$p <= alpha), p_in_2 = pchisq(-2 * sum(log(kept)), 2 * (k - 1), lower.tail = FALSE))
    primaries[[length(primaries) + 1]] <- data.frame(part = pt_, dataset = "", family = "aging", p = tt$p,
      contradicted = tt$contradicted)
  } else if (k == 1) {
    primaries[[length(primaries) + 1]] <- data.frame(part = pt_, dataset = rows$dataset, family = "aging", p = rows$p,
      contradicted = rows$contradicted)
  }
  mat <- d[d$part == pt_ & d$kind == "maturation" & d$independence == "independent", ]
  for (i in seq_len(nrow(mat)))
    primaries[[length(primaries) + 1]] <- data.frame(part = pt_, dataset = mat$dataset[i], family = "maturation",
      p = mat$p[i], contradicted = mat$contradicted[i])
}
pooled <- do.call(rbind, pooled)
prim <- do.call(rbind, primaries)
prim$p_adjusted <- NA
for (fam in unique(prim$family)) prim$p_adjusted[prim$family == fam] <- p.adjust(prim$p[prim$family == fam], "BH")
prim$supported <- as.integer(prim$p_adjusted <= alpha & !prim$contradicted)
aging <- prim[prim$family == "aging", ]
R <- sum(aging$supported); m <- nrow(aging); level <- alpha * R / m

bb <- list()
for (pt_ in aging$part[aging$dataset == "" & aging$supported == 1]) {
  rows <- d[d$part == pt_ & d$poolable, ]
  bb[[pt_]] <- data.frame(part = pt_, dataset = rows$dataset, adjusted = p.adjust(rows$p, "BH"), level = level)
}
bb <- do.call(rbind, bb)

w(fmt(d[, c("part", "prediction", "dataset", "kind", "independence", "b", "se", "df", "mrna", "p", "ci_low", "ci_high")]),
  "settest_pooling_datasets.tsv")
w(fmt(pooled), "settest_pooling_pooled.tsv")
w(fmt(prim[, c("part", "dataset", "family", "p", "p_adjusted", "supported")]), "settest_pooling_family.tsv")
w(fmt(bb), "settest_pooling_bogomolov.tsv")
w(fmt(data.frame(aging_family_size = m, aging_supported = R, bb_level = level,
                 maturation_family_size = sum(prim$family == "maturation"))), "settest_pooling_summary.tsv")

writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("metafor", as.character(packageVersion("metafor"))),
  "seed none (fixed values)",
  "script make_settest_pooling_fixtures.R"
), file.path(out, "PROVENANCE_settest_pooling.txt"))
