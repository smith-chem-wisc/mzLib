# Reference values for two real proteins whose full nested model (1 | individual / sample) lme4 refuses, because the
# inner term (a channel within a plex run) holds one value per level: each protein has one ion in every plex run it is
# in. NestedMixedModel with estimableRandomTermsOnly leaves out what cannot be estimated (QuantProject STAT1 ruling 16):
#   - "individual": seen in 2+ runs of one mixture, so the individual (mixture:channel) is observed twice or more:
#     y ~ ion + Condition + (1 | individual), REML (lme4::lmer), Satterthwaite df (lmerTest);
#   - "none": seen at most once per individual, so neither random term can be estimated: y ~ ion + Condition (lm).
# Data: the MSstatsTMT controlled mixtures, SPS-MS3 (MassIVE MSV000084264; Huang et al. 2020, MCP 19:1706), searched by
# MetaMorpheus and converted to MSstatsTMT's PSM table (QuantProject results/reqqp15/p2), prepared as QuantProject's model
# A reference does: log2, each (run, channel) median moved to the median of those medians, pooled-reference channels
# dropped. Each kind's first protein by name with all four conditions and 30-120 values is taken.
#
#   Rscript make_estimable_terms_lmm_fixtures.R <msstatstmt_input.csv>
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)

d <- read.csv(args[1], stringsAsFactors = FALSE, na.strings = "NA")
d <- d[!is.na(d$Intensity) & d$Intensity > 0, ]
d$y <- log2(d$Intensity)
med <- tapply(d$y, paste(d$Run, d$Channel), median)
d$y <- d$y - med[paste(d$Run, d$Channel)] + median(med)
d <- d[d$Condition != "Norm" & d$Condition != "Empty", ]
d$ion <- paste(d$PeptideSequence, d$Charge, d$Run, sep = "@")
d$individual <- paste(d$Mixture, d$Channel, sep = ":")
d$sample <- paste(d$Run, d$Channel, sep = ":")
d$Condition <- factor(d$Condition, levels = c("0.125", "0.5", "0.667", "1"))

kind <- function(x) {
  if (length(unique(x$Condition)) < 4 || nrow(x) < 30 || nrow(x) > 120) return(NA)
  if (max(tapply(x$ion, x$Run, function(v) length(unique(v)))) != 1) return(NA)
  if (length(unique(x$sample)) != nrow(x)) return(NA)
  if (length(unique(x$individual)) < nrow(x)) "individual" else "none"
}
proteins <- sort(unique(d$ProteinName), method = "radix")
chosen <- c(individual = NA, none = NA)
for (p in proteins) {
  k <- kind(d[d$ProteinName == p, ])
  if (!is.na(k) && is.na(chosen[[k]])) chosen[[k]] <- p
  if (!any(is.na(chosen))) break
}
cat("chosen:", paste(names(chosen), chosen, sep = "=", collapse = ", "), "\n")

cn <- c("Condition0.5", "Condition0.667", "Condition1")
values <- NULL; results <- NULL; refused <- NULL
for (k in names(chosen)) {
  p <- chosen[[k]]
  x <- d[d$ProteinName == p, ]
  x$ion <- factor(x$ion, levels = sort(unique(x$ion), method = "radix"))
  full <- tryCatch({ lmer(y ~ ion + Condition + (1 | individual / sample), data = x, REML = TRUE, control = ctrl); "fitted" },
                   error = function(e) conditionMessage(e))
  refused <- c(refused, paste0(p, " (", k, "), full model: ", full))
  if (k == "individual") {
    m <- lmer(y ~ ion + Condition + (1 | individual), data = x, REML = TRUE, control = ctrl)
    vc <- as.data.frame(VarCorr(m))
    b <- fixef(m); v <- as.matrix(vcov(m))
    df <- sapply(cn, function(n) { L <- as.numeric(names(b) == n); contest(m, L, joint = FALSE)$df })
    reml <- REMLcrit(m); sigma2 <- sigma(m)^2; varInd <- vc$vcov[vc$grp == "individual"]
  } else {
    m <- lm(y ~ ion + Condition, data = x)
    b <- coef(m); v <- vcov(m)
    df <- setNames(rep(df.residual(m), 3), cn)
    reml <- NA; sigma2 <- summary(m)$sigma^2; varInd <- NA
  }
  values <- rbind(values, data.frame(protein = p, kind = k, ion = as.character(x$ion), sample = x$sample,
                                     individual = x$individual, condition = as.character(x$Condition), y = fmt(x$y)))
  results <- rbind(results, data.frame(protein = p, kind = k, n = nrow(x), reml = fmt(reml), sigma2 = fmt(sigma2),
    var_individual = fmt(varInd), coef_0.5 = fmt(b[cn[1]]), coef_0.667 = fmt(b[cn[2]]), coef_1 = fmt(b[cn[3]]),
    se_0.5 = fmt(sqrt(v[cn[1], cn[1]])), se_0.667 = fmt(sqrt(v[cn[2], cn[2]])), se_1 = fmt(sqrt(v[cn[3], cn[3]])),
    df_0.5 = fmt(df[[cn[1]]]), df_0.667 = fmt(df[[cn[2]]]), df_1 = fmt(df[[cn[3]]])))
}
write.table(values, file.path(here, "lmm_estimable_terms_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results, file.path(here, "lmm_estimable_terms_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4"))),
  paste("lmerTest", as.character(packageVersion("lmerTest"))),
  "script make_estimable_terms_lmm_fixtures.R",
  "data: MSstatsTMT controlled mixtures SPS-MS3 (MassIVE MSV000084264), prepared as QuantProject's model A",
  paste("chosen:", paste(names(chosen), chosen, sep = "=", collapse = ", ")),
  refused,
  "individual: lmer REML (1 | individual), bobyqa rhoend 1e-12, lmerTest Satterthwaite df; none: lm, residual df"
), file.path(here, "PROVENANCE_estimable_terms_lmm.txt"))
cat("written\n")
