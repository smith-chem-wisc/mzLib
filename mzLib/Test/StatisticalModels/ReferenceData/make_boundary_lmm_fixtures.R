# Reference values for two real proteins whose REML optimum is interior but lies near the boundary where one variance
# component is zero: NestedMixedModel's simplex used to stop on that boundary (a worse criterion) for them. Data: the
# MSstatsTMT controlled mixtures, SPS-MS3 (MassIVE MSV000084264; Huang et al. 2020, MCP 19:1706), searched by MetaMorpheus
# and converted to MSstatsTMT's PSM table (QuantProject results/reqqp15/p2), prepared as QuantProject's model A reference
# does: log2, each (run, channel) median moved to the median of those medians, pooled-reference channels dropped.
#
#   Rscript make_boundary_lmm_fixtures.R <msstatstmt_input.csv>
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# Model per protein: y ~ ion + Condition + (1 | individual / sample), REML (lme4::lmer), with
#   ion = peptide ion (sequence and charge) in its run; individual = mixture:channel; sample = run:channel.
suppressPackageStartupMessages({ library(lme4) })
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

values <- NULL; results <- NULL
for (p in c("O75569", "Q9H479")) {
  x <- d[d$ProteinName == p, ]
  x$ion <- factor(x$ion, levels = sort(unique(x$ion), method = "radix"))
  m <- lmer(y ~ ion + Condition + (1 | individual / sample), data = x, REML = TRUE, control = ctrl)
  vc <- as.data.frame(VarCorr(m))
  b <- fixef(m); v <- as.matrix(vcov(m))
  cn <- c("Condition0.5", "Condition0.667", "Condition1")
  values <- rbind(values, data.frame(protein = p, ion = as.character(x$ion), sample = x$sample, individual = x$individual,
                                     condition = as.character(x$Condition), y = fmt(x$y)))
  results <- rbind(results, data.frame(protein = p, reml = fmt(REMLcrit(m)), sigma2 = fmt(sigma(m)^2),
    var_sample = fmt(vc$vcov[vc$grp == "sample:individual"]), var_individual = fmt(vc$vcov[vc$grp == "individual"]),
    coef_0.5 = fmt(b[cn[1]]), coef_0.667 = fmt(b[cn[2]]), coef_1 = fmt(b[cn[3]]),
    se_0.5 = fmt(sqrt(v[cn[1], cn[1]])), se_0.667 = fmt(sqrt(v[cn[2], cn[2]])), se_1 = fmt(sqrt(v[cn[3], cn[3]]))))
}
write.table(values, file.path(here, "lmm_boundary_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results, file.path(here, "lmm_boundary_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4"))),
  "script make_boundary_lmm_fixtures.R",
  "data: MSstatsTMT controlled mixtures SPS-MS3 (MassIVE MSV000084264), proteins O75569 and Q9H479, prepared as QuantProject's model A",
  "lmer REML, bobyqa rhoend 1e-12; lme4's default, bobyqa and Nelder_Mead optimizers all reach the same optimum for both proteins"
), file.path(here, "PROVENANCE_boundary_lmm.txt"))
cat("written\n")
