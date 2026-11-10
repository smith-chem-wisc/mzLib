# Reference values for NestedMixedModelFit.NuisanceCoefficients: each protein's peptide coefficients (treatment-coded
# against its first peptide in ordinal order) from the same lme4 fits as make_nested_lmm_fixtures.R, on that script's
# committed data (lmm_<scenario>_values.tsv, lmm_<scenario>_samples.tsv), with its formula and optimizer settings.
#
# Run ONCE, by hand, beside the data:
#   Rscript make_nuisance_coefficient_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
fmt <- function(v) ifelse(is.na(v), "NaN", sprintf("%.17g", v))
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(rhobeg = 0.2, rhoend = 1e-12, maxfun = 1e5),
                    check.conv.singular = "ignore", check.conv.grad = "ignore", calc.derivs = FALSE)

for (name in c("nested", "sample_only")) {
  samples <- read.delim(file.path(here, sprintf("lmm_%s_samples.tsv", name)), na.strings = c("NA", ""))
  values <- read.delim(file.path(here, sprintf("lmm_%s_values.tsv", name)), na.strings = c("NaN", "NA"))
  common <- colnames(samples)[which(colnames(samples) == "intercept"):ncol(samples)]
  random <- if (name == "nested") "(1 | individual) + (1 | individual:sample)" else "(1 | sample)"
  out <- NULL
  for (prot in unique(values$protein)) {
    d <- merge(values[values$protein == prot & !is.na(values$y), ], samples, by = "sample")
    d$peptide <- factor(d$peptide, levels = sort(unique(d$peptide), method = "radix"))
    if (nlevels(d$peptide) < 2) next
    terms <- paste(common[-1], collapse = " + ")
    f <- as.formula(sprintf("y ~ %s + peptide + %s", terms, random))
    fit <- tryCatch(lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl), error = function(e) NULL)
    if (is.null(fit)) next
    b <- fixef(fit)
    pep <- grep("^peptide", names(b))
    out <- rbind(out, data.frame(protein = prot, peptide = sub("^peptide", "", names(b)[pep]), estimate = fmt(b[pep])))
  }
  write.table(out, file.path(here, sprintf("lmm_%s_nuisance.tsv", name)), sep = "\t", quote = FALSE, row.names = FALSE)
}
writeLines(c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", R.version.string),
  paste("lme4", as.character(packageVersion("lme4")), "; lmerTest", as.character(packageVersion("lmerTest"))),
  "script make_nuisance_coefficient_fixtures.R (inputs: lmm_nested_*, lmm_sample_only_* values and samples)",
  "fits as make_nested_lmm_fixtures.R: y ~ common + peptide (treatment-coded, ordinal first level) + random, REML, bobyqa rhoend 1e-12"
), file.path(here, "PROVENANCE_nuisance_coefficients.txt"))
