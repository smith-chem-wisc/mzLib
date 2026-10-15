# Reference design matrices for Quantification.Differential.AnalysisDesign (STAT1 milestone M1).
#
# Run ONCE, by hand, to (re)generate the fixtures beside this script:
#   Rscript make_design_fixtures.R
# R is never a build, runtime or test dependency: the C# tests read the committed TSVs only.
#
# For each design this writes
#   design_<name>_samples.tsv  one row per sample, ordered by sample_id (ordinal), the columns the C# test reads
#   design_<name>_matrix.tsv   stats::model.matrix for that design, same row order, R's column names as header
#
# Level order matches AnalysisDesign: the reference level first, then the other levels in ordinal (C-locale)
# order; a factor with no declared reference is all levels in ordinal order. Batch levels are ordinal.
# Terms appear in the order AnalysisDesign builds them: intercept, condition factors (spec order), batch,
# covariates, interactions.

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
here <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()

lev <- function(x, ref = NULL) {
  others <- sort(unique(x[x != if (is.null(ref)) "" else ref]), method = "radix")
  factor(x, levels = c(ref, others))
}

fmt <- function(v) if (is.numeric(v)) sprintf("%.17g", v) else v
# model.matrix needs the factors typed; build the typed frame from the written columns.
mm_frame <- function(s, refs = list(), batch = FALSE) {
  d <- s
  for (col in names(refs)) d[[col]] <- lev(s[[col]], refs[[col]])
  if (batch) d$batch <- lev(as.character(s$batch))
  d
}

go <- function(name, s, refs, formula, batch = FALSE) {
  s <- s[order(s$sample_id, method = "radix"), , drop = FALSE]
  rownames(s) <- NULL
  d <- mm_frame(s, refs, batch)
  mm <- model.matrix(formula, data = d)
  write.table(s, file.path(here, sprintf("design_%s_samples.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  m <- as.data.frame(lapply(as.data.frame(mm, check.names = FALSE), fmt), check.names = FALSE)
  write.table(m, file.path(here, sprintf("design_%s_matrix.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(name, ": ", nrow(mm), " x ", ncol(mm), ", rank ", qr(mm)$rank, "\n", sep = "")
}

ids <- function(n) sprintf("s%02d", seq_len(n))

# 1. balanced: two groups of three
go("balanced",
   data.frame(sample_id = ids(6), condition = rep(c("A", "B"), each = 3)),
   list(condition = "A"), ~ condition)

# 2. unbalanced: two against four
go("unbalanced",
   data.frame(sample_id = ids(6), condition = c("A", "A", "B", "B", "B", "B")),
   list(condition = "A"), ~ condition)

# 3. a lost sample: three levels, one of them down to two
go("lost_sample",
   data.frame(sample_id = ids(8), condition = c("ctrl", "ctrl", "ctrl", "x", "x", "x", "y", "y")),
   list(condition = "ctrl"), ~ condition)

# 4. paired: four individuals before and after; individual is a random effect, so not in the fixed matrix
go("paired",
   data.frame(sample_id = ids(8), time = rep(c("pre", "post"), 4),
              individual = rep(c("m1", "m2", "m3", "m4"), each = 2)),
   list(time = "pre"), ~ time)

# 5. within and between individuals: genotype differs between animals, time within each animal
go("within_between",
   data.frame(sample_id = ids(12), genotype = rep(c("wt", "ko"), each = 6), time = rep(c("d0", "d7"), 6),
              individual = rep(sprintf("m%d", 1:6), each = 2)),
   list(genotype = "wt", time = "d0"), ~ genotype + time)

# 6. batch: two conditions spread over two batches
go("batch",
   data.frame(sample_id = ids(8), condition = rep(c("A", "B"), 4), batch = rep(c("1", "2"), each = 4)),
   list(condition = "A"), ~ condition + batch, batch = TRUE)

# 7. age slope, centre 50 years, per decade
s7 <- data.frame(sample_id = ids(8), age = c(22, 31, 38, 47, 55, 63, 71, 84))
go("age_slope", s7, list(), ~ I((age - 50) / 10))

# 8. age with no stated centre: centred at the stratum mean, plus sex
s8 <- data.frame(sample_id = ids(8), sex = rep(c("female", "male"), 4), age = c(3, 6, 9, 12, 15, 18, 24, 27))
s8 <- s8[order(s8$sample_id, method = "radix"), ]
m8 <- mean(s8$age)
go("age_default_centre", s8, list(sex = "female"), as.formula(sprintf("~ sex + I((age - %.17g) / 10)", m8)))

# 9. two factors, additive
s9 <- data.frame(sample_id = ids(12), age_group = rep(c("young", "old"), each = 6),
                 treatment = rep(c("vehicle", "drug"), 6))
go("two_factors", s9, list(age_group = "young", treatment = "vehicle"), ~ age_group + treatment)

# 10. two factors with their interaction
go("interaction", s9, list(age_group = "young", treatment = "vehicle"), ~ age_group + treatment + age_group:treatment)

# 11. a reference that is not first alphabetically
go("nonalpha_reference",
   data.frame(sample_id = ids(9), age_group = rep(c("old", "middle", "young"), 3)),
   list(age_group = "young"), ~ age_group)

# 12. confounded: every A is batch 1 and every B batch 2
go("confounded",
   data.frame(sample_id = ids(6), condition = rep(c("A", "B"), each = 3), batch = rep(c("1", "2"), each = 3)),
   list(condition = "A"), ~ condition + batch, batch = TRUE)

# 13. no reference declared: all levels in ordinal order, the first is the coding baseline
go("no_reference",
   data.frame(sample_id = ids(9), condition = rep(c("C", "A", "B"), 3)),
   list(condition = NULL), ~ condition)

# 14. factor, batch and a covariate together
go("factor_batch_covariate",
   data.frame(sample_id = ids(10), condition = rep(c("ctrl", "treated"), 5), batch = rep(c("b1", "b2"), each = 5),
              age = c(40, 45, 52, 58, 61, 39, 47, 50, 66, 70)),
   list(condition = "ctrl"), ~ condition + batch + I((age - 50) / 10), batch = TRUE)

writeLines(c(
  sprintf("generated_utc %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R %s", R.version.string),
  "script make_design_fixtures.R (base R stats::model.matrix only; no packages)",
  sprintf("age_default_centre mean age %.17g", m8)
), file.path(here, "PROVENANCE_design.txt"))
