# Activity Log

This file tracks progress across agent sessions. Each session should append an entry below.

---

## 2026-02-19 — Session 0 (Architect)

**What was done:**
- Read and analyzed the MSstatsTMT paper (main text + supplemental) for the SpikeIn-5mix-MS3 dataset
- Explored the full Quantification project codebase (interfaces, strategies, engine, tests)
- Explored FlashLFQ normalization and median polish algorithms (future porting targets)
- Created `TMT_Spike-In_Info/TMT_SpikeIn_Reference.md` with complete experimental design, file mappings, expected fold changes, and data format documentation
- Designed 8-task implementation plan with detailed instruction documents in `.claude/plans/tasks/`
- Created `plan.md` (this project's orchestration document) and `Activity.md` (this file)

**Task instruction files created:**
- `.claude/plans/tasks/01-parse-tmt-columns.md`
- `.claude/plans/tasks/02-psm-tsv-adapter.md`
- `.claude/plans/tasks/03-experimental-design.md`
- `.claude/plans/tasks/04-pivot-by-file.md`
- `.claude/plans/tasks/05-combine-peptide-matrices.md`
- `.claude/plans/tasks/06-run-tmt-pipeline.md`
- `.claude/plans/tasks/07-synthetic-tmt-tests.md`
- `.claude/plans/tasks/08-spike-in-test-harness.md`

**Key decisions made with user:**
- TMT engine path is SEPARATE from LFQ path (do not modify existing Pivot/Run)
- Extend existing SpectrumMatchFromTsv reader + adapter in tests (both approach)
- ExperimentalDesign is hardcoded in test code (not a config file parser)
- Start with NoNormalization; global median normalization comes later
- SILAC labels: K8 (+8.014168 Da), R10 (+10.008269 Da) on the HeLa background

**Known issues:**
- Current .psmtsv search results have all-zero TMT reporter ion intensities (data needs re-searching)
- Mixture4_01 raw file is missing from the dataset (only 14 of 15 expected files)

**Next session should:** Start with Task 1 (parse TMT columns), Task 3 (experimental design), or Task 4 (PivotByFile) — these have no dependencies.

---

## 2026-02-20 — Session 1 (Implementation)

**What was done:**
- Implemented Tasks 1, 3, 4, 5, 6 in a single session. All build and pass existing tests.

**Task 1 — Parse TMT reporter ion columns (SpectrumMatchFromTsv):**
- Added `TmtChannelNames` string array to `SpectrumMatchFromTsvHeader.cs` (TMT10 + TMT16+ channels)
- Added TMT channel names to `ParseHeader()` in `SpectrumMatchTsvReader.cs` using `parsedHeader[name] = Array.IndexOf(spl, name)`
- Added `QuantValues` property (`double[]?`) to `SpectrumMatchFromTsv.cs`
- Added `ParseReporterIonColumns()` static helper that returns `null` if no TMT columns found
- Called it at the end of the main constructor: `QuantValues = ParseReporterIonColumns(spl, parsedHeader);`
- Fixed pre-existing build error: removed non-existent `using Quantification.Strategies.Collapse;` from `QuantificationParameters.cs`

**Task 3 — SpikeInExperimentalDesign:**
- Created `Test/Quantification/TestHelpers/SpikeInExperimentalDesign.cs`
- Implements `IExperimentalDesign` for the 14-file TMT10 SpikeIn-5mix-MS3 dataset
- Encodes all channel-to-condition mappings for Mixtures 1–5 (Mixture4_01 is absent from dataset)
- Dictionary keys are file names WITHOUT extension (no `.raw`), to match `FileNameWithoutExtension` from the reader

**Task 4 — PivotByFile():**
- Added to `QuantificationEngine.cs` as a `public static` method
- Groups spectral matches by `FullFilePath`, creates one dense `SpectralMatchMatrix` per file
- Uses direct positional copy of `QuantValues` into the matrix rows

**Task 5 — CombinePeptideMatrices():**
- Added to `QuantificationEngine.cs` as a `public static` method
- Collects union of all peptides (ordered by `FullSequence`) and all columns (ordered by file path then channel)
- Copies per-file peptide values into correct column offsets in the combined matrix

**Task 6 — RunTmt() and RunTmtAndReturnProteinMatrix():**
- Both delegate to a private `RunTmtCore(out proteinMatrixNorm)` to avoid code duplication
- Pipeline: PivotByFile → per-file PSM normalize → per-file PSM-to-peptide roll-up → CombinePeptideMatrices → peptide normalize → collapse → protein roll-up → protein normalize

**Commit:** `b22e2d05` — "Implement Tasks 1, 3, 4, 5, 6: TMT quantification pipeline foundation"

**Build/Test status:** Build succeeds (0 errors), 14/14 quantification tests pass.

**Remaining tasks:**
- Task 2 (PSM adapter: `PsmFromTsv → ISpectralMatch`) — needed before Task 7 and 8
- Task 7 (Synthetic data tests for TMT pipeline)
- Task 8 (Real spike-in data test harness)

**Next session should:** Implement Task 2 (PSM adapter), then Task 7 (synthetic TMT tests). Task 8 depends on re-searched psmtsv files with real TMT reporter ion intensities (data was all-zero in current searches).

---

## 2026-02-20 — Session 2 (Implementation)

**What was done:**
- Implemented Tasks 2 and 7. All build and pass (18/18 quantification tests pass).

**Task 2 — PsmTsvQuantAdapter (Test/Quantification/TestHelpers/PsmTsvQuantAdapter.cs):**
- `LoadSpectralMatches()`: reads `.psmtsv` via `SpectrumMatchTsvReader.ReadPsmTsv`, filters by q-value and target/decoy flag, wraps each record in `BaseSpectralMatch` with `QuantValues` copied from the parsed TMT reporter ion columns (Task 1 work)
- `GetUniqueAccessions()`: collects unique Accession strings from a psmtsv file
- Key design choice: `FullFilePath = FileNameWithoutExtension` (no extension) to match `SpikeInExperimentalDesign` dict keys (which also have no extension)
- Fixed namespace ambiguity: `Readers.ISpectralMatch` and `Omics.ISpectralMatch` both exist; used `using ISpectralMatch = Omics.ISpectralMatch;` alias and fully-qualified `Readers.SpectrumMatchTsvReader`

**Task 7 — TmtSpikeInTests (Test/Quantification/TmtSpikeInTests.cs):**
- `PivotByFile_CreatesCorrectPerFileMatrices`: creates 2 files × 2 PSMs each, verifies per-file matrix dimensions and exact intensity values; PSMs ordered by FullSequence within each file
- `CombinePeptideMatrices_MergesCorrectly`: manually constructs per-file PeptideMatrix objects, verifies union of 3 peptides in 6 columns (3 per file), correct values and zeros for missing peptides
- `RunTmt_FullPipeline_ProducesCorrectProteinMatrix`: end-to-end test with 2 proteins, 2 files, 3 channels, known QuantValues; verifies protein matrix has 2 rows and 6 columns with correct sums
- `RunTmt_WithSumCollapse_CombinesTechnicalReplicates`: verifies SumCollapse reduces 3 files × 3 channels (9 cols) to 3 collapsed condition columns with summed values
- Defined `TmtTestExperimentalDesign` private inner class (mirrors `TestExperimentalDesign` pattern from QuantificationTests.cs)

**Commit:** `3872147c` — "Implement Tasks 2 and 7: PsmTsvQuantAdapter and synthetic TMT tests"

**Build/Test status:** Build succeeds (0 errors), 18/18 quantification tests pass (14 existing + 4 new TmtSpikeInTests).

**Remaining tasks:**
- Task 8 (Real spike-in data test harness) — depends on re-searched psmtsv files with real TMT reporter ion intensities

**Next session should:** Implement Task 8 once re-searched psmtsv files (with actual TMT reporter ion values) are available. The PsmTsvQuantAdapter from Task 2 is ready to load those files.

---

## 2026-02-20 — Session 3 (Implementation)

**What was done:**
- Implemented Task 8. All 20 quantification tests pass (6 TmtSpikeIn + 14 existing).

**Bug Fix — `131N` channel name:**
- The actual psmtsv files use `131N` as the last TMT10 channel name, but `TmtChannelNames` only had `"131"` → `ParseReporterIonColumns` was producing 9 values instead of 10, causing `SetRow` dimension mismatch
- Added `"131N"` to `TmtChannelNames` after `"131"` in `SpectrumMatchFromTsvHeader.cs`
- Updated `SpikeInExperimentalDesign.cs` to use `"131N"` as the last channel label

**`PsmTsvQuantAdapter.BuildQuantificationInputs()` (new method):**
- Groups PSMs by `BaseSeq` → creates one `Protein(BaseSeq, accession)` per unique base sequence
- Digests each protein with `maxMissedCleavages: 100` to recover the full peptide even if it has internal K/R
- Calls `AddIdentifiedBioPolymer(peptide)` on each `BaseSpectralMatch` so `GetPsmToPeptideMap` can map PSMs to peptides
- Groups by first-non-decoy accession → creates one `BioPolymerGroup` per protein
- Filters: target PSMs only, q-value ≤ 0.01, non-zero QuantValues

**`QuantificationEvaluator.cs` (new helper in `TestHelpers/`):**
- `GetMeanIntensityByCondition()`: groups columns by `Condition` (skipping reference channels), averages non-zero intensities
- `CalculateFoldChange()`: numerator/denominator with NaN guards
- `ComputeFoldChanges()`: loops over all proteins, returns list of (accession, foldChange) pairs
- `Median()`: returns median of a sequence

**Task 8 Tests (added to `TmtSpikeInTests.cs`):**
- `LoadAndRunSpikeInData_BasicPipeline`: loads UPS PSMs, runs pipeline, asserts 140 columns + > 0 proteins (38s runtime due to 401k lines in psmtsv)
- `EvaluateFoldChanges_UPSProteins`: checks median fold change "1" vs "0.125" > 1.5 (true value is 8.0; passes with 39s runtime)
- `EvaluateFoldChanges_HelaBackground`: `[Explicit]` — processes HeLa file (~1.35M lines), checks median fold change "1" vs "0.5" is in [0.3, 3.0]

**Path resolution fix:**
- `TestContext.CurrentContext.TestDirectory` resolves to `Test/bin/Debug/net8.0-windows/` (framework subfolder exists even with `AppendTargetFrameworkToOutputPath=false`)
- Used directory-walking approach (`GetSolutionDir()` walks up until it finds `TMT_Spike-In_Info/`) — robust to any output depth

**Commit:** `852aab68` — "Implement Task 8: Real spike-in data test harness"

**Build/Test status:** Build succeeds (0 errors), 20/20 quantification tests pass. 6 TmtSpikeIn tests pass in ~90s total (2 real-data tests each take ~38s).

**All tasks DONE.** The TMT quantification testing harness is complete:
- Tasks 1–6: TMT pipeline foundation (reader, adapter, engine methods)
- Task 7: Synthetic TMT pipeline tests (fast)
- Task 8: Real SpikeIn-5mix-MS3 integration tests

---

## 2026-02-20 — Session 4 (Strategy Optimization)

**What was done:**
- Implemented Tasks 9–14 from plan2.md in a single session. All 26 quantification tests pass.

**Task 10 — GlobalMedianNormalization (Quantification/Strategies/Normalization/GlobalMedianNormalization.cs):**
- Equalizes median log2 intensity across all columns by shifting each column so its median matches the global median
- Zero values (missing data) excluded from median computation and preserved unchanged
- Columns with ≤1 positive value are copied as-is (no meaningful shift)

**Task 11 — ReferenceChannelNormalization (Quantification/Strategies/Normalization/ReferenceChannelNormalization.cs):**
- Divides each channel by the mean of the reference channels (126, 131N) within the same file
- Reference channels become 1.0; non-reference channels become fold-change-relative-to-reference
- If both reference channels are zero for a row, all channels for that row/file are zeroed (missing)

**Task 12 — MedianRollUp (Quantification/Strategies/RollUp/MedianRollUp.cs):**
- Rolls up PSMs→peptides or peptides→proteins using per-column median (instead of sum)
- Zero values excluded from median; result is 0 if all values are zero
- More robust to outlier PSMs/peptides than SumRollUp

**Task 13 — MeanCollapse (Quantification/Strategies/Collapse/MeanCollapse.cs):**
- Collapses technical replicates by averaging (sum / group count) instead of summing
- Same grouping logic as SumCollapse (Condition + '_' + BiologicalReplicate)
- Denominator is total group size (treats zero as "observed zero", not missing)

**Task 9 — Baseline Metrics + RunPipelineWithStrategies helper:**
- Added MAE, MAE(log2), FractionWithinFactor metrics to QuantificationEvaluator.cs
- Added RunPipelineWithStrategies(path, ...) and RunPipelineWithStrategies(preloaded, ...) helpers to TmtSpikeInTests.cs
- Added BaselineMetrics_NoNormalization_UPS test (prints metrics, no hard assertions)
- Added BaselineMetrics_NoNormalization_HeLa test (Explicit, for HeLa background)

**Task 14 — StrategyEvaluation_AllCombinations_UPS:**
- Tests all 48 strategy combinations (3 psmNorm × 2 rollUp × 2 pepNorm × 2 collapse × 2 protNorm)
- PSMs loaded once and reused; each pipeline run completes in ~1–2s after loading
- Full evaluation completes in ~3 minutes
- Ranked results printed sorted by combined MAE(log2)

**Strategy Evaluation Results (UPS spike-in, SpikeIn-5mix-MS3):**
- **Best strategy**: No Normalization | Sum Roll-Up | No Normalization | Mean Collapse | Sum Roll-Up | No Normalization
  - MAE_log2(8x)=1.430, MAE_log2(2x)=0.367, MedFC(8x)=3.154, MedFC(2x)=1.594
- **Worst strategy**: combined MAE_log2 = 4.157 (protein-level GlobalMedianNormalization hurts accuracy)
- **Key finding**: ratio compression is severe — median fold change for 8x comparison is only ~3.15 (vs true 8.0)
  This is a known artifact of TMT-SPS-MS3 co-isolation interference
- **Normalization effects**: GlobalMedianNormalization at peptide/protein level hurts accuracy by equalizing
  out the real signal differences; ReferenceChannelNormalization at PSM level provides mild improvement
- **Collapse**: MeanCollapse slightly outperforms NoCollapse (rank 1 vs rank 2)

**Commit:** `374dcfd5` — "Implement Tasks 9-14: normalization/rollup/collapse strategies + evaluation"

**Build/Test status:** Build succeeds (0 errors), 26/26 quantification tests pass.
- 14 QuantificationTests (existing)
- 12 TmtSpikeInTests (6 existing + 6 new: unit tests for each new strategy + baseline + evaluation)

**All tasks in plan2.md are DONE.**
