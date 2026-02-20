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
