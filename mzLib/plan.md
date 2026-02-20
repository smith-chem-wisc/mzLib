# Quantification Project: TMT Testing Harness Build Plan

## Session Startup Routine (READ THIS FIRST)

1. Run `pwd` to confirm you are in `C:/Users/Alex/Source/Repos/mzLib/mzLib/`
2. Read `Activity.md` in this directory to see what previous sessions accomplished
3. Run `git log --oneline -10` to see recent commits
4. Check the task list below - find the first task with status `TODO` whose dependencies are all `DONE`
5. Read the task's instruction file in `tasks/`
6. Work on that task. Build and test after each meaningful change.
7. Commit completed work with descriptive messages.

## End-of-Session Protocol (DO THIS BEFORE YOUR CONTEXT FILLS UP)

When you are nearing the end of your context window or finishing a task:
1. **Build**: `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"`
2. **Test**: `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"`
3. **Commit** any working changes: `git add <specific files> && git commit -m "descriptive message"`
4. **Update Activity.md**: Append a dated entry describing what you accomplished, what's working, and what remains
5. **Update this file**: Change the task status from `TODO` to `DONE` for any completed tasks. If a task is partially done, change to `IN PROGRESS` and note what remains in Activity.md
6. **Do NOT mark a task DONE unless**: it builds, tests pass, and you have verified the output

## Task List

Only change the **Status** field. Do not edit task descriptions or dependencies.

| # | Status | Instruction File | Description | Depends On |
|---|--------|-----------------|-------------|-----------|
| 1 | TODO | `tasks/01-parse-tmt-columns.md` | Parse TMT reporter ion columns in SpectrumMatchFromTsv | None |
| 2 | TODO | `tasks/02-psm-tsv-adapter.md` | Create adapter: PsmFromTsv → ISpectralMatch | Task 1 |
| 3 | TODO | `tasks/03-experimental-design.md` | Hardcoded SpikeIn ExperimentalDesign (14 files × 10 channels) | None |
| 4 | TODO | `tasks/04-pivot-by-file.md` | Add PivotByFile() to QuantificationEngine | None |
| 5 | TODO | `tasks/05-combine-peptide-matrices.md` | Add CombinePeptideMatrices() method | Task 4 |
| 6 | TODO | `tasks/06-run-tmt-pipeline.md` | Add RunTmt() pipeline method | Tasks 4, 5 |
| 7 | TODO | `tasks/07-synthetic-tmt-tests.md` | Synthetic data tests for TMT pipeline | Tasks 4, 5, 6 |
| 8 | TODO | `tasks/08-spike-in-test-harness.md` | Real spike-in data test harness | Tasks 1, 2, 3, 6 |

**Recommended execution order**: 1 → 3 → 4 → 5 → 6 → 2 → 7 → 8

## Architecture Rules (DO NOT VIOLATE)

- **Do NOT modify** existing `Pivot()`, `Run()`, or `RunAndReturnProteinMatrix()` methods in QuantificationEngine.cs
- **Add new methods alongside them**: `PivotByFile()`, `CombinePeptideMatrices()`, `RunTmt()`, `RunTmtAndReturnProteinMatrix()`
- **Use NoNormalization** for all normalization steps initially
- **Use SumRollUp** for all roll-up steps
- Follow existing code patterns from `Test/Quantification/QuantificationTests.cs`

## Key File Locations

```
Quantification/QuantificationEngine.cs        ← Tasks 4, 5, 6 ADD methods here
Readers/.../SpectrumMatchFromTsv.cs            ← Task 1 modifies
Readers/.../SpectrumMatchFromTsvHeader.cs       ← Task 1 modifies
Test/Quantification/QuantificationTests.cs     ← READ for patterns, do not modify
Test/Quantification/TestHelpers/               ← Tasks 2, 3 create files here
Test/Quantification/TmtSpikeInTests.cs         ← Tasks 7, 8 create/extend
TMT_Spike-In_Info/TMT_SpikeIn_Reference.md    ← Full experimental design reference
```

## Build and Test Commands

```bash
# Build
dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"

# Run quantification tests only
dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"

# Run specific test class
dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~TmtSpikeIn"
```

## Domain Context (for understanding the data)

This project quantifies TMT (Tandem Mass Tag) proteomics data. Key concepts:
- **TMT10-plex**: 10 isobaric labels (channels: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131) that allow multiplexing 10 samples in one MS run
- **Reporter ion intensities**: Stored in the last 10 columns of .psmtsv files, one per channel
- **SILAC HeLa background**: Heavy-labeled HeLa peptides (SILAC K8 = +8.014168 Da on K, SILAC R10 = +10.008269 Da on R) serve as constant background across all channels
- **UPS1 spike-in**: 48 universal standard proteins at known varying concentrations across channels
- **Reference channels**: Channels 126 and 131 contain pooled reference samples for normalization
- **Per-file processing**: TMT data must be normalized within each file before combining across files, because each PSM only has intensities for channels within its own file
