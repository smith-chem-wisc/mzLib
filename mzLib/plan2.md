# Quantification Project: Strategy Optimization for TMT Accuracy

## Goal

Iteratively implement and evaluate normalization, roll-up, and collapse strategies to find the combination that most accurately recovers known fold changes from the SpikeIn-5mix-MS3 dataset.

**Success criteria:**
- UPS proteins: median fold change for "1" vs "0.125" approaches 8.0 (currently unknown with NoNormalization baseline)
- UPS proteins: median fold change for "1" vs "0.5" approaches 2.0
- HeLa background proteins: median fold change between any two conditions approaches 1.0
- Low coefficient of variation (CV) across technical replicates

## Session Startup Routine (READ THIS FIRST)

1. Run `pwd` to confirm you are in `C:/Users/Alex/Source/Repos/mzLib/mzLib/`
2. Read `Activity.md` to see what previous sessions accomplished
3. Run `git log --oneline -10` to see recent commits
4. Check the task list below — find the first task with status `TODO` whose dependencies are all `DONE`
5. Read the task's instruction file in `tasks/`
6. Work on that task. Build and test after each meaningful change.
7. Commit completed work with descriptive messages.

## End-of-Session Protocol (DO THIS BEFORE YOUR CONTEXT FILLS UP)

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
| 1 | DONE | `tasks/09-baseline-metrics.md` | Establish baseline accuracy metrics with current NoNormalization pipeline | None |
| 2 | DONE | `tasks/10-global-median-normalization.md` | Implement GlobalMedianNormalization strategy | None |
| 3 | DONE | `tasks/11-reference-channel-normalization.md` | Implement ReferenceChannelNormalization strategy | None |
| 4 | DONE | `tasks/12-median-rollup.md` | Implement MedianRollUp strategy | None |
| 5 | DONE | `tasks/13-mean-collapse.md` | Implement MeanCollapse strategy | None |
| 6 | DONE | `tasks/14-strategy-evaluation.md` | Combinatorial evaluation of all strategies against spike-in data | Tasks 1–5 |

**Recommended execution order**: 1 → 2 → 3 → 4 → 5 → 6

## Architecture Rules (DO NOT VIOLATE)

- **Do NOT modify** existing strategies (NoNormalization, SumRollUp, NoCollapse, SumCollapse)
- **Do NOT modify** existing tests in TmtSpikeInTests.cs or QuantificationTests.cs
- **Add new strategy classes** alongside existing ones in `Quantification/Strategies/`
- **New strategies must implement the existing interfaces**: `INormalizationStrategy`, `IRollUpStrategy`, `ICollapseStrategy`
- Follow the patterns established by `NoNormalization`, `SumRollUp`, and `SumCollapse`
- All new strategies go in the `Quantification.Strategies` namespace

## Key File Locations

```
Quantification/Strategies/Normalization/       ← Tasks 2, 3 add files here
Quantification/Strategies/RollUp/              ← Task 4 adds file here
Quantification/Strategies/Collapse/            ← Task 5 adds file here
Quantification/QuantificationEngine.cs         ← DO NOT MODIFY
Quantification/QuantificationParameters.cs     ← DO NOT MODIFY
Quantification/Interfaces/                     ← DO NOT MODIFY
Test/Quantification/TmtSpikeInTests.cs         ← Tasks 1, 6 add tests here
Test/Quantification/TestHelpers/QuantificationEvaluator.cs  ← Task 1 may extend
Test/Quantification/TestHelpers/PsmTsvQuantAdapter.cs       ← Existing adapter
Test/Quantification/TestHelpers/SpikeInExperimentalDesign.cs ← Existing design
TMT_Spike-In_Info/TMT_SpikeIn_Reference.md     ← Full experimental design reference
TMT_Spike-In_Info/UPS_TMT3_Search/Task1-SearchTask/AllPSMs.psmtsv  ← UPS spike-in data
TMT_Spike-In_Info/TMT3_HeLa_Search/Task1-SearchTask/AllPSMs.psmtsv ← HeLa background data
```

## Build and Test Commands

```bash
# Build
dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"

# Run quantification tests only
dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"

# Run TMT spike-in tests only
dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~TmtSpikeIn"
```

## Data Context

### Two AllPSMs files (searched separately, must be concatenated for full picture)

1. **UPS search** (`UPS_TMT3_Search/Task1-SearchTask/AllPSMs.psmtsv`): ~401k lines, 47 protein groups
   - UPS1 proteins spiked at known varying concentrations — these are the **true positives**
   - Expected fold changes: 1 vs 0.125 = 8.0x, 1 vs 0.5 = 2.0x, 1 vs 0.667 = 1.5x

2. **HeLa search** (`TMT3_HeLa_Search/Task1-SearchTask/AllPSMs.psmtsv`): ~1.35M lines, 3472 protein groups
   - SILAC-labeled HeLa background at constant concentration — these are the **negative controls**
   - Expected fold change between any conditions = 1.0 (no differential abundance)

### TMT10 channels (columns 58-67 in psmtsv)
126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N

### Experimental design
- 14 files x 10 channels = 140 samples
- Channels 126 and 131N are always reference channels
- Conditions: "Reference", "0.125", "0.5", "0.667", "1"
- See `SpikeInExperimentalDesign.cs` for the full mapping

### Existing pipeline (QuantificationEngine.RunTmt)
1. PivotByFile -> per-file PSM matrices
2. Per-file PSM normalization (currently NoNormalization)
3. Per-file PSM->peptide roll-up (currently SumRollUp)
4. CombinePeptideMatrices -> single peptide matrix (140 columns)
5. Peptide normalization (currently NoNormalization)
6. Collapse samples (currently NoCollapse)
7. Peptide->protein roll-up (currently SumRollUp)
8. Protein normalization (currently NoNormalization)

### Where each new strategy applies in the pipeline

| Strategy | Pipeline step(s) | Parameter slot(s) |
|----------|-----------------|-------------------|
| GlobalMedianNormalization | 2, 5, or 8 | SpectralMatchNormalizationStrategy, PeptideNormalizationStrategy, or ProteinNormalizationStrategy |
| ReferenceChannelNormalization | 2, 5, or 8 | Same as above |
| MedianRollUp | 3 or 7 | SpectralMatchToPeptideRollUpStrategy or PeptideToProteinRollUpStrategy |
| MeanCollapse | 6 | CollapseStrategy |
