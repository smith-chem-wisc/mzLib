# Task 8: Write Tests Using Real Spike-In Data

## Objective
Create integration tests that load the real SpikeIn-5mix-MS3 .psmtsv files, run the TMT quantification pipeline, and evaluate the results against known ground truth fold changes.

## Prerequisites
- Tasks 1, 2, 3, 6 must all be complete
- Task 7 should pass (synthetic tests validate pipeline mechanics)

## Background
The SpikeIn-5mix-MS3 dataset has known UPS1 protein fold changes and constant HeLa background. By running the quantification pipeline on this real data, we can measure how well different strategies recover the expected ratios.

The current .psmttsv files have non-zero TMT reporter ion intensities. Previously the data had all zero intensities, but this issue has been resolved. Data can be found in the TMT_Spike-In_Info folder

## File to Modify

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Quantification/TmtSpikeInTests.cs`
(Created in Task 7 - add new tests to this file)

## Data Locations
- UPS search results: `TMT_Spike-In_Info/UPS_TMT3_Search/Task1-SearchTask/AllPSMs.psmtsv`
- HeLa search results: `TMT_Spike-In_Info/TMT3_HELA_Search/Task1-SearchTask/AllPSMs.psmtsv`
- Reference doc: `TMT_Spike-In_Info/TMT_SpikeIn_Reference.md`

## Expected Fold Changes (UPS1 proteins)

| Comparison | True Fold Change |
|------------|-----------------|
| 1.0 vs 0.667 | 1.5 |
| 1.0 vs 0.5 | 2.0 |
| 1.0 vs 0.125 | 8.0 |
| 0.667 vs 0.5 | 1.33 |
| 0.667 vs 0.125 | 5.328 |
| 0.5 vs 0.125 | 4.0 |

For HeLa proteins: all fold changes should be ~1.0.

## Tests to Implement

### Test 1: LoadAndRunSpikeInData_BasicPipeline
```csharp
[Test]
// [Ignore("TMT intensities are zero")] // Uncomment if data has zero intensities
public void LoadAndRunSpikeInData_BasicPipeline()
{
    // 1. Load PSMs from UPS search psmtsv
    //    Use PsmTsvQuantAdapter.LoadSpectralMatches(path)

    // 2. Create experimental design
    //    Use SpikeInExperimentalDesign from Task 3

    // 3. Build peptide/protein lists from the loaded PSMs
    //    - Extract unique FullSequence → peptide objects
    //    - Extract unique Accession → protein groups
    //    NOTE: Since we're reading from TSV, we don't have actual Protein objects.
    //    We need to create lightweight protein/peptide representations.
    //    See "Creating Protein/Peptide Objects from PSM Data" section below.

    // 4. Create QuantificationEngine with simple parameters
    //    NoNormalization + SumRollUp + NoCollapse

    // 5. Run TMT pipeline
    //    var proteinMatrix = engine.RunTmtAndReturnProteinMatrix();

    // 6. Basic sanity checks:
    //    - proteinMatrix is not null
    //    - Number of protein rows > 0
    //    - Number of columns = 14 files × 10 channels = 140
    //    - At least some values are non-zero

    Assert.IsNotNull(proteinMatrix);
    Assert.That(proteinMatrix.RowKeys.Count, Is.GreaterThan(0));
}
```

### Test 2: EvaluateFoldChanges_UPSProteins
```csharp
[Test]
// [Ignore("TMT intensities are zero")]
public void EvaluateFoldChanges_UPSProteins()
{
    // ... (load data, run pipeline as above) ...

    // For each protein in the result:
    // 1. Group columns by condition using the experimental design
    // 2. Calculate mean intensity per condition (across biological replicates and technical replicates)
    // 3. Calculate fold change = mean(condition1) / mean(condition2)
    // 4. Compare to expected fold change

    // Expected: UPS proteins show fold changes close to expected values
    // Expected: direction of fold changes is correct (higher concentration = higher intensity)
}
```

### Test 3: EvaluateFoldChanges_HelaBackground
```csharp
[Test]
// [Ignore("TMT intensities are zero")]
public void EvaluateFoldChanges_HelaBackground()
{
    // ... (load HeLa data, run pipeline) ...

    // For HeLa proteins, fold changes between conditions should be ~1.0
    // Calculate CV (coefficient of variation) across conditions
    // Assert that median fold change is close to 1.0
}
```

## Creating Protein/Peptide Objects from PSM Data

Since we're reading from .psmtsv files (not re-running a search), we don't have actual `Protein` objects. We need to create them from the PSM metadata.

**Approach**:
1. Parse PSMs from TSV
2. Group PSMs by Accession to identify protein groups
3. Group PSMs by FullSequence to identify unique peptides
4. Create `Protein` objects with dummy sequences (or use BaseSequence)
5. Create `PeptideWithSetModifications` by digesting the Protein
6. Create `BioPolymerGroup` objects

**Alternatively** (simpler for testing):
- Create lightweight mock objects that implement the required interfaces
- The key requirement is that `ISpectralMatch.GetIdentifiedBioPolymersWithSetMods()` returns a peptide that exists in the `modifiedBioPolymers` list, and that peptide maps to a protein in the `bioPolymerGroups` list

**Practical implementation**:
```csharp
// Group PSMs by base sequence → these become our "peptides"
// Group peptides by accession → these become our "protein groups"
// For each accession:
//   Create a Protein with the accession
//   For each unique base sequence mapping to that accession:
//     Create a PeptideWithSetModifications
//   Create a BioPolymerGroup with the protein and its peptides
// Then go back and set each BaseSpectralMatch's identified biopolymer
```

This is the most complex part of this task. Study how `CreateTestProteins()` and `CreateTestSpectralMatches()` work in `QuantificationTests.cs` for patterns.

## Helper: QuantificationEvaluator

Create a static helper class for computing evaluation metrics:

```csharp
public static class QuantificationEvaluator
{
    /// <summary>
    /// Calculates mean intensity per condition from a protein matrix.
    /// </summary>
    public static Dictionary<string, double> GetMeanIntensityByCondition(
        QuantMatrix<IBioPolymerGroup> proteinMatrix,
        IBioPolymerGroup protein)
    {
        // Get the row for this protein
        // Group columns by Condition
        // Calculate mean intensity per condition
    }

    /// <summary>
    /// Calculates fold change between two conditions for a protein.
    /// </summary>
    public static double CalculateFoldChange(
        Dictionary<string, double> meanIntensities,
        string numeratorCondition,
        string denominatorCondition)
    {
        return meanIntensities[numeratorCondition] / meanIntensities[denominatorCondition];
    }
}
```

## File Paths for Test Data

The test needs to find the .psmtsv files. Use a pattern like:
```csharp
// Navigate from test execution directory to the data
var testDir = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
var solutionDir = Path.GetFullPath(Path.Combine(testDir, "..", "..", "..", ".."));
var dataDir = Path.Combine(solutionDir, "TMT_Spike-In_Info");
var upsPsmPath = Path.Combine(dataDir, "UPS_Search", "Task1-SearchTask", "AllPSMs.psmtsv");
```

Or check how other tests in the project find test data files.

## Verification
1. All tests pass (or are correctly [Ignore]'d if data has zero intensities)
2. `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~TmtSpikeIn"`
3. No regression in existing tests
