# Task 7: Write Synthetic Data Tests for TMT Pipeline

## Objective
Create unit tests using small synthetic data that validate the TMT pipeline mechanics: PivotByFile, per-file roll-up, CombinePeptideMatrices, and the full RunTmt pipeline.

## Prerequisites
- Tasks 4, 5, 6 must be complete (PivotByFile, CombinePeptideMatrices, RunTmt all exist)

## Background
The existing test file `Test/Quantification/QuantificationTests.cs` has helper methods and patterns for creating synthetic test data. Study this file to understand:
- How `TestExperimentalDesign` implements `IExperimentalDesign`
- How `IsobaricQuantSampleInfo` objects are created for TMT channels
- How `BaseSpectralMatch` objects are created with `QuantValues`
- How `Protein` objects are created and digested to get peptides
- How `BioPolymerGroup` objects wrap proteins

## File to Create

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Quantification/TmtSpikeInTests.cs`

## Test Class Structure

```csharp
using NUnit.Framework;
using Quantification;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using Omics;
using Omics.BioPolymerGroup;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.Quantification;

[TestFixture]
public class TmtSpikeInTests
{
    // Test helper methods and test cases below
}
```

## Tests to Implement

### Test 1: PivotByFile_CreatesCorrectPerFileMatrices

**Setup**: 2 files, 3 TMT channels each, 2-3 PSMs per file
**Assertions**:
- Result dictionary has 2 entries (one per file)
- Each matrix has correct number of rows (PSMs in that file)
- Each matrix has 3 columns (channels)
- Values match the synthetic QuantValues exactly
- No zeros in any cell (fully dense)

**Synthetic data pattern**:
```
File1:
  PSM1: QuantValues = [100, 200, 300]
  PSM2: QuantValues = [400, 500, 600]

File2:
  PSM3: QuantValues = [700, 800, 900]
  PSM4: QuantValues = [1000, 1100, 1200]
```

### Test 2: CombinePeptideMatrices_MergesCorrectly

**Setup**: Create 2 per-file peptide matrices manually (not through the full pipeline)
- File1 has Peptide A (100, 200, 300) and Peptide B (400, 500, 600)
- File2 has Peptide A (700, 800, 900) and Peptide C (1000, 1100, 1200)

**Assertions**:
- Combined matrix has 3 rows (Peptide A, B, C)
- Combined matrix has 6 columns (3 from file1 + 3 from file2)
- Peptide A: [100, 200, 300, 700, 800, 900]
- Peptide B: [400, 500, 600, 0, 0, 0]
- Peptide C: [0, 0, 0, 1000, 1100, 1200]

### Test 3: RunTmt_FullPipeline_ProducesCorrectProteinMatrix

**Setup**: 2 files, 3 channels each, 2 proteins with 2 unique peptides each, 1-2 PSMs per peptide
- Use NoNormalization, SumRollUp, NoCollapse
- Assign known QuantValues so expected protein-level sums are calculable

**Assertions**:
- Protein matrix has 2 rows (proteins)
- Protein matrix has 6 columns (3 channels × 2 files)
- Each protein's values = sum of its unique peptide values
- Each peptide's values = sum of its PSM values
- Verify a few specific cells against manually calculated expected values

### Test 4: RunTmt_WithSumCollapse_CombinesTechnicalReplicates

**Setup**: 3 files from the same mixture (same condition assignments, different technical replicates)
- Use SumCollapse to combine technical replicates
- 3 channels: Reference, ConditionA, ConditionB

**Assertions**:
- After collapse, the number of columns is reduced (from 9 to 3)
- Collapsed values = sum across technical replicates for each condition

## Helper Methods to Create

### CreateTmtExperimentalDesign
```csharp
/// <summary>
/// Creates a simple TMT experimental design for testing.
/// </summary>
/// <param name="fileNames">File names (without extension)</param>
/// <param name="channelLabels">TMT channel labels per file</param>
/// <param name="conditions">Condition string per channel</param>
/// <param name="bioReps">Biological replicate per channel</param>
/// <param name="techReps">Technical replicate per file</param>
private static IExperimentalDesign CreateTmtExperimentalDesign(
    string[] fileNames,
    string[] channelLabels,
    string[][] conditions,  // [fileIndex][channelIndex]
    int[][] bioReps,
    int[] techReps)
```

### CreateSyntheticTmtPsms
```csharp
/// <summary>
/// Creates synthetic PSMs with known QuantValues for TMT testing.
/// </summary>
/// <param name="filePath">File path for these PSMs</param>
/// <param name="numChannels">Number of TMT channels</param>
/// <param name="peptideSequences">Sequences for each PSM</param>
/// <param name="quantValues">Intensity values per PSM</param>
private static List<ISpectralMatch> CreateSyntheticTmtPsms(
    string filePath,
    string[] peptideSequences,
    double[][] quantValues,
    List<IBioPolymerWithSetMods> peptideObjects)
```

## Creating Test Objects

### Proteins and Peptides
Follow the pattern from existing `QuantificationTests.cs`:
```csharp
// Create proteins
var protein1 = new Protein("PEPTIDEK", "P001");
var protein2 = new Protein("ANOTHERPEPTIDEK", "P002");

// Digest to get peptides
var digestionParams = new DigestionParams(minPeptideLength: 5);
var peptides1 = protein1.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
var peptides2 = protein2.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

// Create protein groups
var group1 = new BioPolymerGroup(protein1, peptides1, new List<IBioPolymerWithSetMods>());
var group2 = new BioPolymerGroup(protein2, peptides2, new List<IBioPolymerWithSetMods>());
```

### PSMs with QuantValues
```csharp
var psm = new BaseSpectralMatch(
    fullFilePath: "file1",
    oneBasedScanNumber: 1,
    score: 100.0,
    fullSequence: "[TMT]PEPTIDEK",
    baseSequence: "PEPTIDEK",
    identifiedBioPolymers: new[] { peptides1[0] });
psm.QuantValues = new double[] { 100.0, 200.0, 300.0 };
```

### ExperimentalDesign
```csharp
var design = new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]>
{
    ["file1"] = new ISampleInfo[]
    {
        new IsobaricQuantSampleInfo("file1", "Reference", 1, 1, 0, 1, "126", 126.12776, true),
        new IsobaricQuantSampleInfo("file1", "CondA", 1, 1, 0, 1, "127N", 127.12476, false),
        new IsobaricQuantSampleInfo("file1", "CondB", 1, 1, 0, 1, "127C", 127.13108, false),
    },
    // ... more files
});
```

Note: The `TestExperimentalDesign` class is defined in `QuantificationTests.cs`. If it's internal/private, you may need to either make it accessible or create your own test implementation.

## Important: Using QuantificationParameters

```csharp
var parameters = QuantificationParameters.GetSimpleParameters();
// This gives: NoNormalization everywhere, SumRollUp, NoCollapse
```

For the collapse test, create custom parameters:
```csharp
var parameters = new QuantificationParameters
{
    SpectralMatchNormalizationStrategy = new NoNormalization(),
    SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
    PeptideNormalizationStrategy = new NoNormalization(),
    CollapseStrategy = new SumCollapse(),  // <-- Use SumCollapse here
    PeptideToProteinRollUpStrategy = new SumRollUp(),
    ProteinNormalizationStrategy = new NoNormalization(),
};
```

## Verification
1. All tests pass: `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~TmtSpikeIn"`
2. Existing tests still pass: `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"`
