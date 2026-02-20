# Task 5: Add CombinePeptideMatrices Method to QuantificationEngine

## Objective
Add a `CombinePeptideMatrices()` static method to `QuantificationEngine` that merges per-file peptide matrices into a single combined peptide matrix spanning all files and channels.

## Prerequisites
- Task 4 must be complete (PivotByFile exists)

## Background
After PivotByFile creates per-file PSM matrices and they are normalized and rolled up to peptides within each file, we need to combine them. The combined matrix has:
- **Rows** = union of all peptides across all files
- **Columns** = all channels from all files (e.g., 14 files × 10 channels = 140 columns)
- **Values** = intensity where the peptide was observed in that file/channel; 0.0 otherwise

## File to Modify

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Quantification/QuantificationEngine.cs`

## Method to Add

```csharp
/// <summary>
/// Combines per-file peptide matrices into a single matrix spanning all files.
/// Rows = union of all peptides across files.
/// Columns = all channels from all files, ordered by file path then channel.
/// Values = peptide intensity in that channel, or 0 if the peptide was not observed in that file.
/// </summary>
/// <param name="perFilePeptideMatrices">Dictionary of file path → peptide matrix for that file</param>
/// <param name="experimentalDesign">The experimental design</param>
/// <returns>A single PeptideMatrix covering all files and channels</returns>
public static PeptideMatrix CombinePeptideMatrices(
    Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>> perFilePeptideMatrices,
    IExperimentalDesign experimentalDesign)
{
    // Implementation steps:
    //
    // 1. Collect all unique peptides across all per-file matrices
    //    Use a HashSet<IBioPolymerWithSetMods> or similar
    //    Then convert to an ordered list (order by FullSequence or similar for determinism)
    //
    // 2. Collect all column keys (ISampleInfo) across all files
    //    Order by file path (alphabetically), then by channel within each file
    //    This matches the pattern used by GetOrderedSampleInfos()
    //
    // 3. Create the combined PeptideMatrix:
    //    new PeptideMatrix(allPeptides, allColumnKeys, experimentalDesign)
    //
    // 4. For each per-file matrix:
    //    a. Determine the column offset for this file in the combined matrix
    //    b. For each peptide in the per-file matrix:
    //       - Get the row values from the per-file matrix
    //       - Place them at the correct column positions in the combined matrix
    //       - Use matrix indexing or SetRow with a pre-built full-width array
    //
    // 5. Return the combined matrix
}
```

## Implementation Details

### Determining Column Offsets
The combined matrix columns are ordered: [file1_ch0, file1_ch1, ..., file1_ch9, file2_ch0, ...].

To place per-file values into the correct positions:
1. Sort file paths alphabetically (same order as `GetOrderedSampleInfos`)
2. Track cumulative column offset
3. For file at offset N with K channels, its channels map to combined columns [N, N+1, ..., N+K-1]

### Handling Peptides Not Seen in All Files
A peptide observed in file A but not file B will have:
- Actual intensities in columns belonging to file A
- Zeros in columns belonging to file B (the default for a new matrix)

Since `DenseMatrix` initializes to 0.0, you only need to set values where the peptide WAS observed.

### PeptideMatrix Type
Located in `Quantification/QuantMatrix.cs`:
```csharp
public class PeptideMatrix : QuantMatrix<IBioPolymerWithSetMods>
{
    public PeptideMatrix(
        ICollection<IBioPolymerWithSetMods> rowKeys,
        ICollection<ISampleInfo> columnKeys,
        IExperimentalDesign experimentalDesign)
        : base(rowKeys, columnKeys, experimentalDesign) { }
}
```

### Setting Values in the Combined Matrix
Two approaches:
1. **Row-at-a-time**: For each peptide, build a full-width double[] with 0s, fill in the per-file values at the correct offsets, then call `SetRow(peptide, fullArray)`.
2. **Direct matrix indexing**: Access `Matrix[rowIndex, colIndex]` directly on the DenseMatrix.

Approach 1 is simpler. Note that if a peptide appears in multiple files, you need to set values from each file. Since `SetRow` overwrites the entire row, you should build the combined row from ALL files before calling SetRow.

**Recommended approach**: Build a `double[]` per peptide, iterate over all files that contain it, fill in the values at the correct column offsets, then SetRow once.

### IBioPolymerWithSetMods Equality
The per-file matrices may contain the same peptide object (same reference) or equivalent peptide objects (same sequence/mods). The roll-up step maps PSMs to peptides using the `modifiedBioPolymers` list, so the same peptide object instance is used across files if the same peptide list is shared. Verify this in tests.

## Verification
1. Build compiles
2. Existing tests pass
3. Will be tested in Task 7 with synthetic data
