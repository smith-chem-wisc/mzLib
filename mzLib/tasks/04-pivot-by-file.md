# Task 4: Add PivotByFile Method to QuantificationEngine

## Objective
Add a new `PivotByFile()` static method to `QuantificationEngine` that creates one `SpectralMatchMatrix` per file. This is the first step in the TMT-specific pipeline that avoids the sparse matrix problem of the existing `Pivot()` method.

## Background
In TMT data, each PSM belongs to exactly one file and has intensities for the channels within that file only. The existing `Pivot()` creates a single matrix where columns span ALL channels across ALL files, making most cells zero. The TMT pipeline needs dense per-file matrices for within-file normalization.

## File to Modify

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Quantification/QuantificationEngine.cs`

**Do NOT modify** the existing `Pivot()` or `Run()` methods. Add the new method alongside them.

## Method to Add

```csharp
/// <summary>
/// Creates one SpectralMatchMatrix per file for TMT/isobaric data.
/// Each matrix has rows = PSMs from that file, columns = channels within that file.
/// This produces dense matrices (no sparse zeros) suitable for within-file normalization.
/// </summary>
/// <param name="spectralMatches">All spectral matches across all files</param>
/// <param name="experimentalDesign">Maps file names to channel ISampleInfo arrays</param>
/// <returns>Dictionary mapping file path to its SpectralMatchMatrix</returns>
public static Dictionary<string, SpectralMatchMatrix> PivotByFile(
    List<ISpectralMatch> spectralMatches,
    IExperimentalDesign experimentalDesign)
{
    // Implementation steps:
    //
    // 1. Filter to spectral matches with non-null QuantValues
    //
    // 2. Group by FullFilePath
    //
    // 3. For each file group:
    //    a. Get the ISampleInfo[] for this file from experimentalDesign
    //       Use Path.GetFileName(filePath) to look up in FileNameSampleInfoDictionary
    //    b. Get the list of PSMs for this file, ordered by FullSequence
    //    c. Create a SpectralMatchMatrix with:
    //       - rowKeys = the PSMs from this file
    //       - columnKeys = the ISampleInfo[] for this file (typically 10 for TMT10)
    //       - experimentalDesign = the shared design object
    //    d. For each PSM in this file:
    //       - Copy QuantValues directly into the matrix row
    //       - QuantValues[i] corresponds to columnKeys[i] (same positional mapping)
    //       - Use smMatrix.SetRow(spectralMatch, quantValues)
    //    e. Add to result dictionary: filePath → matrix
    //
    // 4. Return the dictionary
}
```

## Key Differences from Existing Pivot()

| Aspect | Existing Pivot() | New PivotByFile() |
|--------|-----------------|-------------------|
| Output | Single SpectralMatchMatrix | Dictionary of SpectralMatchMatrix per file |
| Columns | All channels from all files | Only channels from one file |
| Density | Sparse (most values 0) | Dense (all values meaningful) |
| QuantValues mapping | Needs filePathToArrayPositionDict | Direct positional copy |

## Implementation Notes

1. **Column keys per file**: For a TMT10 experiment, each file has exactly 10 ISampleInfo objects (one per channel). The `experimentalDesign.FileNameSampleInfoDictionary[fileName]` gives these in the correct order matching `QuantValues`.

2. **QuantValues direct mapping**: Since `QuantValues[0]` corresponds to `ISampleInfo[0]` for the file, and we're building a per-file matrix where columns ARE the per-file ISampleInfo array, the `QuantValues` can be used directly as the row values. No index remapping needed.

3. **File name matching**: The `SpectrumMatchFromTsv` stores `FileNameWithoutExtension`. The `IExperimentalDesign` may use file names with or without extensions. Use `Path.GetFileName()` and handle both cases. Look at how the existing `GetOrderedSampleInfos()` handles this.

4. **SpectralMatchMatrix constructor**: Located in `Quantification/QuantMatrix.cs`:
   ```csharp
   public SpectralMatchMatrix(
       ICollection<ISpectralMatch> rowKeys,
       ICollection<ISampleInfo> columnKeys,
       IExperimentalDesign experimentalDesign)
       : base(rowKeys, columnKeys, experimentalDesign) { }
   ```

5. **Existing reference code**: Study the existing `Pivot()` method at lines 213-249 and `GetOrderedSampleInfos()` at lines 258-288 for patterns on file name resolution.

## Verification
1. Build compiles: `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"`
2. Existing tests pass: `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"`
3. The method should be testable in Task 7 (synthetic TMT tests)
