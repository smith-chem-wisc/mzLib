# Task 2: Create PsmTsvQuantAdapter - Convert PsmFromTsv to ISpectralMatch

## Objective
Create an adapter class that reads a .psmtsv file and converts the records into `ISpectralMatch` objects (specifically `BaseSpectralMatch`) with `QuantValues` populated from TMT reporter ion columns.

## Prerequisites
- Task 1 must be complete (TMT columns parsed into `SpectrumMatchFromTsv.QuantValues`)

## Background

### BaseSpectralMatch (in `Omics/BaseSpectralMatch.cs`)
```csharp
public class BaseSpectralMatch : ISpectralMatch
{
    public BaseSpectralMatch(
        string fullFilePath,
        int oneBasedScanNumber,
        double score,
        string fullSequence,
        string baseSequence,
        IEnumerable<IBioPolymerWithSetMods>? identifiedBioPolymers = null)

    public double[]? QuantValues { get; set; }  // Can be set after construction
    public void AddIdentifiedBioPolymer(IBioPolymerWithSetMods bioPolymer)
}
```

### SpectrumMatchFromTsv properties we need
- `FileNameWithoutExtension` → maps to `FullFilePath`
- `Ms2ScanNumber` → maps to `OneBasedScanNumber`
- `Score` → maps to `Score`
- `FullSequence` → maps to `FullSequence`
- `BaseSeq` → maps to `BaseSequence`
- `Accession` → used to map to protein objects
- `DecoyContamTarget` → used to filter targets only
- `QValue` → used to filter by FDR
- `QuantValues` → copied to `BaseSpectralMatch.QuantValues` (from Task 1)

## File to Create

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Quantification/TestHelpers/PsmTsvQuantAdapter.cs`

```csharp
using Omics;
using Omics.BioPolymerGroup;
using Readers;

namespace Test.Quantification.TestHelpers;

/// <summary>
/// Reads a .psmtsv file and converts records to ISpectralMatch objects
/// suitable for the Quantification pipeline.
/// </summary>
public static class PsmTsvQuantAdapter
{
    /// <summary>
    /// Reads a psmtsv file and returns ISpectralMatch objects with QuantValues populated.
    /// </summary>
    /// <param name="psmtsvFilePath">Path to the AllPSMs.psmtsv file</param>
    /// <param name="qValueCutoff">Maximum q-value for filtering (default 0.01)</param>
    /// <param name="includeDecoys">Whether to include decoy matches (default false)</param>
    /// <returns>List of ISpectralMatch objects with QuantValues set</returns>
    public static List<ISpectralMatch> LoadSpectralMatches(
        string psmtsvFilePath,
        double qValueCutoff = 0.01,
        bool includeDecoys = false)
    {
        // 1. Read the psmtsv file using the existing PsmFromTsv infrastructure
        //    Look at how other tests or tools read psmtsv files in the codebase
        //    The pattern is typically:
        //    - Read header line, build parsedHeader dictionary
        //    - For each subsequent line, create PsmFromTsv(line, split, parsedHeader)

        // 2. Filter: targets only (DecoyContamTarget contains "T"), QValue <= cutoff

        // 3. For each PsmFromTsv, create BaseSpectralMatch:
        //    - fullFilePath = record.FileNameWithoutExtension (note: this needs to include
        //      enough of the path to match the ExperimentalDesign file name keys)
        //    - Set QuantValues from the parsed TMT columns

        // 4. Return the list

        throw new NotImplementedException("Implement this method");
    }

    /// <summary>
    /// Extracts unique protein accessions from PSM records and creates a mapping.
    /// </summary>
    public static HashSet<string> GetUniqueAccessions(string psmtsvFilePath)
    {
        // Read the file and collect all unique Accession values
        throw new NotImplementedException("Implement this method");
    }
}
```

## Implementation Details

### Reading PsmFromTsv files
Search the codebase for how psmtsv files are read. The typical pattern involves:
1. Reading the first line as header
2. Splitting by `\t` to get column names
3. Building a `Dictionary<string, int>` mapping column name → index
4. For each data line, creating a `PsmFromTsv(line, new[] { '\t' }, parsedHeader)`

Look in the `Readers` project for a `SpectrumMatchTsvReader` or similar class, or search for `PsmFromTsv` constructor calls in tests.

### File Name Handling
The `SpectrumMatchFromTsv.FileNameWithoutExtension` strips the extension. The `ISpectralMatch.FullFilePath` is used by the QuantificationEngine to match against `IExperimentalDesign.FileNameSampleInfoDictionary` keys. The design dictionary uses file names WITH extension. So when creating the BaseSpectralMatch, you need to decide on a convention.

**Recommended approach**: Use the `FileNameWithoutExtension` value as the `FullFilePath`. In the ExperimentalDesign (Task 3), use matching keys. The engine calls `Path.GetFileName(filePath)` to do the lookup.

### Protein Mapping
For the initial implementation, we don't need to map PSMs to actual Protein objects with digested peptides. The adapter just needs to create BaseSpectralMatch objects. The protein/peptide mapping will be handled separately in the test harness (Task 8). For now, leave `identifiedBioPolymers` as null in the BaseSpectralMatch constructor - the test will add them later.

## Verification
1. Build compiles: `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"`
2. Write a small unit test that reads the first 10 lines of `TMT_Spike-In_Info/UPS_Search/Task1-SearchTask/AllPSMs.psmtsv` and verifies that:
   - The correct number of PSMs are returned
   - Each has a non-null QuantValues array of length 10
   - FileNameWithoutExtension is populated correctly
