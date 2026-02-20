# Task 1: Parse TMT Reporter Ion Columns in SpectrumMatchFromTsv

## Objective
Extend the `SpectrumMatchFromTsv` class to detect and parse TMT reporter ion intensity columns from .psmtsv files. Store the parsed values in a new `QuantValues` property.

## Background
MetaMorpheus writes .psmtsv files with optional TMT reporter ion columns at the end. The column names for TMT10-plex are: `126`, `127N`, `127C`, `128N`, `128C`, `129N`, `129C`, `130N`, `130C`, `131`. These columns contain `double` intensity values (or 0 if not detected). Currently, these columns are ignored by the reader.

## Files to Modify

### 1. `C:/Users/Alex/Source/Repos/mzLib/mzLib/Readers/InternalResults/IndividualResultRecords/SpectrumMatchFromTsvHeader.cs`

Add TMT channel name constants. Add them in a new section at the bottom, before the closing brace:

```csharp
// TMT/Isobaric reporter ion channels
public static readonly string[] TmtChannelNames = new[]
{
    "126", "127N", "127C", "128N", "128C",
    "129N", "129C", "130N", "130C", "131",
    "131C", "132N", "132C", "133N", "133C", "134N" // TMT16+ channels
};
```

### 2. `C:/Users/Alex/Source/Repos/mzLib/mzLib/Readers/InternalResults/IndividualResultRecords/SpectrumMatchFromTsv.cs`

Add a `QuantValues` property and parse TMT columns in the constructor.

**Add property** (near the other properties, around line 63):
```csharp
public double[]? QuantValues { get; protected set; }
```

**Add parsing logic** at the end of the constructor (the one that takes `string line, char[] split, Dictionary<string, int> parsedHeader`, around line 183, after `SpectralAngle` parsing):

```csharp
// Parse TMT/isobaric reporter ion columns if present
QuantValues = ParseReporterIonColumns(spl, parsedHeader);
```

**Add the helper method** as a private static method:

```csharp
/// <summary>
/// Detects TMT reporter ion columns in the header and parses their values.
/// Returns null if no TMT columns are found.
/// </summary>
private static double[]? ParseReporterIonColumns(string[] spl, Dictionary<string, int> parsedHeader)
{
    // Find which TMT channels are present in the header, preserving order
    var presentChannels = new List<(string name, int index)>();
    foreach (var channelName in SpectrumMatchFromTsvHeader.TmtChannelNames)
    {
        if (parsedHeader.TryGetValue(channelName, out int colIndex) && colIndex >= 0 && colIndex < spl.Length)
        {
            presentChannels.Add((channelName, colIndex));
        }
    }

    if (presentChannels.Count == 0)
        return null;

    double[] values = new double[presentChannels.Count];
    for (int i = 0; i < presentChannels.Count; i++)
    {
        if (double.TryParse(spl[presentChannels[i].index].Trim(),
            System.Globalization.NumberStyles.Any,
            System.Globalization.CultureInfo.InvariantCulture,
            out double val))
        {
            values[i] = val;
        }
        // else stays 0.0
    }
    return values;
}
```

## Important Notes

- The `parsedHeader` dictionary is populated by the file reader from the header row. It maps column name → column index. If a column is not present, it maps to -1.
- The existing pattern for optional columns uses `(parsedHeader[key] < 0) ? null : ...`. However, TMT columns are not in the header dictionary by default - they would only be present if the file has them. Use `TryGetValue` to check.
- The `parsedHeader` dictionary is built externally. You need to verify that the dictionary builder includes ALL columns from the header, not just the known ones. Check how `parsedHeader` is constructed - if it only includes predefined column names, the TMT columns won't appear. In that case, the dictionary builder also needs to be updated to include any unrecognized columns.

## How parsedHeader is Built
Look at how PsmFromTsv files are read. The header parsing typically happens in a method that reads the first line and creates the dictionary. Search for where `parsedHeader` is created (likely in a reader class or a static method). The TMT column names ("126", "127N", etc.) need to be present as keys in this dictionary for the parsing to work.

## Verification
1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` - must compile
2. `dotnet test "C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Test.csproj" --filter "FullyQualifiedName~Quantification"` - existing tests must pass
3. Verify manually: Read a .psmtsv file from `TMT_Spike-In_Info/UPS_Search/Task1-SearchTask/AllPSMs.psmtsv`, check that the header contains columns "126" through "131", and confirm the parser detects them
