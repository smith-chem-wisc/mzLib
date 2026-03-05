# Task: Integrate SequenceConverter with ProteinDbWriter

## Objective
Integrate the `SequenceConverter` with `ProteinDbWriter` to ensure all modifications are converted to UniProt-compatible format when writing XML databases for use with other search engines.

## Background

`ProteinDbWriter.WriteXmlDatabase()` writes protein databases in mzLibProteinDb XML format. When proteins have modifications from MetaMorpheus searches (GPTMD results), these modifications need to be converted to UniProt format so other search engines can understand them.

Current behavior: Writes modifications as-is with their original names
Desired behavior: Option to convert all mods to UniProt naming convention

## Implementation Steps

- [ ] Review `UsefulProteomicsDatabases/ProteinDbWriter.cs`
- [ ] Add optional parameter for target `ModificationNamingConvention`
- [ ] Integrate `SequenceConverter` for modification conversion
- [ ] Update `GetModsForThisBioPolymer()` to convert mods
- [ ] Add fallback handling for unconvertible modifications
- [ ] Ensure modification cross-references are preserved
- [ ] Update both `WriteXmlDatabase` overloads (Protein and RNA)
- [ ] Add tests for converted output

## File Location
`UsefulProteomicsDatabases/ProteinDbWriter.cs`

## Key Changes

```csharp
public static class ProteinDbWriter
{
    /// <summary>
    /// Writes an XML database with optional modification conversion
    /// </summary>
    public static Dictionary<string, int> WriteXmlDatabase(
        Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
        List<Protein> proteinList,
        string outputFileName,
        bool updateTimeStamp = false,
        ModificationNamingConvention? targetConvention = null,  // NEW
        ConversionFailureMode failureMode = ConversionFailureMode.KeepOriginal) // NEW
    {
        // If targetConvention specified, convert modifications before writing
        var modsToWrite = targetConvention.HasValue
            ? ConvertModifications(additionalModsToAddToProteins, targetConvention.Value, failureMode)
            : additionalModsToAddToProteins;
            
        // ... existing write logic
    }
    
    private static Dictionary<string, HashSet<Tuple<int, Modification>>> ConvertModifications(
        Dictionary<string, HashSet<Tuple<int, Modification>>> mods,
        ModificationNamingConvention targetConvention,
        ConversionFailureMode failureMode)
    {
        var converter = SequenceConverter.Default;
        // Convert each modification...
    }
}
```

## Use Cases

1. **Export for MaxQuant/Proteome Discoverer**: Convert MetaMorpheus results to UniProt mods
2. **Database standardization**: Ensure all mods in database use consistent naming
3. **Cross-engine validation**: Create databases that work with multiple search engines

## Modification Format in XML

Current format:
```xml
<modification>
ID   Carbamidomethyl on C
MT   Common Fixed
...
</modification>
<feature type="modified residue" description="Carbamidomethyl on C">
```

UniProt format:
```xml
<modification>
ID   Carboxymethyl cysteine on C
MT   UniProt
AC   PTM-0100
...
</modification>
<feature type="modified residue" description="Carboxymethyl cysteine on C">
```

## Acceptance Criteria
- [ ] `WriteXmlDatabase` accepts optional `targetConvention` parameter
- [ ] Modifications are correctly converted to target convention
- [ ] Unconvertible mods handled per `failureMode`
- [ ] Database references (DR lines) preserved
- [ ] Output XML is valid and parseable
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~ProteinDbWriter"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~Database"
```

## Test Cases

1. Write database with MetaMorpheus mods -> read back and verify
2. Convert to UniProt convention -> verify mod names changed
3. Handle mod with no UniProt equivalent -> verify fallback behavior
4. Round-trip: Write -> Read -> Compare mods
