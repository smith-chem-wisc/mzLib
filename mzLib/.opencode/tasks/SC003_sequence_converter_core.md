# Task: Create SequenceConverter Core Class

## Objective
Create a `SequenceConverter` class that converts `IBioPolymerWithSetMods` full sequences between different modification naming conventions. This is the main entry point for sequence conversion operations.

## Background

Different tools use different modification naming conventions:
- **MetaMorpheus**: `[Common Fixed:Carbamidomethyl on C]`
- **UniProt**: `[UniProt:Carboxymethyl cysteine on C]`
- **Unimod**: `[Unimod:Carbamidomethyl on C]`

Users need to:
1. Convert sequences for use with external tools (other search engines)
2. Standardize sequences from different sources
3. Display sequences with mass shifts instead of mod names

## Implementation Steps

- [ ] Create `Omics/Modifications/SequenceConverter.cs`
- [ ] Add `ConvertFullSequence(string fullSequence, ModificationNamingConvention source, ModificationNamingConvention target)` method
- [ ] Add `ConvertFullSequence(IBioPolymerWithSetMods bioPolymer, ModificationNamingConvention target)` method
- [ ] Add `ConvertModifications(Dictionary<int, Modification> mods, ModificationNamingConvention target)` method
- [ ] Handle unconvertible modifications (no equivalent in target convention)
- [ ] Add options for handling conversion failures (throw, skip, use mass shift)
- [ ] Add extension methods on `IBioPolymerWithSetMods` for convenient access

## File Location
`Omics/Modifications/SequenceConverter.cs`

## Key Code Structure

```csharp
namespace Omics.Modifications;

public enum ConversionFailureMode
{
    Throw,              // Throw exception if mod can't be converted
    Skip,               // Skip unconvertible mods (dangerous - changes mass)
    UseMassShift,       // Replace with mass shift notation [+123.456]
    KeepOriginal        // Keep original mod if no conversion found
}

public class SequenceConverter
{
    private readonly ModificationCrossRefIndex _crossRefIndex;
    
    public SequenceConverter(ModificationCrossRefIndex? crossRefIndex = null);
    
    public static SequenceConverter Default { get; }
    
    /// <summary>
    /// Convert a full sequence string to a different naming convention
    /// </summary>
    public string ConvertFullSequence(
        string fullSequence,
        ModificationNamingConvention sourceConvention,
        ModificationNamingConvention targetConvention,
        ConversionFailureMode failureMode = ConversionFailureMode.UseMassShift);
    
    /// <summary>
    /// Convert modifications on a BioPolymer to a different naming convention
    /// </summary>
    public Dictionary<int, Modification> ConvertModifications(
        Dictionary<int, Modification> modifications,
        ModificationNamingConvention targetConvention,
        ConversionFailureMode failureMode = ConversionFailureMode.UseMassShift);
    
    /// <summary>
    /// Get full sequence with mods from a specific convention
    /// </summary>
    public string GetFullSequence(
        IBioPolymerWithSetMods bioPolymer,
        ModificationNamingConvention targetConvention,
        ConversionFailureMode failureMode = ConversionFailureMode.UseMassShift);
}

// Extension methods for convenience
public static class SequenceConverterExtensions
{
    public static string ToUniProtSequence(this IBioPolymerWithSetMods bioPolymer);
    public static string ToUnimodSequence(this IBioPolymerWithSetMods bioPolymer);
    public static string ToMetaMorpheusSequence(this IBioPolymerWithSetMods bioPolymer);
}
```

## Sequence Format

Input/Output format follows `IBioPolymerWithSetMods.DetermineFullSequence()`:
```
PEPTIDE[ModType:ModName on X]SEQUENCE
```

Examples:
```
Input (MetaMorpheus):  PEPT[Common Fixed:Carbamidomethyl on C]IDE
Output (UniProt):      PEPT[UniProt:Carboxymethyl cysteine on C]IDE
Output (Mass Shift):   PEPT[+57.021464]IDE
```

## Acceptance Criteria
- [ ] `ConvertFullSequence()` correctly converts between conventions
- [ ] Handles N-terminal and C-terminal modifications
- [ ] `ConversionFailureMode` options work as expected
- [ ] Extension methods provide convenient access
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~SequenceConverter"
```

## Edge Cases to Handle

1. Modifications with no equivalent in target convention
2. Terminal modifications (N-term, C-term)
3. Multiple modifications on same position (shouldn't happen but handle gracefully)
4. Empty/null sequences
5. RNA vs Protein modifications (different conventions)
6. Modifications that exist in both conventions with same name
