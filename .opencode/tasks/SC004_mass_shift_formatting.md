# Task: Add Mass Shift Sequence Formatting

## Objective
Add methods to format sequences with mass shifts instead of modification names. This is useful for display, comparison, and tools that expect mass-based notation.

## Background

The existing `IBioPolymerWithSetMods` interface has:
- `FullSequence` - with full mod names: `PEPT[Common Fixed:Carbamidomethyl on C]IDE`
- `SequenceWithChemicalFormulas` - with formulas: `PEPT[C2H3NO]IDE`

We need to add:
- `FullSequenceWithMassShifts` - with masses: `PEPT[+57.0215]IDE`

This already exists as an extension method that needs to be verified/improved.

## Implementation Steps

- [ ] Review existing `FullSequenceWithMassShift()` extension in `IBioPolymerWithSetMods`
- [ ] Add `GetMassShiftNotation(Modification mod)` helper
- [ ] Add precision parameter for mass formatting (default 4 decimal places)
- [ ] Add option for signed notation (+57.0215 vs 57.0215)
- [ ] Add round-trip capability (parse mass shift back to find closest mod)
- [ ] Ensure consistent behavior for both peptides and oligos
- [ ] Add to `SequenceConverter` class

## File Locations
- `Omics/IBioPolymerWithSetMods.cs` - Extension methods
- `Omics/Modifications/SequenceConverter.cs` - Add methods

## Key Code

```csharp
// Extension method on IBioPolymerWithSetMods
public static string FullSequenceWithMassShift(
    this IBioPolymerWithSetMods bioPolymer,
    int decimalPlaces = 4,
    bool signed = true);

// In SequenceConverter
public static string ToMassShiftNotation(
    string fullSequence,
    Dictionary<string, Modification> modDictionary,
    int decimalPlaces = 4);

public static Dictionary<int, Modification> FromMassShiftNotation(
    string massShiftSequence,
    string baseSequence,
    char targetResidue,
    IList<Modification> candidateMods,
    double tolerance = 0.001);
```

## Mass Shift Format

```
Standard:    [+57.0215] or [-18.0106]
Unsigned:    [57.0215] or [18.0106]
Scientific:  [+5.7E+01] (not implemented, future enhancement)
```

## Acceptance Criteria
- [ ] Mass shift formatting works for all modification types
- [ ] Precision parameter controls decimal places
- [ ] Signed notation is correct for losses (negative mass)
- [ ] Works for terminal modifications
- [ ] Round-trip parsing finds correct mods within tolerance
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~SequenceConverter"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~MassShift"
```

## Test Cases

1. Single modification: `PEPTCIDE` + Carbamidomethyl on C -> `PEPT[+57.0215]IDE`
2. Multiple modifications: Phosphorylation + Oxidation
3. Negative mass modification (water loss, etc.)
4. N-terminal modification
5. C-terminal modification
6. Very small mass changes (label swaps)
7. Round-trip: Parse mass shift -> find closest mod -> verify same mod
