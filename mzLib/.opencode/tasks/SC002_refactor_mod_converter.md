# Task: Refactor ModificationConverter with Database IDs

## Objective
Refactor `ModificationConverter.cs` to use the `ModificationCrossRefIndex` for accurate modification matching before falling back to name/mass heuristics. This will dramatically improve conversion accuracy.

## Background

The current `ModificationConverter.GetClosestMod()` method uses a heuristic approach:
1. Find mods where name contains the search term
2. Find mods where target residue matches
3. Intersect and score by "overlap" (longest common substring)

This is error-prone because:
- "Oxidation" could match "Dioxidation" 
- Similar names don't guarantee same modification
- Mass tolerance matching is imprecise

Better approach: Use database cross-references first, then fall back to heuristics.

## Implementation Steps

- [ ] Add dependency on `ModificationCrossRefIndex` in `ModificationConverter`
- [ ] Create new method `TryGetExactMatch()` that uses cross-reference IDs
- [ ] Modify `GetClosestMod()` to try exact match first
- [ ] Add method to extract database IDs from modification names/strings
- [ ] Update `GetModificationDictionaryFromFullSequence()` to use improved matching
- [ ] Add unit tests for the improved matching logic
- [ ] Ensure backward compatibility - existing behavior should still work

## File Location
`Readers/ExternalResults/SupportClasses/ModificationConverter.cs`

## Key Changes

```csharp
public static class ModificationConverter
{
    // Add reference to cross-ref index
    private static readonly Lazy<ModificationCrossRefIndex> _crossRefIndex = 
        new(() => ModificationCrossRefIndex.Global);
    
    /// <summary>
    /// Try to find an exact modification match using database cross-references
    /// </summary>
    public static bool TryGetExactMatch(
        Modification source, 
        ModificationNamingConvention targetConvention,
        out Modification? target)
    {
        return _crossRefIndex.Value.TryGetEquivalent(source, targetConvention, out target);
    }
    
    /// <summary>
    /// Enhanced GetClosestMod that tries exact match first
    /// </summary>
    public static Modification GetClosestMod(string name, char modifiedResidue, 
        IList<Modification>? allKnownMods = null,
        ModificationNamingConvention? preferredConvention = null)
    {
        // 1. Try to find exact match via cross-references
        // 2. Fall back to existing heuristic matching
    }
}
```

## Matching Priority

1. **By Database ID** - If the mod name contains a known database ID format:
   - "UNIMOD:35" -> Look up in Unimod index
   - "MOD:00046" -> Look up in PSI-MOD index
   - "PTM-0001" -> Look up in UniProt accession index
   
2. **By Cross-Reference** - If source mod has DatabaseReference, find target with same refs

3. **By Mass + Residue** - Match by monoisotopic mass (within tolerance) and target residue

4. **By Name Heuristic** - Existing overlap scoring as last resort

## Acceptance Criteria
- [ ] `TryGetExactMatch()` correctly finds equivalent mods via database IDs
- [ ] `GetClosestMod()` tries exact match before heuristics
- [ ] Existing functionality preserved (backward compatible)
- [ ] Unit tests pass
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~ModificationConverter"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~ModificationTest"
```

## Test Cases to Add

1. Convert "Phosphorylation on S" (MetaMorpheus) -> "Phosphoserine" (UniProt)
2. Convert "Oxidation on M" (MetaMorpheus) -> "Oxidation" (Unimod)
3. Convert UniProt mod back to MetaMorpheus format
4. Handle mod with no cross-references (should fall back to heuristics)
5. Handle ambiguous matches (multiple mods with same cross-ref)
