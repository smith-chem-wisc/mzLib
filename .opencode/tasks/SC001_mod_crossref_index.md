# Task: Build Modification Cross-Reference Index

## Objective
Create a `ModificationCrossRefIndex` class that builds lookup tables allowing modifications to be found by their database cross-reference IDs (RESID, PSI-MOD, Unimod accession numbers). This is the foundation for accurate modification conversion between naming conventions.

## Background

The `Modification.DatabaseReference` property contains external database IDs:
```csharp
public Dictionary<string, IList<string>> DatabaseReference { get; protected set; }
```

Example entries from ptmlist.txt:
```
DR   RESID; AA0441.
DR   PSI-MOD; MOD:01624.
```

Example from Unimod (added by ModificationLoader):
```csharp
dblinks = new Dictionary<string, IList<string>>
{
    { "Unimod",  new List<string>{ac.ToString() } },
};
```

These cross-references allow us to definitively match "Phosphorylation on S" (MetaMorpheus) with "Phosphoserine" (UniProt) because they share the same RESID or PSI-MOD ID.

## Implementation Steps

- [ ] Create `Omics/Modifications/ModificationCrossRefIndex.cs`
- [ ] Add dictionary properties for each database type:
  - `ByResidId` - maps "AA0441" -> List<Modification>
  - `ByPsiModId` - maps "MOD:01624" -> List<Modification>
  - `ByUnimodId` - maps "35" -> List<Modification>
  - `ByAccession` - maps "PTM-0450" -> Modification (for UniProt accessions)
- [ ] Add constructor that builds indices from a list of modifications
- [ ] Add factory method `Build(IEnumerable<Modification> mods)` 
- [ ] Add static instance that indexes all known mods from `Mods.AllKnownMods`
- [ ] Add `TryGetEquivalent(Modification source, ModificationNamingConvention targetConvention, out Modification target)` method
- [ ] Handle edge cases: mods with no cross-refs, multiple mods sharing same ID

## File Location
`Omics/Modifications/Conversion/ModificationCrossRefIndex.cs`

## Key Code Structure

```csharp
namespace Omics.Modifications;

public class ModificationCrossRefIndex
{
    public Dictionary<string, List<Modification>> ByResidId { get; }
    public Dictionary<string, List<Modification>> ByPsiModId { get; }
    public Dictionary<string, List<Modification>> ByUnimodId { get; }
    public Dictionary<string, Modification> ByAccession { get; }
    
    // Lazy-initialized global instance
    public static ModificationCrossRefIndex Global { get; }
    
    public ModificationCrossRefIndex(IEnumerable<Modification> modifications);
    
    public bool TryGetEquivalent(
        Modification source, 
        ModificationNamingConvention targetConvention, 
        out Modification? target);
        
    public IEnumerable<Modification> GetEquivalentMods(Modification source);
}
```

## Acceptance Criteria
- [ ] Code compiles without errors
- [ ] Index correctly maps RESID, PSI-MOD, and Unimod IDs to modifications
- [ ] `TryGetEquivalent` successfully finds equivalent mods across conventions
- [ ] Edge cases handled (null refs, empty refs, multiple matches)
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test modifications
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~ModificationTest"
```

## Notes
- The ptmlist.txt format uses "RESID; AA0441" format in DR lines
- PSI-MOD IDs have format "MOD:XXXXX"
- Unimod IDs are just the record_id number as string
- Some mods may not have any cross-references - these can only be matched by mass/name heuristics
- Multiple UniProt mods may map to the same RESID (e.g., phosphorylation on different residues)
