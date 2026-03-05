# Task: Add Unit Tests for SequenceConverter

## Objective
Create comprehensive unit tests for the `SequenceConverter` and `ModificationCrossRefIndex` classes to ensure correct behavior across all use cases.

## Implementation Steps

- [ ] Create `Test/Omics/SequenceConverterTests.cs`
- [ ] Add tests for `ModificationCrossRefIndex`:
  - Index building from modification lists
  - Lookup by RESID, PSI-MOD, Unimod IDs
  - `TryGetEquivalent()` functionality
  - Edge cases (null refs, no matches)
- [ ] Add tests for `SequenceConverter`:
  - Full sequence conversion between conventions
  - Mass shift notation formatting
  - Round-trip conversions
  - Failure mode handling
- [ ] Add integration tests with real modifications
- [ ] Add performance tests for large-scale conversions

## File Location
`Test/Omics/SequenceConverterTests.cs`

## Test Structure

```csharp
namespace Test.Omics;

[TestFixture]
public class ModificationCrossRefIndexTests
{
    [Test]
    public void BuildIndex_FromAllKnownMods_IndexesCorrectly() { }
    
    [Test]
    public void LookupByResidId_ExistingMod_FindsCorrectly() { }
    
    [Test]
    public void LookupByPsiModId_ExistingMod_FindsCorrectly() { }
    
    [Test]
    public void LookupByUnimodId_ExistingMod_FindsCorrectly() { }
    
    [Test]
    public void TryGetEquivalent_MetaMorpheusToUniProt_ConvertsCorrectly() { }
    
    [Test]
    public void TryGetEquivalent_NoEquivalent_ReturnsFalse() { }
}

[TestFixture]
public class SequenceConverterTests
{
    [Test]
    public void ConvertFullSequence_MetaMorpheusToUniProt_ConvertsCorrectly() { }
    
    [Test]
    public void ConvertFullSequence_UniProtToUnimod_ConvertsCorrectly() { }
    
    [Test]
    public void ConvertFullSequence_WithTerminalMods_ConvertsCorrectly() { }
    
    [Test]
    public void ConvertFullSequence_UnknownMod_ThrowsOnThrowMode() { }
    
    [Test]
    public void ConvertFullSequence_UnknownMod_UsesMassShiftOnMassShiftMode() { }
    
    [Test]
    public void ToMassShiftNotation_StandardMod_FormatsCorrectly() { }
    
    [Test]
    public void ToMassShiftNotation_NegativeMass_IncludesSign() { }
    
    [Test]
    public void FromMassShiftNotation_ValidMass_FindsClosestMod() { }
    
    [Test]
    [TestCase("Phosphorylation on S", ModificationNamingConvention.MetaMorpheus)]
    [TestCase("Phosphoserine", ModificationNamingConvention.UniProt)]
    public void RoundTrip_KnownMod_Succeeds(string modName, ModificationNamingConvention source) { }
}

[TestFixture]
public class SequenceConverterIntegrationTests
{
    [Test]
    public void ConvertPeptide_FromPsmTsv_ToUniProtFormat() { }
    
    [Test]
    public void ConvertOligo_FromMetaMorpheus_ToMassShift() { }
}
```

## Key Test Cases

### Cross-Reference Index Tests
1. **Phosphorylation mapping**: MetaMorpheus "Phosphorylation on S" <-> UniProt "Phosphoserine" via RESID
2. **Carbamidomethyl mapping**: MetaMorpheus <-> Unimod via Unimod ID
3. **Multiple matches**: When multiple mods share same RESID, pick correct by target residue
4. **No cross-ref**: Mod with no database references returns false

### Sequence Converter Tests
1. **Simple conversion**: `PEPT[Common Fixed:Carbamidomethyl on C]IDE` -> UniProt format
2. **Terminal mods**: N-terminal acetylation conversion
3. **Multiple mods**: Sequence with phosphorylation + oxidation
4. **Mass shift**: `PEPTIDE` with 79.966 mass mod -> `PEPT[+79.9663]IDE`
5. **Failure modes**: Verify each `ConversionFailureMode` works correctly

### Integration Tests
1. Create `PeptideWithSetModifications`, convert to UniProt, verify round-trip
2. Read PSM from TSV, convert sequence, verify mods match
3. Large-scale: Convert 1000 peptides, verify performance

## Acceptance Criteria
- [ ] All test cases pass
- [ ] Code coverage >90% for SequenceConverter and ModificationCrossRefIndex
- [ ] Edge cases are tested
- [ ] Integration tests verify real-world scenarios
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Run all new tests
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~SequenceConverter"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~CrossRefIndex"

# Run with coverage (if set up)
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --collect:"XPlat Code Coverage"
```
