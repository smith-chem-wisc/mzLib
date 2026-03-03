# Task: Consolidate Parsing Methods in SpectrumMatchFromTsv

## Objective
Refactor the generic parsing methods in `SpectrumMatchFromTsv` to support both nullable and non-nullable types from a single unified method, eliminating code duplication and improving maintainability.

## Background
Currently, we have separate methods for parsing optional values:
- `GetOptionalValue<TNumber>` - returns non-nullable with default
- `GetOptionalValueNullable<TNumber>` - returns nullable

This creates confusion and duplicate code. We need a single method that handles both cases.

## Implementation Steps
- [x] Create unified `GetOptionalValue<TNumber>` that returns nullable
- [x] Update `SpectrumMatchFromTsv` constructor to use new method
- [x] Update `PsmFromTsv` constructor to use new method
- [x] Update `GlycoPsmFromTsv` constructor to use new method
- [x] Add `GetOptionalChemicalFormula` helper method
- [x] Update `OsmFromTsv` constructor to use new method
- [x] Remove duplicate code and old methods
- [x] Build and verify all changes compile

## Acceptance Criteria
- [x] Single `GetOptionalValue<TNumber>` method returns `TNumber?`
- [x] All derived classes use the unified method
- [x] No duplicate parsing code remains
- [x] Solution builds without errors
- [x] All existing tests pass
- [x] Changes are committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test Readers specifically
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~Readers"

# Run all tests
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj
```

## Files Modified
- `Readers/InternalResults/IndividualResultRecords/SpectrumMatchFromTsv.cs`
- `Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`
- `Readers/InternalResults/IndividualResultRecords/GlycoPsmFromTsv.cs`
- `Readers/InternalResults/IndividualResultRecords/OsmFromTsv.cs`

## Implementation Notes
- Use `.GetValueOrDefault()` when assigning nullable return to non-nullable property
- For properties that should be nullable (like `PrecursorIntensity`, `TotalIonCurrent`), use the nullable return directly
- The `GetOptionalChemicalFormula` method wraps ChemicalFormula parsing with null handling
- All parsing follows the same pattern: check header exists ? check value is not empty ? try parse ? return default on failure

## Status
? COMPLETED
- All parsing methods consolidated
- Build successful
- Tests pass
- Changes committed
