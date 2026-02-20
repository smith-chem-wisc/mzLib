# Task 6: Add RunTmt Pipeline Method to QuantificationEngine

## Objective
Add `RunTmt()` and `RunTmtAndReturnProteinMatrix()` methods that implement the full TMT-specific quantification pipeline using per-file processing.

## Prerequisites
- Task 4 (PivotByFile) and Task 5 (CombinePeptideMatrices) must be complete

## Background
The TMT pipeline differs from the LFQ pipeline because normalization must happen within each file before combining across files. The steps are:

1. **PivotByFile** → one SpectralMatchMatrix per file
2. **Per-file PSM normalization** → normalize within each file
3. **Per-file PSM-to-peptide roll-up** → peptide matrix per file
4. **Combine** → single peptide matrix across all files
5. **Peptide normalization** → normalize combined peptide matrix
6. **Collapse** → combine technical replicates
7. **Peptide-to-protein roll-up** → protein matrix
8. **Protein normalization** → normalize protein matrix

## File to Modify

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Quantification/QuantificationEngine.cs`

**Do NOT modify** the existing `Run()` or `RunAndReturnProteinMatrix()` methods.

## Methods to Add

### RunTmt()
```csharp
/// <summary>
/// Runs the TMT-specific quantification pipeline.
/// Processes each file independently (pivot, normalize, roll-up) before combining.
/// </summary>
public QuantificationResults RunTmt()
{
    // 0) Validate engine state
    if (!ValidateEngine(out QuantificationResults badResults))
    {
        return badResults;
    }

    // 1) PivotByFile - one matrix per file
    var perFileMatrices = PivotByFile(SpectralMatches, ExperimentalDesign);

    // 2) Per-file PSM normalization
    var perFileNormalized = new Dictionary<string, SpectralMatchMatrix>();
    foreach (var kvp in perFileMatrices)
    {
        perFileNormalized[kvp.Key] = (SpectralMatchMatrix)Parameters
            .SpectralMatchNormalizationStrategy.NormalizeIntensities(kvp.Value);
    }

    // 3) Per-file roll-up to peptides
    var perFilePeptideMatrices = new Dictionary<string, QuantMatrix<IBioPolymerWithSetMods>>();
    foreach (var kvp in perFileNormalized)
    {
        var peptideMap = GetPsmToPeptideMap(kvp.Value, ModifiedBioPolymers);
        perFilePeptideMatrices[kvp.Key] = Parameters
            .SpectralMatchToPeptideRollUpStrategy.RollUp(kvp.Value, peptideMap);
    }

    // 4) Combine per-file peptide matrices
    var combinedPeptideMatrix = CombinePeptideMatrices(perFilePeptideMatrices, ExperimentalDesign);

    // 5) Normalize combined peptide matrix
    var peptideMatrixNorm = Parameters.PeptideNormalizationStrategy
        .NormalizeIntensities(combinedPeptideMatrix);

    // 6) Collapse samples (technical replicates, fractions)
    peptideMatrixNorm = Parameters.CollapseStrategy.CollapseSamples(peptideMatrixNorm);

    // 7) Roll up to proteins
    var proteinMap = Parameters.UseSharedPeptidesForProteinQuant
        ? GetAllPeptideToProteinMap(peptideMatrixNorm)
        : GetUniquePeptideToProteinMap(peptideMatrixNorm, BioPolymerGroups);
    var proteinMatrix = Parameters.PeptideToProteinRollUpStrategy
        .RollUp(peptideMatrixNorm, proteinMap);

    // 8) Normalize protein matrix
    var proteinMatrixNorm = Parameters.ProteinNormalizationStrategy
        .NormalizeIntensities(proteinMatrix);

    return new QuantificationResults
    {
        Summary = "TMT Quantification completed successfully."
    };
}
```

### RunTmtAndReturnProteinMatrix()
```csharp
/// <summary>
/// Runs the TMT pipeline and returns the final protein matrix for testing/inspection.
/// </summary>
internal QuantMatrix<IBioPolymerGroup> RunTmtAndReturnProteinMatrix()
{
    // Same as RunTmt() but returns proteinMatrixNorm instead of QuantificationResults
    // (Copy the logic, return the matrix instead of wrapping in results)
}
```

## Implementation Notes

### Cast Issue
The `INormalizationStrategy.NormalizeIntensities<T>()` returns `QuantMatrix<T>`. When `T = ISpectralMatch`, it returns `QuantMatrix<ISpectralMatch>`, not `SpectralMatchMatrix`. The cast to `SpectralMatchMatrix` may fail.

**Solution**: Either:
- Don't cast - use `QuantMatrix<ISpectralMatch>` throughout the per-file section
- Or modify the code to work with the base type

The `GetPsmToPeptideMap` accepts `QuantMatrix<ISpectralMatch>`, so no cast needed.

### Peptide Map Across Files
`GetPsmToPeptideMap()` maps PSMs to peptides from the `ModifiedBioPolymers` list. It iterates over the matrix's RowKeys (PSMs) and finds which peptide each PSM maps to. This works correctly per-file because:
- The peptide list is shared across all files
- Each file's PSMs will only map to peptides that appear in that file
- Peptides not found in a file will have empty lists in the map

The roll-up strategy will then create rows only for peptides that had mapped PSMs.

### DRY Principle
`RunTmt()` and `RunTmtAndReturnProteinMatrix()` share almost identical code. Consider extracting a private helper that returns the protein matrix, then:
- `RunTmt()` calls the helper and wraps in QuantificationResults
- `RunTmtAndReturnProteinMatrix()` calls the helper and returns the matrix

## Verification
1. Build compiles
2. Existing tests pass
3. Will be fully tested in Task 7
