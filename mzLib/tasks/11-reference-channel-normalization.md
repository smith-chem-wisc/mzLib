# Task 11: Implement ReferenceChannelNormalization Strategy

## Objective
Implement a normalization strategy that uses the reference channels (126 and 131N) to normalize each non-reference channel within each file. This converts raw intensities to ratios relative to the reference, removing file-to-file variation.

## Background

From the MSstatsTMT paper (Table I): "Local ratio-based normalization" computes log2 ratios of each channel to the reference channel within each spectrum, then normalizes at the protein level by equalizing reference channel summaries across runs.

At the PSM/spectrum level, the simplest form is:
- For each row (PSM/peptide/protein), divide each non-reference channel by the average of the two reference channels in that row
- This converts absolute intensities to fold-change-relative-to-reference

**Algorithm:**
1. Identify reference columns (where `ISampleInfo` is `IsobaricQuantSampleInfo` with `IsReferenceChannel == true`)
2. Group columns by file (using the file name from `ISampleInfo`)
3. For each row:
   a. For each file, compute the mean of the reference channel values for that row in that file
   b. Divide each non-reference channel value by its file's reference mean
   c. Set reference channel values to 1.0 (or leave as the ratio, which is 1.0 by construction)
4. If both reference channels are zero for a file/row, leave all values for that file as zero (missing data)

## File to Create

### `Quantification/Strategies/Normalization/ReferenceChannelNormalization.cs`

```csharp
using MassSpectrometry.ExperimentalDesign;
using Quantification.Interfaces;

namespace Quantification.Strategies;

public class ReferenceChannelNormalization : INormalizationStrategy
{
    public string Name => "Reference Channel Normalization";

    public QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
    {
        // 1. Identify reference and non-reference column indices
        //    Reference columns: IsobaricQuantSampleInfo.IsReferenceChannel == true
        // 2. Group reference columns by file name
        // 3. For each row:
        //    - For each file group, compute mean of reference channel values
        //    - Divide each non-reference column by its file's reference mean
        //    - If reference mean is 0, set all columns for that file to 0
        // 4. Return new matrix with ratio values
    }
}
```

## Key Design Decisions

- **Per-file reference**: Each file has its own pair of reference channels (126, 131N). The normalization is done within each file independently.
- **Two reference channels**: Average both reference values for a more robust denominator. If one is zero and the other is non-zero, use only the non-zero one.
- **Output scale**: Values become ratios (fold change relative to reference). A value of 2.0 means "twice the reference intensity."
- **Zero handling**: If both reference channels are zero for a given row in a given file, all channels for that file in that row should be zero (can't normalize without a reference).
- **Works at any level**: This strategy works on PSM matrices (per-file, step 2), peptide matrices (combined, step 5), or protein matrices (step 8), but is most meaningful at steps 2 or 5 where per-file structure is preserved.

## Unit Test

Add a test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void ReferenceChannelNormalization_ConvertsToRatios()
{
    // Create a small matrix where:
    //   - Column 0 is reference (value = 100 for row 0, 200 for row 1)
    //   - Column 1 is non-reference (value = 200, 400)
    //   - Column 2 is reference (value = 100, 200)
    //
    // After normalization:
    //   - Reference columns become 1.0
    //   - Column 1 row 0: 200 / mean(100,100) = 2.0
    //   - Column 1 row 1: 400 / mean(200,200) = 2.0
    //
    // Assert values match expected ratios within tolerance
}
```

## Verification

1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` succeeds
2. New unit test passes
3. No regression in existing tests
