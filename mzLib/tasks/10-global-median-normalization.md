# Task 10: Implement GlobalMedianNormalization Strategy

## Objective
Implement a global median normalization strategy that equalizes the median log2 intensity across all columns (channels/samples) in a QuantMatrix. This is the spectrum-level normalization used by MSstatsTMT.

## Background

From the MSstatsTMT paper: "Global median normalization equalizes the median log2 intensity across all spectra, channels, and runs." The goal is to remove systematic biases between channels within each file and across files.

**Algorithm:**
1. Log2-transform all positive intensities in the matrix
2. Compute the median log2 intensity for each column (channel)
3. Compute the global median across all column medians
4. For each column, compute shift = global_median - column_median
5. Add the shift to all log2 values in that column
6. Transform back: intensity = 2^(shifted_log2_value)
7. Zero values remain zero (they are not log-transformed or shifted)

## File to Create

### `Quantification/Strategies/Normalization/GlobalMedianNormalization.cs`

```csharp
using Quantification.Interfaces;

namespace Quantification.Strategies;

public class GlobalMedianNormalization : INormalizationStrategy
{
    public string Name => "Global Median Normalization";

    public QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
    {
        // 1. Clone the matrix data
        // 2. For each column, collect all positive (non-zero) values, log2-transform, compute median
        // 3. Compute global median of all column medians
        // 4. For each column, shift = globalMedian - columnMedian
        // 5. Apply shift: for each positive value in column, newVal = 2^(log2(oldVal) + shift)
        // 6. Return new QuantMatrix with shifted values
    }
}
```

## Key Design Decisions

- **Zero handling**: Zeros represent missing data (PSM not observed in that channel). Do NOT include zeros when computing column medians. Do NOT modify zero values.
- **Negative values**: Should not occur in raw intensity data. If encountered, treat as missing (skip).
- **Single-value columns**: If a column has only one non-zero value, skip normalization for that column (its median equals itself, no meaningful shift).
- **Empty columns**: If a column has no non-zero values, skip it entirely.

## Unit Test

Add a test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void GlobalMedianNormalization_EqualizesColumnMedians()
{
    // Create a small matrix (2 rows x 3 columns) with known values
    // Column medians before: e.g., [100, 200, 400] -> log2 = [6.64, 7.64, 8.64]
    // Global median of column medians = 7.64
    // After normalization, all column medians should be equal (or very close)

    // Assert that column medians are within 1% of each other after normalization
    // Assert that zeros remain zeros
}
```

## Verification

1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` succeeds
2. New unit test passes
3. No regression in existing tests
