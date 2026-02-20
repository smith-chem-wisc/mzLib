# Task 12: Implement MedianRollUp Strategy

## Objective
Implement a roll-up strategy that uses the median (instead of sum) to aggregate lower-level entities into higher-level entities. Median is more robust to outliers than sum and is commonly used in proteomics quantification.

## Background

When rolling up PSMs to peptides or peptides to proteins, using the median instead of the sum:
- Is more robust to outlier PSMs/peptides with extreme intensities
- Better represents the "typical" intensity for a protein
- Is less sensitive to the number of PSMs identified (sum is biased toward proteins with more PSMs)

## File to Create

### `Quantification/Strategies/RollUp/MedianRollUp.cs`

```csharp
using Quantification.Interfaces;

namespace Quantification.Strategies;

public class MedianRollUp : IRollUpStrategy
{
    public string Name => "Median Roll-Up";

    public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map)
        where TLow : IEquatable<TLow>
        where THigh : IEquatable<THigh>
    {
        // Same structure as SumRollUp, but instead of summing values across rows,
        // take the median of non-zero values for each column position.
        //
        // For each higher-level key (e.g., protein):
        //   For each column (sample/channel):
        //     Collect all non-zero values from the mapped lower-level rows
        //     Result = median of those values (or 0 if no non-zero values)
    }
}
```

## Key Design Decisions

- **Zero handling**: Zeros represent missing data. When computing the median for a column, only consider non-zero values. If all values are zero, the result is zero.
- **Single value**: If only one non-zero value exists, the median is that value.
- **Even count**: For an even number of values, use the average of the two middle values (standard median definition).

## Unit Test

Add a test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void MedianRollUp_ComputesMedianPerColumn()
{
    // Create a matrix with 3 PSM rows mapping to 1 protein, 2 columns:
    //   PSM1: [100, 200]
    //   PSM2: [300, 400]
    //   PSM3: [500, 100]
    //
    // After MedianRollUp:
    //   Protein: [300, 200]  (median of [100,300,500] = 300; median of [200,400,100] = 200)
    //
    // Also test with zeros:
    //   PSM1: [100, 0]
    //   PSM2: [300, 400]
    //   PSM3: [0, 100]
    //
    // After MedianRollUp (ignoring zeros):
    //   Protein: [200, 250]  (median of [100,300] = 200; median of [400,100] = 250)
}
```

## Verification

1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` succeeds
2. New unit test passes
3. No regression in existing tests
