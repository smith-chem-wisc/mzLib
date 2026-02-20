# Task 13: Implement MeanCollapse Strategy

## Objective
Implement a collapse strategy that averages (instead of sums) technical replicates. This produces values on the same scale as individual measurements, making fold change calculations more interpretable.

## Background

The existing `SumCollapse` sums intensities across technical replicates and fractions grouped by `Condition_BiologicalReplicate`. This inflates values proportionally to the number of replicates, which doesn't affect fold change ratios but can make absolute values harder to interpret.

`MeanCollapse` averages instead of sums, producing values that represent the mean intensity per condition per biological replicate.

## File to Create

### `Quantification/Strategies/Collapse/MeanCollapse.cs`

```csharp
using Quantification.Interfaces;

namespace Quantification.Strategies;

public class MeanCollapse : ICollapseStrategy
{
    public string Name => "Mean Collapse";

    public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
    {
        // Same grouping logic as SumCollapse (group by Condition + '_' + BiologicalReplicate)
        // But divide the sum by the number of non-zero values in the group (or total count)
        //
        // Design decision: divide by number of members in the group (not number of non-zero values)
        // This is consistent with treating zeros as "observed zero" rather than "missing"
        // For TMT data, a zero in one tech rep is meaningful (the peptide was identified but had zero reporter ion intensity)
    }
}
```

## Key Design Decisions

- **Denominator**: Divide by total group size (number of technical replicates), NOT by number of non-zero values. If a protein has intensity [100, 0, 200] across 3 tech reps, the mean is 100, not 150. This avoids inflating estimates when some replicates have missing data.
- **Follow SumCollapse pattern**: Use the exact same grouping and column selection logic as `SumCollapse`. The only difference is dividing by the group count at the end.

## Unit Test

Add a test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void MeanCollapse_AveragesTechnicalReplicates()
{
    // Same setup as RunTmt_WithSumCollapse_CombinesTechnicalReplicates:
    // 3 files (tech reps), 3 channels, 1 protein with values [100, 200, 300] per file
    //
    // After MeanCollapse:
    //   Reference_1: mean(100, 100, 100) = 100  (SumCollapse gives 300)
    //   CondA_1:     mean(200, 200, 200) = 200  (SumCollapse gives 600)
    //   CondB_1:     mean(300, 300, 300) = 300  (SumCollapse gives 900)
}
```

## Verification

1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` succeeds
2. New unit test passes
3. No regression in existing tests
