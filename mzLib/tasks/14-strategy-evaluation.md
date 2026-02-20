# Task 14: Combinatorial Strategy Evaluation

## Objective
Systematically test all meaningful combinations of normalization, roll-up, and collapse strategies against the real spike-in data. Produce a ranked comparison table showing which combination best recovers the known fold changes.

## Prerequisites
- Tasks 9-13 must all be complete (baseline metrics + all new strategies implemented)

## What to Do

### 1. Create a comprehensive evaluation test

Add a test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void StrategyEvaluation_AllCombinations_UPS()
{
    string dataDir = Path.Combine(GetSolutionDir(), "TMT_Spike-In_Info");
    string upsPsmPath = Path.Combine(dataDir, "UPS_TMT3_Search", "Task1-SearchTask", "AllPSMs.psmtsv");

    // Define strategy options
    var psmNorms = new INormalizationStrategy[]
    {
        new NoNormalization(),
        new GlobalMedianNormalization(),
        new ReferenceChannelNormalization()
    };

    var rollUps = new IRollUpStrategy[]
    {
        new SumRollUp(),
        new MedianRollUp()
    };

    var peptideNorms = new INormalizationStrategy[]
    {
        new NoNormalization(),
        new GlobalMedianNormalization()
    };

    var collapses = new ICollapseStrategy[]
    {
        new NoCollapse(),
        new MeanCollapse()
    };

    var proteinNorms = new INormalizationStrategy[]
    {
        new NoNormalization(),
        new GlobalMedianNormalization()
    };

    // Evaluate each combination
    var results = new List<(string config, double mae_8x, double mae_2x, double medianFc_8x, double medianFc_2x)>();

    foreach (var psmNorm in psmNorms)
    foreach (var rollUp in rollUps)
    foreach (var pepNorm in peptideNorms)
    foreach (var collapse in collapses)
    foreach (var protNorm in proteinNorms)
    {
        string config = $"{psmNorm.Name} | {rollUp.Name} | {pepNorm.Name} | {collapse.Name} | {rollUp.Name} | {protNorm.Name}";

        try
        {
            var proteinMatrix = RunPipelineWithStrategies(
                upsPsmPath, psmNorm, rollUp, pepNorm, collapse, rollUp, protNorm);

            var fc_8x = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.125")
                .Select(x => x.foldChange).ToList();
            var fc_2x = QuantificationEvaluator.ComputeFoldChanges(proteinMatrix, "1", "0.5")
                .Select(x => x.foldChange).ToList();

            if (fc_8x.Count > 0 && fc_2x.Count > 0)
            {
                double mae_8x = QuantificationEvaluator.MeanAbsoluteLog2Error(fc_8x, 8.0);
                double mae_2x = QuantificationEvaluator.MeanAbsoluteLog2Error(fc_2x, 2.0);
                double medFc_8x = QuantificationEvaluator.Median(fc_8x);
                double medFc_2x = QuantificationEvaluator.Median(fc_2x);

                results.Add((config, mae_8x, mae_2x, medFc_8x, medFc_2x));
            }
        }
        catch (Exception ex)
        {
            TestContext.WriteLine($"FAILED: {config} -> {ex.Message}");
        }
    }

    // Sort by combined MAE (sum of log2 errors for both comparisons)
    results = results.OrderBy(r => r.mae_8x + r.mae_2x).ToList();

    // Print ranked results
    TestContext.WriteLine("=== Strategy Evaluation Results (sorted by combined MAE_log2) ===");
    TestContext.WriteLine($"{"Rank",-5} {"MAE_log2(8x)",-14} {"MAE_log2(2x)",-14} {"MedFC(8x)",-12} {"MedFC(2x)",-12} {"Config"}");
    TestContext.WriteLine(new string('-', 120));

    for (int i = 0; i < results.Count; i++)
    {
        var r = results[i];
        TestContext.WriteLine($"{i + 1,-5} {r.mae_8x,-14:F3} {r.mae_2x,-14:F3} {r.medianFc_8x,-12:F3} {r.medianFc_2x,-12:F3} {r.config}");
    }

    // Assert that at least some combinations were evaluated
    Assert.That(results.Count, Is.GreaterThan(0), "No strategy combinations could be evaluated");

    // Assert that the best combination is better than the worst
    // (sanity check that strategies actually make a difference)
    if (results.Count > 1)
    {
        double bestScore = results.First().mae_8x + results.First().mae_2x;
        double worstScore = results.Last().mae_8x + results.Last().mae_2x;
        TestContext.WriteLine($"\nBest combined MAE_log2: {bestScore:F3}");
        TestContext.WriteLine($"Worst combined MAE_log2: {worstScore:F3}");
    }
}
```

### 2. Practical considerations

- **Runtime**: Each pipeline run takes ~38s for UPS data. With ~48 combinations (3 x 2 x 2 x 2 x 2), expect ~30 minutes total.
  - To reduce runtime, pre-load the PSMs once and reuse them across combinations. Modify `RunPipelineWithStrategies` to accept pre-loaded data:
    ```csharp
    private static QuantMatrix<IBioPolymerGroup> RunPipelineWithStrategies(
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups,
        INormalizationStrategy psmNorm, ...)
    ```
  - Load PSMs once at the start of the test, then pass them to each combination.

- **Note on data mutation**: Check whether `RunTmtAndReturnProteinMatrix()` mutates the input lists. If it does, you'll need to clone them for each run. If not (likely since matrices are created fresh), reuse is safe.

### 3. Optional: HeLa background evaluation

If time permits, add an explicit test that evaluates the best few UPS configurations against HeLa data to check for false positives:

```csharp
[Test, Explicit("Long-running HeLa evaluation")]
public void StrategyEvaluation_TopConfigs_HeLa()
{
    // Run the top 3-5 configurations from the UPS evaluation against HeLa data
    // For each: compute median fold change and fraction within 1.5x of 1.0
    // The best overall strategy should be good on BOTH UPS (high accuracy) and HeLa (low false positives)
}
```

## Output

The test output should clearly show:
1. A ranked table of all strategy combinations sorted by accuracy
2. The winning combination
3. How much improvement the best combination provides over the baseline (NoNormalization everywhere)

## Verification

1. `dotnet build "C:/Users/Alex/Source/Repos/mzLib/mzLib/mzLib.sln"` succeeds
2. The evaluation test runs to completion and produces a ranked table
3. No regression in existing tests
4. Commit the results and update Activity.md with the winning strategy combination
