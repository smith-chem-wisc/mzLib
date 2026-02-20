# Task 9: Establish Baseline Accuracy Metrics

## Objective
Run the existing NoNormalization + SumRollUp + NoCollapse pipeline on both UPS and HeLa data, and record precise baseline metrics. These numbers will be the benchmark against which all new strategies are compared.

## What to Do

### 1. Extend QuantificationEvaluator with additional metrics

Add the following methods to `Test/Quantification/TestHelpers/QuantificationEvaluator.cs`:

```csharp
/// <summary>
/// Computes the Mean Absolute Error (MAE) between observed and expected fold changes.
/// </summary>
public static double MeanAbsoluteError(List<double> observed, double expected)
{
    return observed.Average(fc => Math.Abs(fc - expected));
}

/// <summary>
/// Computes the Mean Absolute Log2 Error between observed and expected fold changes.
/// Log2 space is symmetric: a 2x overestimate and 2x underestimate are equidistant from truth.
/// </summary>
public static double MeanAbsoluteLog2Error(List<double> observed, double expected)
{
    double expectedLog2 = Math.Log2(expected);
    return observed.Average(fc => Math.Abs(Math.Log2(fc) - expectedLog2));
}

/// <summary>
/// Computes what fraction of observed fold changes are within a tolerance factor of the expected value.
/// E.g., tolerance=2.0 means within [expected/2, expected*2].
/// </summary>
public static double FractionWithinFactor(List<double> observed, double expected, double factor)
{
    return (double)observed.Count(fc => fc >= expected / factor && fc <= expected * factor) / observed.Count;
}
```

### 2. Create a baseline metrics test

Add a new test to `Test/Quantification/TmtSpikeInTests.cs`:

```csharp
[Test]
public void BaselineMetrics_NoNormalization_UPS()
{
    // Load UPS data, run pipeline with NoNormalization + SumRollUp + NoCollapse
    // Compute and print (via TestContext.WriteLine) the following metrics:
    //
    // For each fold change comparison (1 vs 0.125, 1 vs 0.5, 1 vs 0.667):
    //   - Number of quantifiable proteins
    //   - Median observed fold change
    //   - Mean observed fold change
    //   - Mean Absolute Error (MAE) vs expected
    //   - Mean Absolute Log2 Error vs expected
    //   - Fraction within 2x of expected
    //
    // Do NOT assert specific values -- just print them.
    // The test should pass as long as it runs without errors.
}
```

### 3. Create a baseline metrics test for HeLa background

```csharp
[Test, Explicit("Large HeLa file")]
public void BaselineMetrics_NoNormalization_HeLa()
{
    // Load HeLa data, run pipeline with NoNormalization + SumRollUp + NoCollapse
    // Compute and print:
    //   - Number of quantifiable proteins
    //   - Median fold change "1" vs "0.5" (expected: 1.0)
    //   - Mean Absolute Error vs 1.0
    //   - Fraction within 1.5x of 1.0 (i.e., between 0.667 and 1.5)
    //   - Fraction within 2x of 1.0 (i.e., between 0.5 and 2.0)
}
```

### 4. Helper method for running a strategy combination

To avoid duplicating the load-configure-run pattern across Tasks 9 and 14, add a shared helper method to `TmtSpikeInTests.cs`:

```csharp
/// <summary>
/// Loads PSMs from a psmtsv file, configures a QuantificationEngine with the given strategies,
/// runs the TMT pipeline, and returns the protein matrix.
/// </summary>
private static QuantMatrix<IBioPolymerGroup> RunPipelineWithStrategies(
    string psmtsvPath,
    INormalizationStrategy psmNorm,
    IRollUpStrategy psmToPeptideRollUp,
    INormalizationStrategy peptideNorm,
    ICollapseStrategy collapse,
    IRollUpStrategy peptideToProteinRollUp,
    INormalizationStrategy proteinNorm)
{
    var (spectralMatches, peptides, proteinGroups) =
        PsmTsvQuantAdapter.BuildQuantificationInputs(psmtsvPath);

    var design = new SpikeInExperimentalDesign();
    var parameters = new QuantificationParameters
    {
        SpectralMatchNormalizationStrategy = psmNorm,
        SpectralMatchToPeptideRollUpStrategy = psmToPeptideRollUp,
        PeptideNormalizationStrategy = peptideNorm,
        CollapseStrategy = collapse,
        PeptideToProteinRollUpStrategy = peptideToProteinRollUp,
        ProteinNormalizationStrategy = proteinNorm,
        OutputDirectory = string.Empty,
        WriteRawInformation = false,
        WritePeptideInformation = false,
        WriteProteinInformation = false
    };

    var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);
    return engine.RunTmtAndReturnProteinMatrix();
}
```

## Output Format

Print results using `TestContext.WriteLine` in a structured format like:

```
=== Baseline Metrics (NoNormalization + SumRollUp + NoCollapse) ===
Comparison: 1 vs 0.125 (expected FC = 8.0)
  Proteins quantified: 42
  Median FC: 3.21
  Mean FC: 4.15
  MAE: 4.79
  MAE(log2): 1.32
  Fraction within 2x: 0.38

Comparison: 1 vs 0.5 (expected FC = 2.0)
  Proteins quantified: 44
  Median FC: 1.45
  ...
```

## Verification

1. `dotnet test --filter "FullyQualifiedName~BaselineMetrics_NoNormalization_UPS"` passes
2. Check the test output for the printed metrics
3. No regression in existing tests
