using BenchmarkDotNet.Attributes;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.SSRCalc;

namespace Benchmarks.Chromatography;

/// <summary>
/// Benchmarks for retention time predictors.
/// Compares SSRCalc3 (algorithmic) vs Chronologer (deep learning) performance.
/// </summary>
[Config(typeof(MzLibBenchmarkConfig))]
[MemoryDiagnoser]
public class RetentionTimePredictorBenchmarks : IDisposable
{
    private SSRCalc3RetentionTimePredictor _ssrCalc3 = null!;
    private ChronologerRetentionTimePredictor _chronologer = null!;
    private List<IRetentionPredictable> _peptides = null!;

    [Params(10, 100, 1000)]
    public int PeptideCount { get; set; }

    [GlobalSetup]
    public void Setup()
    {
        _ssrCalc3 = new SSRCalc3RetentionTimePredictor();
        _chronologer = new ChronologerRetentionTimePredictor();
        _peptides = TestDataGenerator.GenerateUnmodifiedPeptides(PeptideCount);
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        _chronologer?.Dispose();
    }

    [Benchmark(Baseline = true)]
    public int SSRCalc3_PredictAll()
    {
        int successCount = 0;
        foreach (var peptide in _peptides)
        {
            var result = _ssrCalc3.PredictRetentionTime(peptide, out _);
            if (result.HasValue) successCount++;
        }
        return successCount;
    }

    [Benchmark]
    public int Chronologer_PredictAll()
    {
        int successCount = 0;
        foreach (var peptide in _peptides)
        {
            var result = _chronologer.PredictRetentionTime(peptide, out _);
            if (result.HasValue) successCount++;
        }
        return successCount;
    }

    public void Dispose()
    {
        Cleanup();
        GC.SuppressFinalize(this);
    }
}

/// <summary>
/// Benchmarks for single peptide prediction latency.
/// Useful for understanding per-call overhead.
/// </summary>
[Config(typeof(MzLibBenchmarkConfig))]
[MemoryDiagnoser]
public class RetentionTimePredictorSinglePeptideBenchmarks : IDisposable
{
    private SSRCalc3RetentionTimePredictor _ssrCalc3 = null!;
    private ChronologerRetentionTimePredictor _chronologer = null!;
    private IRetentionPredictable _shortPeptide = null!;
    private IRetentionPredictable _longPeptide = null!;

    [GlobalSetup]
    public void Setup()
    {
        _ssrCalc3 = new SSRCalc3RetentionTimePredictor();
        _chronologer = new ChronologerRetentionTimePredictor();

        var peptides = TestDataGenerator.GeneratePeptidesWithVaryingLengths(7, 7, 1);
        _shortPeptide = peptides[0];

        var longPeptides = TestDataGenerator.GeneratePeptidesWithVaryingLengths(40, 40, 1);
        _longPeptide = longPeptides[0];
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        _chronologer?.Dispose();
    }

    [Benchmark(Baseline = true)]
    public double? SSRCalc3_ShortPeptide() => _ssrCalc3.PredictRetentionTime(_shortPeptide, out _);

    [Benchmark]
    public double? SSRCalc3_LongPeptide() => _ssrCalc3.PredictRetentionTime(_longPeptide, out _);

    [Benchmark]
    public double? Chronologer_ShortPeptide() => _chronologer.PredictRetentionTime(_shortPeptide, out _);

    [Benchmark]
    public double? Chronologer_LongPeptide() => _chronologer.PredictRetentionTime(_longPeptide, out _);

    public void Dispose()
    {
        Cleanup();
        GC.SuppressFinalize(this);
    }
}