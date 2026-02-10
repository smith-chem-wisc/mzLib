using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Diagnosers;
using BenchmarkDotNet.Exporters;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Reports;

namespace Benchmarks;

/// <summary>
/// Shared benchmark configuration for all mzLib benchmarks.
/// </summary>
public class MzLibBenchmarkConfig : ManualConfig
{
    public MzLibBenchmarkConfig()
    {
        // Use short run for quick iteration during development
        // Switch to Job.Default for final measurements
        AddJob(Job.ShortRun
            .WithWarmupCount(3)
            .WithIterationCount(5));

        // Memory allocation diagnostics
        AddDiagnoser(MemoryDiagnoser.Default);

        // Export results to multiple formats
        AddExporter(MarkdownExporter.GitHub);

        // Summary style
        WithSummaryStyle(SummaryStyle.Default.WithRatioStyle(BenchmarkDotNet.Columns.RatioStyle.Trend));
    }
}

/// <summary>
/// Longer-running configuration for comprehensive benchmarks.
/// </summary>
public class MzLibBenchmarkConfigFull : ManualConfig
{
    public MzLibBenchmarkConfigFull()
    {
        AddJob(Job.Default);
        AddDiagnoser(MemoryDiagnoser.Default);
        AddExporter(MarkdownExporter.GitHub);
        WithSummaryStyle(SummaryStyle.Default.WithRatioStyle(BenchmarkDotNet.Columns.RatioStyle.Trend));
    }
}