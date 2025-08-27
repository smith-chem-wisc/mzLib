using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace MassSpectrometry.FlashDeconvRuntime;

/// <summary>
/// Strongly-typed representation of FLASHDeconv command line options (full --helphelp set).
/// Only includes flags/values; does NOT execute the process.
/// Use ToArgumentList() to serialize to arguments list (excluding -in / -out which are passed separately).
/// </summary>
public sealed class FlashDeconvOptions
{
    // Mandatory (handled externally, but kept for convenience)
    public string? InputMzMlPath { get; set; }            // -in
    public string? OutputTsvPath { get; set; }            // -out

    // Optional direct output variants
    public List<string> OutSpectraPerMsLevelTsv { get; } = new();   // -out_spec (multiple file list)
    public string? OutMzmlDeconvolved { get; set; }                 // -out_mzml
    public string? OutAnnotatedMzml { get; set; }                   // -out_annotated_mzml
    public string? OutPromex { get; set; }                          // -out_promex
    public List<string> OutTopFDPerMsLevelMsalign { get; } = new(); // -out_topFD
    public List<string> OutTopFDFeaturePerMsLevel { get; } = new(); // -out_topFD_feature

    // Coupled FLASHIda log
    public string? FlashIdaLog { get; set; }                        // -in_log

    // Targeted precursor overrides
    public double TargetPrecursorSnrMin { get; set; } = 1.0;        // -min_precursor_snr
    public int TargetPrecursorCharge { get; set; } = 0;             // -target_precursor_charge
    public double TargetPrecursorMz { get; set; } = 0.0;            // -target_precursor_mz

    // Output mass/charge formatting
    public int MzmlMassChargeMode { get; set; } = 0;                // -mzml_mass_charge ( -1, 0, 1 )
    public int PrecedingMs1Count { get; set; } = 3;                 // -preceding_MS1_count
    public bool WriteDetail { get; set; } = false;                  // -write_detail (0/1)
    public int MaxMsLevel { get; set; } = 3;                        // -max_MS_level
    public int ForcedMsLevel { get; set; } = 0;                     // -forced_MS_level (0=off)
    public MergingMethodEnum MergingMethod { get; set; } = MergingMethodEnum.None; // -merging_method
    public bool ReportFdr { get; set; } = false;                    // -report_FDR
    public bool UseRnaAveragine { get; set; } = false;              // -use_RNA_averagine

    // Nested groups
    public AlgorithmGroup Algorithm { get; } = new();
    public FeatureTracingGroup FeatureTracing { get; } = new();
    public CommonTopGroup Common { get; } = new();

    #region Nested Types

    public enum MergingMethodEnum { None = 0, GaussianAveraging = 1, Block = 2 }

    public sealed class AlgorithmGroup
    {
        // Multi-valued tolerances (MS1, MS2, ...)
        public List<double> PpmTolerancePerMsLevel { get; } = new() { 10.0, 10.0, 10.0 }; // -Algorithm:tol
        public double MinMass { get; set; } = 50.0;                // -Algorithm:min_mass
        public double MaxMass { get; set; } = 1.0e5;               // -Algorithm:max_mass
        public int MinCharge { get; set; } = 1;                    // -Algorithm:min_charge
        public int MaxCharge { get; set; } = 100;                  // -Algorithm:max_charge
        public double MinMz { get; set; } = -1.0;                  // -Algorithm:min_mz (<=0 disables)
        public double MaxMz { get; set; } = -1.0;                  // -Algorithm:max_mz
        public double MinRt { get; set; } = -1.0;                  // -Algorithm:min_rt
        public double MaxRt { get; set; } = -1.0;                  // -Algorithm:max_rt
        public double IsolationWindowWidth { get; set; } = 5.0;    // -Algorithm:isolation_window
        public List<double> MinIsotopeCosinePerMsLevel { get; } = new() { 0.85, 0.85, 0.85 }; // -Algorithm:min_isotope_cosine
        public int AllowedIsotopeError { get; set; } = 1;          // -Algorithm:allowed_isotope_error
        public double MinIntensity { get; set; } = 0.0;            // -Algorithm:min_intensity
    }

    public sealed class FeatureTracingGroup
    {
        public double MassErrorPpm { get; set; } = -1.0;                 // -FeatureTracing:mass_error_ppm
        public double IonMobilityTolerance { get; set; } = 0.01;          // -FeatureTracing:ion_mobility_tolerance
        public QuantMethodEnum QuantMethod { get; set; } = QuantMethodEnum.area; // -FeatureTracing:quant_method
        public double MinSampleRate { get; set; } = 0.05;                 // -FeatureTracing:min_sample_rate
        public double MinTraceLengthSec { get; set; } = 10.0;             // -FeatureTracing:min_trace_length
        public double MaxTraceLengthSec { get; set; } = -1.0;             // -FeatureTracing:max_trace_length
        public double MinIsotopeCosineForFeature { get; set; } = -1.0;    // -FeatureTracing:min_isotope_cosine (-1 => use Algorithm)
        public enum QuantMethodEnum { area, median, max_height }
    }

    public sealed class CommonTopGroup
    {
        public string? IniFile { get; set; }                 // -ini
        public string? LogFile { get; set; }                 // -log
        public int Instance { get; set; } = 1;               // -instance
        public int DebugLevel { get; set; } = 0;             // -debug
        public int Threads { get; set; } = 1;                // -threads
        public string? WriteIni { get; set; }                // -write_ini
        public string? WriteCtdDir { get; set; }             // -write_ctd
        public string? WriteNestedCwlDir { get; set; }       // -write_nested_cwl
        public string? WriteCwlDir { get; set; }             // -write_cwl
        public string? WriteNestedJsonDir { get; set; }      // -write_nested_json
        public string? WriteJsonDir { get; set; }            // -write_json
        public bool NoProgress { get; set; } = false;        // -no_progress
        public bool Force { get; set; } = false;             // -force
        public bool TestMode { get; set; } = false;          // -test
    }

    #endregion

    #region Serialization

    /// <summary>
    /// Validates required fields; throws on missing mandatory arguments.
    /// </summary>
    public void Validate()
    {
        if (string.IsNullOrWhiteSpace(InputMzMlPath))
            throw new ArgumentException("InputMzMlPath (-in) is required.");
        if (string.IsNullOrWhiteSpace(OutputTsvPath))
            throw new ArgumentException("OutputTsvPath (-out) is required.");
    }

    /// <summary>
    /// Builds the full argument list (excluding the executable path).
    /// Includes -in and -out.
    /// </summary>
    public IReadOnlyList<string> ToArgumentList()
    {
        Validate();
        var args = new List<string>
        {
            "-in", Quote(InputMzMlPath!),
            "-out", Quote(OutputTsvPath!)
        };

        if (!string.IsNullOrWhiteSpace(FlashIdaLog))
            args.AddRange(new[] { "-in_log", Quote(FlashIdaLog!) });

        AddIfFiles(args, "-out_spec", OutSpectraPerMsLevelTsv);
        AddIfFiles(args, "-out_topFD", OutTopFDPerMsLevelMsalign);
        AddIfFiles(args, "-out_topFD_feature", OutTopFDFeaturePerMsLevel);

        AddIfFile(args, "-out_mzml", OutMzmlDeconvolved);
        AddIfFile(args, "-out_annotated_mzml", OutAnnotatedMzml);
        AddIfFile(args, "-out_promex", OutPromex);

        // Scalar simple
        args.AddRange(new[]
        {
            "-min_precursor_snr", F(TargetPrecursorSnrMin),
            "-target_precursor_charge", TargetPrecursorCharge.ToString(CultureInfo.InvariantCulture),
            "-target_precursor_mz", F(TargetPrecursorMz),
            "-mzml_mass_charge", MzmlMassChargeMode.ToString(CultureInfo.InvariantCulture),
            "-preceding_MS1_count", PrecedingMs1Count.ToString(CultureInfo.InvariantCulture),
            "-write_detail", (WriteDetail ? 1 : 0).ToString(),
            "-max_MS_level", MaxMsLevel.ToString(),
            "-forced_MS_level", ForcedMsLevel.ToString(),
            "-merging_method", ((int)MergingMethod).ToString(),
            "-report_FDR", (ReportFdr ? 1 : 0).ToString(),
            "-use_RNA_averagine", (UseRnaAveragine ? 1 : 0).ToString()
        });

        // Algorithm group
        if (Algorithm.PpmTolerancePerMsLevel.Count > 0)
            args.AddRange(new[] { "-Algorithm:tol" }.Concat(Algorithm.PpmTolerancePerMsLevel.Select(F)));
        args.AddRange(new[]
        {
            "-Algorithm:min_mass", F(Algorithm.MinMass),
            "-Algorithm:max_mass", F(Algorithm.MaxMass),
            "-Algorithm:min_charge", Algorithm.MinCharge.ToString(),
            "-Algorithm:max_charge", Algorithm.MaxCharge.ToString(),
            "-Algorithm:min_mz", F(Algorithm.MinMz),
            "-Algorithm:max_mz", F(Algorithm.MaxMz),
            "-Algorithm:min_rt", F(Algorithm.MinRt),
            "-Algorithm:max_rt", F(Algorithm.MaxRt),
            "-Algorithm:isolation_window", F(Algorithm.IsolationWindowWidth)
        });
        if (Algorithm.MinIsotopeCosinePerMsLevel.Count > 0)
            args.AddRange(new[] { "-Algorithm:min_isotope_cosine" }.Concat(Algorithm.MinIsotopeCosinePerMsLevel.Select(F)));
        args.AddRange(new[]
        {
            "-Algorithm:allowed_isotope_error", Algorithm.AllowedIsotopeError.ToString(),
            "-Algorithm:min_intensity", F(Algorithm.MinIntensity)
        });

        // FeatureTracing
        args.AddRange(new[]
        {
            "-FeatureTracing:mass_error_ppm", F(FeatureTracing.MassErrorPpm),
            "-FeatureTracing:ion_mobility_tolerance", F(FeatureTracing.IonMobilityTolerance),
            "-FeatureTracing:quant_method", FeatureTracing.QuantMethod.ToString(),
            "-FeatureTracing:min_sample_rate", F(FeatureTracing.MinSampleRate),
            "-FeatureTracing:min_trace_length", F(FeatureTracing.MinTraceLengthSec),
            "-FeatureTracing:max_trace_length", F(FeatureTracing.MaxTraceLengthSec),
            "-FeatureTracing:min_isotope_cosine", F(FeatureTracing.MinIsotopeCosineForFeature)
        });

        // Common TOPP
        AddIfFile(args, "-ini", Common.IniFile);
        AddIfFile(args, "-log", Common.LogFile);
        args.AddRange(new[]
        {
            "-instance", Common.Instance.ToString(),
            "-debug", Common.DebugLevel.ToString(),
            "-threads", Common.Threads.ToString()
        });
        AddIfDir(args, "-write_ini", Common.WriteIni);
        AddIfDir(args, "-write_ctd", Common.WriteCtdDir);
        AddIfDir(args, "-write_nested_cwl", Common.WriteNestedCwlDir);
        AddIfDir(args, "-write_cwl", Common.WriteCwlDir);
        AddIfDir(args, "-write_nested_json", Common.WriteNestedJsonDir);
        AddIfDir(args, "-write_json", Common.WriteJsonDir);

        if (Common.NoProgress) args.Add("-no_progress");
        if (Common.Force) args.Add("-force");
        if (Common.TestMode) args.Add("-test");

        return args;
    }

    public string ToCommandLineString() => string.Join(" ", ToArgumentList());

    #endregion

    #region Helpers

    private static void AddIfFile(List<string> args, string switchName, string? path)
    {
        if (!string.IsNullOrWhiteSpace(path))
            args.AddRange(new[] { switchName, Quote(path!) });
    }

    private static void AddIfFiles(List<string> args, string switchName, List<string> paths)
    {
        if (paths.Count > 0)
        {
            args.Add(switchName);
            args.AddRange(paths.Select(Quote));
        }
    }

    private static void AddIfDir(List<string> args, string switchName, string? path)
    {
        if (!string.IsNullOrWhiteSpace(path))
            args.AddRange(new[] { switchName, Quote(path!) });
    }

    private static string Quote(string s)
        => s.Contains(' ') || s.Contains('"') ? "\"" + s.Replace("\"", "\\\"") + "\"" : s;

    private static string F(double d) => d.ToString("G", CultureInfo.InvariantCulture);

    #endregion
}

/// <summary>
/// Convenience extensions for runner integration.
/// </summary>
public static class FlashDeconvOptionsExtensions
{
    public static (int ExitCode, string StdOut, string StdErr) RunWith(
        this FlashDeconvOptions opts,
        FlashDeconvRunner runner,
        int timeoutMs = 0,
        string extraArgs = "")
    {
        var list = opts.ToArgumentList();
        // Split out the required -in and -out so we can pass plain paths to existing runner signature.
        // Assumes standard ordering produced by ToArgumentList()
        string input = GetValueAfter(list, "-in") ?? throw new InvalidOperationException("Missing -in in serialized args.");
        string output = GetValueAfter(list, "-out") ?? throw new InvalidOperationException("Missing -out in serialized args.");

        // Remove -in <file> -out <file> from list for extraArgs concatenation
        var filtered = FilterOut(list, "-in");
        filtered = FilterOut(filtered, "-out");

        var finalExtra = string.Join(" ", filtered.Concat(new[] { extraArgs }).Where(s => !string.IsNullOrWhiteSpace(s)));
        return runner.Run(input.Trim('"'), output.Trim('"'), finalExtra, timeoutMs);
    }

    private static List<string> FilterOut(IReadOnlyList<string> args, string key)
    {
        var result = new List<string>();
        for (int i = 0; i < args.Count; i++)
        {
            if (string.Equals(args[i], key, StringComparison.Ordinal))
            {
                i++; // skip value
                continue;
            }
            result.Add(args[i]);
        }
        return result;
    }

    private static string? GetValueAfter(IReadOnlyList<string> args, string key)
    {
        for (int i = 0; i < args.Count - 1; i++)
            if (args[i] == key)
                return args[i + 1];
        return null;
    }
}