using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace MassSpectrometry.FlashDeconvRuntime
{
    /// <summary>
    /// Strongly-typed representation of the (verbose) FLASHDeconv command line surface.
    /// 
    /// WHAT:
    /// - Encapsulates every documented switch from the --helphelp output (including advanced / per?MS level options).
    /// - Provides serialization (ToArgumentList / ToCommandLineString) that converts the object graph into
    ///   a flat sequence of CLI tokens suitable for spawning the process.
    /// 
    /// WHY:
    /// - Prevent fragile string concatenations scattered across calling code.
    /// - Allow compile?time discoverability and IntelliSense for available options.
    /// - Centralize quoting, numeric formatting, and conditional argument emission.
    /// - Enforce required arguments (-in / -out) via Validate().
    /// 
    /// Non-goals:
    /// - This class deliberately does NOT execute the process (separation of concerns).
    /// - It does not perform semantic validation beyond presence of mandatory paths (tool retains authority).
    /// - No attempt is made to infer defaults or skip emitting values equal to tool defaults (kept explicit
    ///   for reproducibility and easier debugging).
    /// </summary>
    public sealed class FlashDeconvOptions
    {
        #region Core Required Paths

        /// <summary>
        /// Input mzML file (maps to -in). Required by Validate().
        /// Keeping it nullable allows staged construction (e.g., DI container sets properties post-creation).
        /// </summary>
        public string? InputMzMlPath { get; set; }            // -in

        /// <summary>
        /// Primary feature TSV output (maps to -out). Required by Validate().
        /// </summary>
        public string? OutputTsvPath { get; set; }            // -out

        #endregion

        #region Optional Direct Output Variants

        /// <summary>
        /// Per-MS-level deconvolved spectra TSVs (-out_spec). Multiple file paths allowed.
        /// WHY List: tool expects one file per MS level in order.
        /// </summary>
        public List<string> OutSpectraPerMsLevelTsv { get; } = new();   // -out_spec

        /// <summary>
        /// Combined deconvolved spectra written as mzML (-out_mzml).
        /// </summary>
        public string? OutMzmlDeconvolved { get; set; }                 // -out_mzml

        /// <summary>
        /// Annotated mzML with isotope/charge annotations added as meta data (-out_annotated_mzml).
        /// </summary>
        public string? OutAnnotatedMzml { get; set; }                   // -out_annotated_mzml

        /// <summary>
        /// Promex-compatible ms1ft export (-out_promex).
        /// </summary>
        public string? OutPromex { get; set; }                          // -out_promex

        /// <summary>
        /// TopFD-compatible msalign output files per MS level (-out_topFD).
        /// </summary>
        public List<string> OutTopFDPerMsLevelMsalign { get; } = new(); // -out_topFD

        /// <summary>
        /// TopFD-compatible feature files per MS level (-out_topFD_feature).
        /// </summary>
        public List<string> OutTopFDFeaturePerMsLevel { get; } = new(); // -out_topFD_feature

        #endregion

        #region FLASHIda Coupling

        /// <summary>
        /// FLASHIda acquisition log (-in_log) used to integrate dynamic acquisition metadata.
        /// </summary>
        public string? FlashIdaLog { get; set; }                        // -in_log

        #endregion

        #region Target Precursor Overrides

        /// <summary>
        /// Minimum precursor SNR for identification (topFD msalign outputs only) (-min_precursor_snr).
        /// WHY double: preserves precision for thresholds like 1.25 etc.
        /// </summary>
        public double TargetPrecursorSnrMin { get; set; } = 1.0;        // -min_precursor_snr

        /// <summary>
        /// Fixed precursor charge for targeted / specialized workflows (-target_precursor_charge).
        /// 0 => disabled (tool infers).
        /// </summary>
        public int TargetPrecursorCharge { get; set; } = 0;             // -target_precursor_charge

        /// <summary>
        /// Fixed precursor m/z to pair with TargetPrecursorCharge (-target_precursor_mz).
        /// </summary>
        public double TargetPrecursorMz { get; set; } = 0.0;            // -target_precursor_mz

        #endregion

        #region Output Formatting / Scope

        /// <summary>
        /// Controls charge assignment to deconvolved masses in mzML output (-mzml_mass_charge) allowed values -1,0,+1.
        /// 0 => uncharged masses.
        /// </summary>
        public int MzmlMassChargeMode { get; set; } = 0;                // -mzml_mass_charge

        /// <summary>
        /// Number of preceding MS1 spectra to consider for precursor determination (-preceding_MS1_count).
        /// WHY: In top-down data a precursor may appear in earlier MS1 than immediately preceding.
        /// </summary>
        public int PrecedingMs1Count { get; set; } = 3;                 // -preceding_MS1_count

        /// <summary>
        /// Emits detailed per-peak info in deconvolved spectra TSV (-write_detail).
        /// Tradeoff: richer diagnostics vs. larger file size.
        /// </summary>
        public bool WriteDetail { get; set; } = false;                  // -write_detail

        /// <summary>
        /// Maximum MS level analyzed (-max_MS_level). Constrains runtime and noise.
        /// </summary>
        public int MaxMsLevel { get; set; } = 3;                        // -max_MS_level

        /// <summary>
        /// Force all spectra to a single MS level (-forced_MS_level) for unusual data sets (e.g. only MS2).
        /// 0 => disabled.
        /// </summary>
        public int ForcedMsLevel { get; set; } = 0;                     // -forced_MS_level

        /// <summary>
        /// Pre-deconvolution spectral merging strategy (-merging_method).
        /// </summary>
        public MergingMethodEnum MergingMethod { get; set; } = MergingMethodEnum.None; // -merging_method

        /// <summary>
        /// Whether to compute/report FDR-like q-values (-report_FDR) (beta feature).
        /// </summary>
        public bool ReportFdr { get; set; } = false;                    // -report_FDR

        /// <summary>
        /// Use RNA averagine model (-use_RNA_averagine) for transcript-related analyses.
        /// </summary>
        public bool UseRnaAveragine { get; set; } = false;              // -use_RNA_averagine

        #endregion

        #region Nested Option Groups (mirroring CLI hierarchy)

        /// <summary>
        /// Algorithm-centric thresholds and search bounds (maps to -Algorithm:*).
        /// </summary>
        public AlgorithmGroup Algorithm { get; } = new();

        /// <summary>
        /// Feature tracing / quantification tuning (maps to -FeatureTracing:*).
        /// </summary>
        public FeatureTracingGroup FeatureTracing { get; } = new();

        /// <summary>
        /// Shared TOPP (OpenMS) generic options (maps to misc root flags).
        /// </summary>
        public CommonTopGroup Common { get; } = new();

        #endregion

        #region Nested Types

        /// <summary>
        /// Enum for merging method selection (explicit ints ensure stable CLI emission).
        /// </summary>
        public enum MergingMethodEnum { None = 0, GaussianAveraging = 1, Block = 2 }

        /// <summary>
        /// Encapsulates core algorithmic constraints/scoring thresholds.
        /// WHY separate class: logically grouped; simplifies property discovery and potential reuse.
        /// </summary>
        public sealed class AlgorithmGroup
        {
            /// <summary>
            /// Per-MS-level ppm tolerances (-Algorithm:tol). A list allows specifying MS1, MS2, ... individually.
            /// Tool expects positional mapping (order matters).
            /// </summary>
            public List<double> PpmTolerancePerMsLevel { get; } = new() { 10.0, 10.0, 10.0 }; // -Algorithm:tol

            public double MinMass { get; set; } = 50.0;                // -Algorithm:min_mass
            public double MaxMass { get; set; } = 1.0e5;               // -Algorithm:max_mass
            public int MinCharge { get; set; } = 1;                    // -Algorithm:min_charge
            public int MaxCharge { get; set; } = 100;                  // -Algorithm:max_charge

            /// <summary>
            /// Minimum m/z bound; <= 0 disables (tool internal default).
            /// </summary>
            public double MinMz { get; set; } = -1.0;                  // -Algorithm:min_mz

            /// <summary>
            /// Maximum m/z bound; <= 0 disables (tool internal default).
            /// </summary>
            public double MaxMz { get; set; } = -1.0;                  // -Algorithm:max_mz

            public double MinRt { get; set; } = -1.0;                  // -Algorithm:min_rt
            public double MaxRt { get; set; } = -1.0;                  // -Algorithm:max_rt
            public double IsolationWindowWidth { get; set; } = 5.0;    // -Algorithm:isolation_window

            /// <summary>
            /// Per-MS-level isotope pattern cosine similarity thresholds (list structure mirrors tolerance list).
            /// </summary>
            public List<double> MinIsotopeCosinePerMsLevel { get; } = new() { 0.85, 0.85, 0.85 }; // -Algorithm:min_isotope_cosine

            public int AllowedIsotopeError { get; set; } = 1;          // -Algorithm:allowed_isotope_error
            public double MinIntensity { get; set; } = 0.0;            // -Algorithm:min_intensity
        }

        /// <summary>
        /// Feature tracing refinement and quantification mode selection group.
        /// </summary>
        public sealed class FeatureTracingGroup
        {
            public double MassErrorPpm { get; set; } = -1.0;                 // -FeatureTracing:mass_error_ppm
            public double IonMobilityTolerance { get; set; } = 0.01;          // -FeatureTracing:ion_mobility_tolerance
            public QuantMethodEnum QuantMethod { get; set; } = QuantMethodEnum.area; // -FeatureTracing:quant_method
            public double MinSampleRate { get; set; } = 0.05;                 // -FeatureTracing:min_sample_rate
            public double MinTraceLengthSec { get; set; } = 10.0;             // -FeatureTracing:min_trace_length
            public double MaxTraceLengthSec { get; set; } = -1.0;             // -FeatureTracing:max_trace_length
            public double MinIsotopeCosineForFeature { get; set; } = -1.0;    // -FeatureTracing:min_isotope_cosine

            /// <summary>
            /// Quantification choices. Kept lower-case to match CLI tokens verbatim.
            /// </summary>
            public enum QuantMethodEnum { area, median, max_height }
        }

        /// <summary>
        /// Generic OpenMS / TOPP options that are not algorithm-specific.
        /// Some of these (e.g., -write_ini) alter control flow or produce side outputs.
        /// </summary>
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
        /// Ensures mandatory parameters are populated prior to serialization.
        /// WHY: Fail here (deterministic) instead of letting the external tool error unpredictably.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown if required path is null/empty.</exception>
        public void Validate()
        {
            if (string.IsNullOrWhiteSpace(InputMzMlPath))
                throw new ArgumentException("InputMzMlPath (-in) is required.");
            if (string.IsNullOrWhiteSpace(OutputTsvPath))
                throw new ArgumentException("OutputTsvPath (-out) is required.");
        }

        /// <summary>
        /// Produces an ordered list of CLI tokens (switch/value tokens kept distinct).
        /// WHY list-of-strings: safer for ProcessStartInfo usage (no ambiguous string concatenation).
        /// 
        /// Implementation notes:
        /// - Always emits -in / -out first (predictable ordering simplifies consumption).
        /// - Only emits optional outputs if non-null or list non-empty.
        /// - Multi-valued options (e.g., tolerances) expand inline following the switch.
        /// - Values are quoted only when necessary (embedded space or quote).
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

            // Outputs (lists first)
            AddIfFiles(args, "-out_spec", OutSpectraPerMsLevelTsv);
            AddIfFiles(args, "-out_topFD", OutTopFDPerMsLevelMsalign);
            AddIfFiles(args, "-out_topFD_feature", OutTopFDFeaturePerMsLevel);

            AddIfFile(args, "-out_mzml", OutMzmlDeconvolved);
            AddIfFile(args, "-out_annotated_mzml", OutAnnotatedMzml);
            AddIfFile(args, "-out_promex", OutPromex);

            // Scalar simple group
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

            // FeatureTracing group
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

            // Common TOPP group (conditionally emitted)
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

        /// <summary>
        /// Convenience method that flattens arguments into a single string.
        /// Useful for logging or diagnostics. Not recommended for direct process spawning (list safer).
        /// </summary>
        public string ToCommandLineString() => string.Join(" ", ToArgumentList());

        #endregion

        #region Helpers (Internal Argument Assembly)

        /// <summary>
        /// Append switch + single file path when present.
        /// WHY: Avoids if clutter in main serialization method.
        /// </summary>
        private static void AddIfFile(List<string> args, string switchName, string? path)
        {
            if (!string.IsNullOrWhiteSpace(path))
                args.AddRange(new[] { switchName, Quote(path!) });
        }

        /// <summary>
        /// Append switch followed by multiple file paths.
        /// Only emits if the list is non-empty.
        /// </summary>
        private static void AddIfFiles(List<string> args, string switchName, List<string> paths)
        {
            if (paths.Count > 0)
            {
                args.Add(switchName);
                args.AddRange(paths.Select(Quote));
            }
        }

        /// <summary>
        /// Semantically identical to AddIfFile; separated for clarity (directory semantics).
        /// </summary>
        private static void AddIfDir(List<string> args, string switchName, string? path)
        {
            if (!string.IsNullOrWhiteSpace(path))
                args.AddRange(new[] { switchName, Quote(path!) });
        }

        /// <summary>
        /// Adds quotes only when necessary (spaces or embedded quotes).
        /// WHY: Minimizes noise while preserving safety.
        /// </summary>
        private static string Quote(string s)
            => s.Contains(' ') || s.Contains('"') ? "\"" + s.Replace("\"", "\\\"") + "\"" : s;

        /// <summary>
        /// Formats doubles with invariant culture and "G" to keep CLI-friendly compactness and precision.
        /// </summary>
        private static string F(double d) => d.ToString("G", CultureInfo.InvariantCulture);

        #endregion
    }

    /// <summary>
    /// Extension helpers bridging FlashDeconvOptions ? FlashDeconvRunner.
    /// WHY separate static class: keeps core options model free of execution semantics.
    /// </summary>
    public static class FlashDeconvOptionsExtensions
    {
        /// <summary>
        /// Executes FLASHDeconv using a pre-configured runner.
        /// 
        /// Flow:
        /// 1. Serialize options to argument list.
        /// 2. Extract -in/-out (runner API expects them separately for clarity).
        /// 3. Remove those tokens from the list.
        /// 4. Concatenate remainder into extraArgs and invoke runner.
        /// 
        /// WHY not re-parse after join: we already have discrete tokens; joining preserves original order and spacing.
        /// </summary>
        public static (int ExitCode, string StdOut, string StdErr) RunWith(
            this FlashDeconvOptions opts,
            FlashDeconvRunner runner,
            int timeoutMs = 0,
            string extraArgs = "")
        {
            var list = opts.ToArgumentList();

            // Locate required tokens; if missing, serialization invariant broken (fail fast).
            string input = GetValueAfter(list, "-in") ?? throw new InvalidOperationException("Missing -in in serialized args.");
            string output = GetValueAfter(list, "-out") ?? throw new InvalidOperationException("Missing -out in serialized args.");

            // Strip -in / -out pairs to avoid duplication when passing back to runner.
            var filtered = FilterOut(list, "-in");
            filtered = FilterOut(filtered, "-out");

            // Append caller-supplied extra args (if any) to the remainder.
            var finalExtra = string.Join(" ", filtered.Concat(new[] { extraArgs }).Where(s => !string.IsNullOrWhiteSpace(s)));
            return runner.Run(input.Trim('"'), output.Trim('"'), finalExtra, timeoutMs);
        }

        /// <summary>
        /// Removes a switch and its immediate value token (assumes well-formed list).
        /// </summary>
        private static List<string> FilterOut(IReadOnlyList<string> args, string key)
        {
            var result = new List<string>();
            for (int i = 0; i < args.Count; i++)
            {
                if (string.Equals(args[i], key, StringComparison.Ordinal))
                {
                    i++; // skip associated value
                    continue;
                }
                result.Add(args[i]);
            }
            return result;
        }

        /// <summary>
        /// Returns the token immediately following the specified key switch.
        /// Returns null if not found or key is last element (defensive).
        /// </summary>
        private static string? GetValueAfter(IReadOnlyList<string> args, string key)
        {
            for (int i = 0; i < args.Count - 1; i++)
                if (args[i] == key)
                    return args[i + 1];
            return null;
        }
    }
}