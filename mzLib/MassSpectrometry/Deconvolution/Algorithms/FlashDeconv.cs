using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Text;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    /// <summary>
    /// Thin wrapper around the external FLASHDeconv executable.
    /// 
    /// Design goals (WHY):
    /// - Keep responsibility limited: locating the exe and invoking it.
    /// - Do NOT embed parsing logic until a stable specification for output → IsotopicEnvelope exists.
    /// - Fail fast with informative diagnostics when the executable is missing or returns non‑zero exit.
    /// 
    /// Scope (WHAT):
    /// - Resolve a usable path to FLASHDeconv via (explicit parameter → environment variable → co-located exe name variants).
    /// - Run the process with caller-supplied arguments.
    /// - Capture stdout/stderr for error transparency.
    /// - Surface errors as MzLibException instead of swallowing them.
    /// 
    /// Out of scope (NOT DONE HERE – reasons):
    /// - Output parsing: Deferred because formats (TSV / mzML annotations) may evolve and
    ///   higher-level code can choose which representation to consume.
    /// - Temporary working directory management: This wrapper assumes caller controls paths
    ///   (prevents hidden cleanup bugs and makes test assertions easier).
    /// - Advanced retry / preflight logic: That belongs in a richer runtime (see FlashDeconvRunner elsewhere).
    /// </summary>
    public sealed class FlashDeconv : DeconvolutionAlgorithm
    {
        /// <summary>
        /// Candidate executable file names searched when a path is not explicitly provided.
        /// Rationale: Case differences or packaging variations (build / deployment scripts) may
        /// produce either capitalization. Keeping both avoids brittle assumptions.
        /// </summary>
        private static readonly string[] CandidateExeNames =
        {
            "FLASHDeconv.exe",
            "FlashDeconv.exe"
        };

        /// <summary>
        /// Resolved absolute path to the executable. Immutable after construction to ensure
        /// stable behavior across invocations and thread safety for read‑only access.
        /// </summary>
        private readonly string _exePath;

        /// <summary>
        /// Exposes the resolved executable path for diagnostics / logging.
        /// </summary>
        public string ExecutablePath => _exePath;

        /// <summary>
        /// Constructor resolves the executable path and validates presence.
        /// WHY: Fail early so downstream code does not proceed under false assumptions.
        /// </summary>
        /// <param name="deconParameters">Parameters (currently unused by this wrapper; passed to base to satisfy abstraction).</param>
        /// <param name="exePath">Optional explicit path; if null, falls back to environment and probing logic.</param>
        /// <exception cref="FileNotFoundException">Thrown when no suitable executable can be found.</exception>
        public FlashDeconv(DeconvolutionParameters deconParameters, string? exePath = null)
            : base(deconParameters)
        {
            _exePath = ResolveExecutable(exePath);
            if (!File.Exists(_exePath))
            {
                // WHY include environment hint: guides users to quickest recovery path (set FLASHDECONV_PATH).
                throw new FileNotFoundException(
                    $"FLASHDeconv executable not found at '{_exePath}'. " +
                    "Provide full path or set FLASHDECONV_PATH.", _exePath);
            }
        }

        /// <summary>
        /// Resolve candidate executable path.
        /// Resolution strategy (WHY this order):
        /// 1. Explicit path (caller intent should always win).
        /// 2. Environment variable FLASHDECONV_PATH (enables deployment flexibility / CI overrides).
        /// 3. Co-located executable names (common for bundling with the library).
        /// 4. Fallback to the first candidate name in assembly directory (gives deterministic message on failure).
        /// 
        /// No validation of existence here (deferred to caller) so path computation remains pure.
        /// </summary>
        private static string ResolveExecutable(string? explicitPath)
        {
            // 1. Caller-supplied explicit path.
            if (!string.IsNullOrWhiteSpace(explicitPath))
                return Path.GetFullPath(explicitPath);

            // 2. Environment override.
            var env = Environment.GetEnvironmentVariable("FLASHDECONV_PATH");
            if (!string.IsNullOrWhiteSpace(env))
                return Path.GetFullPath(env);

            // 3. Probe alongside the executing assembly (robust for xcopy deployments / self-contained distributions).
            string assemblyDir = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
                                 ?? AppDomain.CurrentDomain.BaseDirectory;

            foreach (var name in CandidateExeNames)
            {
                var candidate = Path.Combine(assemblyDir, name);
                if (File.Exists(candidate))
                    return candidate;
            }

            // 4. Deterministic fallback (will fail later with a clear exception if not present).
            return Path.Combine(assemblyDir, CandidateExeNames[0]);
        }

        /// <summary>
        /// Future hook: Convert FLASHDeconv output to IsotopicEnvelope collection.
        /// WHY not implemented:
        /// - Avoid premature coupling to output format before stable internal consumption spec.
        /// - Encourages explicit parsing modules (e.g., a FlashDeconvResultParser) for SRP compliance.
        /// </summary>
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            throw new NotImplementedException("FLASHDeconv result parsing to IsotopicEnvelope not implemented.");
        }

        /// <summary>
        /// Runs the external FLASHDeconv process with basic argument composition.
        /// WHAT:
        /// - Ensures input exists, output directory is created.
        /// - Streams stdout/stderr asynchronously to avoid deadlocks (WHY: large buffers can block if not drained).
        /// - Throws on non-zero exit with full diagnostic context (command + captured IO).
        /// 
        /// Caller responsibility (WHY left external):
        /// - Selecting additional arguments (tuning / filtering / output format).
        /// - Post-run validation or parsing.
        /// - Idempotent cleanup (temp files or large intermediates).
        /// </summary>
        /// <param name="inputFile">Input mzML path consumed by FLASHDeconv.</param>
        /// <param name="outputFile">Primary TSV (or feature) output path expected from FLASHDeconv.</param>
        /// <param name="extraArgs">Raw argument tail appended verbatim (caller escapes quoting as needed).</param>
        /// <exception cref="ArgumentException">Thrown if input/output paths are null/empty.</exception>
        /// <exception cref="FileNotFoundException">Thrown if input file cannot be located.</exception>
        /// <exception cref="MzLibException">Thrown if process returns non‑zero exit code.</exception>
        public void Run(string inputFile, string outputFile, string? extraArgs = null)
        {
            // Validate inputs early for predictable failure.
            if (string.IsNullOrWhiteSpace(inputFile))
                throw new ArgumentException("Input file path is null/empty.", nameof(inputFile));
            if (!File.Exists(inputFile))
                throw new FileNotFoundException("Input file not found.", inputFile);
            if (string.IsNullOrWhiteSpace(outputFile))
                throw new ArgumentException("Output file path is null/empty.", nameof(outputFile));

            // Ensure output directory exists (WHY: FLASHDeconv does not create deep directories itself).
            Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(outputFile))!);
            extraArgs ??= string.Empty;

            // Build ProcessStartInfo:
            // - UseShellExecute = false + redirection → enables controlled capture.
            // - WorkingDirectory set to executable folder (WHY: some tools load relative assets / INI defaults).
            var psi = new ProcessStartInfo
            {
                FileName = _exePath,
                Arguments = $"-in \"{inputFile}\" -out \"{outputFile}\" {extraArgs}".Trim(),
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                WorkingDirectory = Path.GetDirectoryName(_exePath) ?? Environment.CurrentDirectory
            };

            using var proc = new Process { StartInfo = psi };
            var stdout = new StringBuilder();
            var stderr = new StringBuilder();

            // Asynchronous handlers prevent potential deadlock if one stream fills before WaitForExit completes.
            proc.OutputDataReceived += (_, e) => { if (e.Data != null) stdout.AppendLine(e.Data); };
            proc.ErrorDataReceived  += (_, e) => { if (e.Data != null) stderr.AppendLine(e.Data); };

            if (!proc.Start())
                throw new InvalidOperationException("Failed to start FLASHDeconv.");

            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();

            // Non‑zero exit → surface full context (WHY: aids reproducibility in CI and remote debugging).
            if (proc.ExitCode != 0)
            {
                throw new MzLibException(
                    $"FLASHDeconv exited with code {proc.ExitCode}.{Environment.NewLine}" +
                    $"Exe: {_exePath}{Environment.NewLine}" +
                    $"Args: {psi.Arguments}{Environment.NewLine}" +
                    $"StdOut:{Environment.NewLine}{stdout}{Environment.NewLine}" +
                    $"StdErr:{Environment.NewLine}{stderr}");
            }
        }
    }
}
