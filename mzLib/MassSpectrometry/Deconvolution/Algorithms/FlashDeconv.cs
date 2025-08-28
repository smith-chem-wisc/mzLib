#nullable enable
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
    /// Thin, responsibility‑focused wrapper around the external FLASHDeconv executable.
    /// 
    /// WHAT this class does:
    /// - Resolves a usable FLASHDeconv executable path (explicit → env var → co‑located names).
    /// - Executes the process with caller‑supplied arguments (-in / -out plus any extra).
    /// - Captures stdout / stderr asynchronously and throws rich diagnostics on non‑zero exit.
    /// 
    /// WHY this exists (and remains intentionally minimal):
    /// - Keeps deconvolution framework decoupled from external binary lifecycle.
    /// - Avoids prematurely locking into an output parsing contract (tsv / mzML may evolve).
    /// - Allows higher-level orchestration (FlashDeconvRunner) to add preflight, timeouts, etc.
    /// 
    /// NOT handled here:
    /// - Parsing of produced feature/spectra files → deferred to a future, spec‑driven parser.
    /// - Preflight dependency probing (done by FlashDeconvRunner / Locator).
    /// - Temporary directory management (caller retains explicit control for testability).
    /// </summary>
    public sealed class FlashDeconv : DeconvolutionAlgorithm
    {
        /// <summary>
        /// Candidate executable file names (WHY two: packaging / OS build scripts sometimes differ in casing).
        /// </summary>
        private static readonly string[] CandidateExeNames =
        {
            "FLASHDeconv.exe",
            "FlashDeconv.exe"
        };

        /// <summary>
        /// Resolved absolute path after constructor logic. Immutable for thread safety & reproducibility.
        /// </summary>
        private readonly string _exePath;

        /// <summary>
        /// Public exposure for logging / diagnostics (e.g., tests printing exact binary used).
        /// </summary>
        public string ExecutablePath => _exePath;

        /// <summary>
        /// Resolve the executable path immediately and fail fast if missing.
        /// WHY fail here: prevents deferred runtime surprises deep inside pipelines.
        /// </summary>
        public FlashDeconv(DeconvolutionParameters deconParameters, string? exePath = null)
            : base(deconParameters)
        {
            _exePath = ResolveExecutable(exePath);
            if (!File.Exists(_exePath))
            {
                // Include remediation hint to minimize user friction.
                throw new FileNotFoundException(
                    $"FLASHDeconv executable not found at '{_exePath}'. Provide full path or set FLASHDECONV_PATH.", _exePath);
            }
        }

        /// <summary>
        /// Path resolution strategy (order reflects increasing generality):
        /// 1. Explicit user-supplied path (highest intent).
        /// 2. Environment variable FLASHDECONV_PATH (allows CI / user profile overrides).
        /// 3. Probe for bundled filenames in the executing assembly directory.
        /// 4. Fallback (first candidate) → will trigger missing file exception if absent.
        /// 
        /// WHY no existence check except for #3: we keep function pure; caller validates final.
        /// </summary>
        private static string ResolveExecutable(string? explicitPath)
        {
            if (!string.IsNullOrWhiteSpace(explicitPath))
                return Path.GetFullPath(explicitPath);

            var env = Environment.GetEnvironmentVariable("FLASHDECONV_PATH");
            if (!string.IsNullOrWhiteSpace(env))
                return Path.GetFullPath(env);

            string assemblyDir = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
                                 ?? AppDomain.CurrentDomain.BaseDirectory;

            foreach (var name in CandidateExeNames)
            {
                var candidate = Path.Combine(assemblyDir, name);
                if (File.Exists(candidate))
                    return candidate;
            }

            return Path.Combine(assemblyDir, CandidateExeNames[0]);
        }

        /// <summary>
        /// Stub required by abstract base. Intentionally throws to highlight that parsing
        /// logic has not yet been integrated (WHY: avoid premature commitment).
        /// </summary>
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range) =>
            throw new NotImplementedException("FLASHDeconv result parsing to IsotopicEnvelope not implemented.");

        /// <summary>
        /// Execute FLASHDeconv with minimal ceremony.
        /// 
        /// WHAT happens:
        /// - Validates input/output parameters (fast deterministic failures).
        /// - Ensures output directory exists (tool does not create nested directories).
        /// - Adds heuristic OPENMS_DATA_PATH detection if user did not supply it (why: reduce CI friction).
        /// - Starts process with asynchronous read of stdout/stderr (prevents deadlock when buffers fill).
        /// - On non‑zero exit code, throws MzLibException containing full captured streams & command line.
        /// 
        /// Caller retains responsibility for:
        /// - Choosing extra arguments (filtering, output expansions, etc.).
        /// - Post-run validation / parsing / cleanup.
        /// - Deciding whether to suppress or propagate exceptions.
        /// </summary>
        public void Run(string inputFile, string outputFile, string? extraArgs = null)
        {
            // Validate early for clearer error messages than a tool crash later.
            if (string.IsNullOrWhiteSpace(inputFile))
                throw new ArgumentException("Input file path is null/empty.", nameof(inputFile));
            if (!File.Exists(inputFile))
                throw new FileNotFoundException("Input file not found.", inputFile);
            if (string.IsNullOrWhiteSpace(outputFile))
                throw new ArgumentException("Output file path is null/empty.", nameof(outputFile));

            Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(outputFile))!);
            extraArgs ??= string.Empty;

            var psi = new ProcessStartInfo
            {
                FileName = _exePath,
                Arguments = $"-in \"{inputFile}\" -out \"{outputFile}\" {extraArgs}".Trim(),
                UseShellExecute = false,             // Needed for stream redirection (& security).
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                WorkingDirectory = Path.GetDirectoryName(_exePath) ?? Environment.CurrentDirectory
            };

            // Heuristic: If OPENMS_DATA_PATH not set, attempt to infer typical share/OpenMS layouts relative to exe.
            if (string.IsNullOrWhiteSpace(psi.Environment["OPENMS_DATA_PATH"]))
            {
                TryInferOpenMsDataPath(Path.GetDirectoryName(_exePath)!, p => psi.Environment["OPENMS_DATA_PATH"] = p);
            }

            using var proc = new Process { StartInfo = psi };
            var stdout = new StringBuilder();
            var stderr = new StringBuilder();

            proc.OutputDataReceived += (_, e) => { if (e.Data != null) stdout.AppendLine(e.Data); };
            proc.ErrorDataReceived += (_, e) => { if (e.Data != null) stderr.AppendLine(e.Data); };

            if (!proc.Start())
                throw new InvalidOperationException("Failed to start FLASHDeconv.");

            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();

            if (proc.ExitCode != 0)
            {
                // Provide exhaustive diagnostics for reproducibility (CI logs, user bug reports).
                throw new MzLibException(
                    $"FLASHDeconv exited with code {proc.ExitCode}.{Environment.NewLine}" +
                    $"Exe: {_exePath}{Environment.NewLine}" +
                    $"Args: {psi.Arguments}{Environment.NewLine}" +
                    $"StdOut:{Environment.NewLine}{stdout}{Environment.NewLine}" +
                    $"StdErr:{Environment.NewLine}{stderr}");
            }
        }

        /// <summary>
        /// Attempt to discover the OpenMS shared data directory relative to the executable.
        /// WHY heuristic: CI or local dev often bundles share/OpenMS adjacent to the binary without env config.
        /// </summary>
        private static void TryInferOpenMsDataPath(string exeDir, Action<string> setter)
        {
            try
            {
                var candidates = new[]
                {
                    Path.Combine(exeDir, "share", "OpenMS"),
                    Path.GetFullPath(Path.Combine(exeDir, "..", "share", "OpenMS")),
                    Path.GetFullPath(Path.Combine(exeDir, "..", "..", "share", "OpenMS"))
                };
                foreach (var c in candidates)
                {
                    if (Directory.Exists(c))
                    {
                        setter(c);
                        break;
                    }
                }
            }
            catch
            {
                // Swallow silently: absence simply allows the tool to throw its own explicit fatal error.
            }
        }
    }
}
