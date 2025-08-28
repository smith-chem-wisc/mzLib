using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Diagnostics;
using System.Threading.Tasks;

namespace MassSpectrometry.FlashDeconvRuntime
{
    /// <summary>
    /// Higher-level orchestration runner for executing FLASHDeconv.
    /// 
    /// WHAT:
    /// - Locates executable (FlashDeconvLocator).
    /// - Optional preflight (help call) to detect missing dependencies early.
    /// - Provides sync & async execution paths with optional timeout & cancellation.
    /// - Adds OPENMS_DATA_PATH inference for environments where it's not configured.
    /// 
    /// WHY this layer (vs. direct Process use per call site):
    /// - Consolidates argument shaping & diagnostics logic.
    /// - Reduces duplication and chance of inconsistent error handling.
    /// - Exposes environment-variable overrides for non-code configuration (CI, containers).
    /// - Leaves parsing and domain interpretation to higher-level code (single responsibility).
    /// </summary>
    public sealed class FlashDeconvRunner
    {
        private readonly string _exe;

        /// <summary>
        /// Construct a runner, optionally performing a preflight.
        /// Preflight is enabled by default to fail fast on missing DLLs / share data, but can be skipped for speed.
        /// </summary>
        public FlashDeconvRunner(
            string? overridePath = null,
            bool preflight = true,
            int? preflightTimeoutMs = 30000,
            string? preflightArg = null)
        {
            // Resolve path first; throw early if absent for easier debugging.
            _exe = FlashDeconvLocator.FindExecutable(overridePath)
                   ?? throw new FileNotFoundException("FLASHDeconv executable not found (env FLASHDECONV_PATH or vendored runtime).");

            // Environment can force skip (typical for iterative dev loops or when dependency warmup is known slow).
            var skip = Environment.GetEnvironmentVariable("FLASHDECONV_SKIP_PREFLIGHT");
            if (string.Equals(skip, "1", StringComparison.OrdinalIgnoreCase) ||
                string.Equals(skip, "true", StringComparison.OrdinalIgnoreCase))
                preflight = false;

            if (!preflight)
                return;

            // Allow external runtime to tune preflight timeout without code change (e.g., GH Actions vs local).
            if (Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_TIMEOUT_MS") is string tov &&
                int.TryParse(tov, out var envTimeout))
                preflightTimeoutMs = envTimeout;

            // Use environment override or default to "-h" (some builds may prefer "--help" or blank).
            preflightArg ??= Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_ARG") ?? "-h";

            if (!FlashDeconvLocator.Preflight(
                    _exe,
                    out var so,
                    out var se,
                    out var ec,
                    preflightTimeoutMs,
                    preflightArg))
            {
                // Provide targeted hints to accelerate user remediation.
                var msg = new StringBuilder()
                    .AppendLine("FLASHDeconv preflight failed.")
                    .AppendLine($"  Path   : {_exe}")
                    .AppendLine($"  Arg    : {preflightArg}")
                    .AppendLine($"  Exit   : {ec}")
                    .AppendLine($"  Timeout: {(preflightTimeoutMs is null or <= 0 ? "None (infinite)" : preflightTimeoutMs + " ms")}")
                    .AppendLine($"  StdOut : {(string.IsNullOrWhiteSpace(so) ? "<empty>" : Trim(so, 400))}")
                    .AppendLine($"  StdErr : {(string.IsNullOrWhiteSpace(se) ? "<empty>" : Trim(se, 400))}")
                    .AppendLine("Hints:")
                    .AppendLine("  - Skip preflight: set FLASHDECONV_SKIP_PREFLIGHT=1 or new FlashDeconvRunner(preflight:false)")
                    .AppendLine("  - Increase timeout: set FLASHDECONV_PREFLIGHT_TIMEOUT_MS=60000")
                    .AppendLine("  - Adjust help arg: set FLASHDECONV_PREFLIGHT_ARG=--help (or blank)")
                    .AppendLine("  - If help hangs, skip preflight and run directly.")
                    .ToString();
                throw new InvalidOperationException(msg);
            }
        }

        /// <summary>
        /// Convenience synchronous API that delegates to RunAsync (single implementation surface).
        /// </summary>
        public (int ExitCode, string StdOut, string StdErr) Run(
            string inputMzml, string outputTsv, string extraArgs = "", int timeoutMs = 0)
            => RunAsync(inputMzml, outputTsv, extraArgs, timeoutMs, CancellationToken.None).GetAwaiter().GetResult();

        /// <summary>
        /// Asynchronously execute FLASHDeconv.
        /// - Always supplies -in and -out.
        /// - Appends any additional caller-specified arguments verbatim (caller handles quoting).
        /// - Captures stdout/stderr (line-buffered) to StringBuilder to avoid blocking on large output.
        /// - Supports either a wall-clock timeout (timeoutMs > 0) or external cancellation token.
        /// - Returns exit code instead of throwing (unlike simpler wrapper) to allow nuanced caller handling.
        /// </summary>
        public async Task<(int ExitCode, string StdOut, string StdErr)> RunAsync(
            string inputMzml,
            string outputTsv,
            string extraArgs = "",
            int timeoutMs = 0,
            CancellationToken ct = default)
        {
            if (!File.Exists(inputMzml))
                throw new FileNotFoundException("Input file missing", inputMzml);

            // Guarantee directory exists (tool does not create multi-level paths).
            Directory.CreateDirectory(Path.GetDirectoryName(outputTsv)!);

            var args = new StringBuilder()
                .Append("-in ").Append('"').Append(inputMzml).Append('"')
                .Append(" -out ").Append('"').Append(outputTsv).Append('"');

            if (!string.IsNullOrWhiteSpace(extraArgs))
                args.Append(' ').Append(extraArgs);

            var psi = new ProcessStartInfo
            {
                FileName = _exe,
                Arguments = args.ToString(),
                UseShellExecute = false,             // Required for redirection & security context.
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true
            };

            // Heuristic addition: autoconfigure OPENMS_DATA_PATH if absent to reduce friction in CI
            // (mirrors logic inside the simpler FlashDeconv wrapper so both paths behave uniformly).
            if (string.IsNullOrWhiteSpace(psi.Environment["OPENMS_DATA_PATH"]))
            {
                TryInferOpenMsDataPath(Path.GetDirectoryName(_exe)!, p => psi.Environment["OPENMS_DATA_PATH"] = p);
            }

            using var p = new Process { StartInfo = psi };
            var so = new StringBuilder();
            var se = new StringBuilder();

            p.OutputDataReceived += (_, e) => { if (e.Data != null) so.AppendLine(e.Data); };
            p.ErrorDataReceived  += (_, e) => { if (e.Data != null) se.AppendLine(e.Data); };

            if (!p.Start())
                throw new InvalidOperationException("Failed to start FLASHDeconv process.");

            p.BeginOutputReadLine();
            p.BeginErrorReadLine();

            // Combine caller cancellation with internal timeout.
            using var linkedCts = CancellationTokenSource.CreateLinkedTokenSource(ct);
            if (timeoutMs > 0)
                linkedCts.CancelAfter(timeoutMs);

            try
            {
                // WaitForExit is blocking; wrap in Task.Run so cancellation can interrupt.
                await Task.Run(() => p.WaitForExit(), linkedCts.Token);
            }
            catch (OperationCanceledException)
            {
                // If we canceled due to timeout (process still running), kill and surface a TimeoutException.
                if (!p.HasExited)
                {
                    try { p.Kill(true); } catch { /* ignore race conditions */ }
                    if (timeoutMs > 0 && !ct.IsCancellationRequested)
                        throw new TimeoutException($"FLASHDeconv run timed out after {timeoutMs} ms");
                    throw; // propagate external cancellation
                }
            }

            return (p.ExitCode, so.ToString(), se.ToString());
        }

        /// <summary>
        /// Attempt to infer OpenMS shared data path relative to the executable (common layouts).
        /// WHY: Many CI / local setups bundle share/OpenMS but do not export OPENMS_DATA_PATH.
        /// If inference fails, we let the tool emit its own fatal error for clarity.
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
                // Silent: best-effort only.
            }
        }

        /// <summary>
        /// Utility for truncating large captured output to keep logs concise.
        /// </summary>
        private static string Trim(string s, int max) =>
            s.Length <= max ? s : s.Substring(0, max) + "...";
    }
}