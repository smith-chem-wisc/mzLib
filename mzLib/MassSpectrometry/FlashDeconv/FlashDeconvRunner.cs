using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Diagnostics;
using System.Threading.Tasks;

namespace MassSpectrometry.FlashDeconvRuntime
{
    /// <summary>
    /// High-level convenience runner for the FLASHDeconv executable.
    /// 
    /// WHAT:
    /// - Locates the executable (delegates to FlashDeconvLocator).
    /// - Optionally performs a "preflight" help invocation to validate environment / dependencies.
    /// - Executes FLASHDeconv with provided arguments (sync or async).
    /// - Supports optional timeout enforcement and cancellation.
    /// 
    /// WHY this abstraction (vs. directly spawning Process everywhere):
    /// - Centralizes preflight logic, error reporting, and timeout handling.
    /// - Reduces duplication in tests and integration code.
    /// - Provides strongly documented environment variable overrides to alter behavior without code changes.
    /// 
    /// Design decisions:
    /// - Preflight is opt-out (default true) because early failure is usually preferable to a late crash.
    /// - Timeout strategy: separate configurable preflight timeout vs. run timeout (both can be infinite).
    /// - Streams (stdout/stderr) are buffered in memory because typical diagnostic output is short; for huge output
    ///   volumes, a future enhancement could stream to disk or expose IProgress handlers.
    /// - No attempt to parse output here; separation of concerns keeps runner focused on orchestration.
    /// </summary>
    public sealed class FlashDeconvRunner
    {
        private readonly string _exe;

        /// <summary>
        /// Creates a runner instance.
        /// 
        /// Parameter semantics (WHAT / WHY):
        /// overridePath:
        ///   Explicit path to the executable (wins over all other resolution strategies). Useful for A/B testing versions.
        /// preflight:
        ///   true  => Run a lightweight help command (-h or custom arg) to ensure binary loads and dependencies resolve.
        ///   false => Skip preflight entirely (reduces startup latency but risks late failure).
        /// preflightTimeoutMs:
        ///   >0    => Maximum time (ms) to wait for preflight before treating it as failure (prevents indefinite hang).
        ///   null or <=0 => Infinite wait (useful in cold container starts where dynamic linking is slow).
        /// preflightArg:
        ///   The argument used during preflight. Some versions may require "--help" instead of "-h"; making this configurable
        ///   avoids recompilation. Null falls back to environment variable or "-h".
        /// 
        /// Environment variable overrides (WHY they matter):
        ///   FLASHDECONV_SKIP_PREFLIGHT          => Force skip regardless of constructor 'preflight' value (CI toggles).
        ///   FLASHDECONV_PREFLIGHT_TIMEOUT_MS    => Adjust preflight timeout externally (no code redeploy).
        ///   FLASHDECONV_PREFLIGHT_ARG           => Provide alternative help flag or blank (some builds behave differently).
        /// </summary>
        public FlashDeconvRunner(
            string? overridePath = null,
            bool preflight = true,
            int? preflightTimeoutMs = 30000,
            string? preflightArg = null)
        {
            // Resolve the executable path (throws if not found) – WHY: fail early for more actionable error reporting.
            _exe = FlashDeconvLocator.FindExecutable(overridePath)
                   ?? throw new FileNotFoundException("FLASHDeconv executable not found (env FLASHDECONV_PATH or vendored runtime).");

            // Environment can enforce skip (WHY: CI / stress tests may prefer speed over early validation).
            var skip = Environment.GetEnvironmentVariable("FLASHDECONV_SKIP_PREFLIGHT");
            if (string.Equals(skip, "1", StringComparison.OrdinalIgnoreCase) ||
                string.Equals(skip, "true", StringComparison.OrdinalIgnoreCase))
                preflight = false;

            if (!preflight)
                return;

            // Allow environment to override the timeout (WHY: flaky shared infra might need longer warmups).
            if (Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_TIMEOUT_MS") is string tov &&
                int.TryParse(tov, out var envTimeout))
                preflightTimeoutMs = envTimeout;

            // Provide overridable help argument (WHY: some distributions require "--help" or a blank invocation).
            preflightArg ??= Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_ARG") ?? "-h";

            // Execute the preflight; on failure gather enough context for straightforward remediation.
            if (!FlashDeconvLocator.Preflight(
                    _exe,
                    out var so,
                    out var se,
                    out var ec,
                    preflightTimeoutMs,
                    preflightArg))
            {
                var msg = new StringBuilder()
                    .AppendLine("FLASHDeconv preflight failed.")
                    .AppendLine($"  Path   : {_exe}")
                    .AppendLine($"  Arg    : {preflightArg}")
                    .AppendLine($"  Exit   : {ec}")
                    .AppendLine($"  Timeout: {(preflightTimeoutMs is null or <= 0 ? "None (infinite)" : preflightTimeoutMs + " ms")}")
                    .AppendLine($"  StdOut : {(string.IsNullOrWhiteSpace(so) ? "<empty>" : Trim(so, 400))}")
                    .AppendLine($"  StdErr : {(string.IsNullOrWhiteSpace(se) ? "<empty>" : Trim(se, 400))}")
                    .AppendLine("Hints:")
                    .AppendLine("  - Skip preflight: set FLASHDECONV_SKIP_PREFLIGHT=1 or use new FlashDeconvRunner(preflight:false)")
                    .AppendLine("  - Increase timeout: set FLASHDECONV_PREFLIGHT_TIMEOUT_MS=60000")
                    .AppendLine("  - Try different help arg: set FLASHDECONV_PREFLIGHT_ARG=--help or empty")
                    .AppendLine("  - If help hangs, skip preflight and run directly.")
                    .ToString();
                throw new InvalidOperationException(msg);
            }
        }

        /// <summary>
        /// Synchronous wrapper over the asynchronous run. WHY:
        /// - Many call sites (e.g., legacy test frameworks) are synchronous; this avoids boilerplate GetAwaiter code.
        /// - Keeps single implementation surface (RunAsync) to avoid duplication.
        /// </summary>
        public (int ExitCode, string StdOut, string StdErr) Run(
            string inputMzml, string outputTsv, string extraArgs = "", int timeoutMs = 0)
            => RunAsync(inputMzml, outputTsv, extraArgs, timeoutMs, CancellationToken.None).GetAwaiter().GetResult();

        /// <summary>
        /// Executes FLASHDeconv asynchronously.
        /// 
        /// WHAT:
        /// - Spawns the process with mandatory -in / -out arguments plus optional extra args appended verbatim.
        /// - Captures stdout and stderr concurrently.
        /// - Supports cooperative cancellation (via provided token) and a simple wall-clock timeout (timeoutMs).
        /// 
        /// WHY (timeout design):
        /// - A positive timeout ensures runaway processes (e.g., hung native calls) do not stall test runs indefinitely.
        /// - A non-positive timeout disables enforcement (caller explicitly accepts indefinite runtime).
        /// 
        /// Failure semantics:
        /// - Throws FileNotFoundException if input is missing (fast path).
        /// - Throws TimeoutException if enforced timeout elapses (with best-effort process kill).
        /// - Does NOT throw on non-zero exit; instead returns exit code to let caller decide (unlike the simpler FlashDeconv wrapper).
        ///   WHY: This runner aims to be more flexible for scenarios where partial failures are inspected programmatically.
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

            // Ensure output directory exists (WHY: tool will not create deep path structures).
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
                UseShellExecute = false,             // WHY: required for redirection; avoids shell quoting ambiguities
                RedirectStandardOutput = true,       // WHY: capture diagnostics & potential version info
                RedirectStandardError = true,        // WHY: capture warnings/errors from native layers
                CreateNoWindow = true                // WHY: suppress window in GUI / CI environments
            };

            using var p = new Process { StartInfo = psi };
            var so = new StringBuilder();
            var se = new StringBuilder();

            // Event handlers accumulate output lines; asynchronous pattern prevents deadlocks if tool writes heavily.
            p.OutputDataReceived += (_, e) => { if (e.Data != null) so.AppendLine(e.Data); };
            p.ErrorDataReceived += (_, e) => { if (e.Data != null) se.AppendLine(e.Data); };

            if (!p.Start())
                throw new InvalidOperationException("Failed to start FLASHDeconv process.");

            p.BeginOutputReadLine();
            p.BeginErrorReadLine();

            // Linked CTS combines caller cancellation + timeout (if any) for unified handling.
            using var linkedCts = CancellationTokenSource.CreateLinkedTokenSource(ct);
            if (timeoutMs > 0)
                linkedCts.CancelAfter(timeoutMs);

            try
            {
                // Offload blocking WaitForExit to thread pool so cancellation token can interrupt.
                await Task.Run(() => p.WaitForExit(), linkedCts.Token);
            }
            catch (OperationCanceledException)
            {
                // Distinguish between timeout vs. external cancellation only by context (both arrive here).
                if (!p.HasExited)
                {
                    try { p.Kill(true); } catch { /* swallow kill exceptions (process might exit between checks) */ }
                    if (timeoutMs > 0 && !ct.IsCancellationRequested)
                        throw new TimeoutException($"FLASHDeconv run timed out after {timeoutMs} ms");
                    throw; // propagate caller cancellation
                }
            }

            return (p.ExitCode, so.ToString(), se.ToString());
        }

        /// <summary>
        /// Utility: trims long stdout/stderr segments for concise diagnostics (WHY: prevents log flooding in CI).
        /// </summary>
        private static string Trim(string s, int max) =>
            s.Length <= max ? s : s.Substring(0, max) + "...";
    }
}