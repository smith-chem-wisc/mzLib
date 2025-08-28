using System;
using System.IO;
using System.Diagnostics;

namespace MassSpectrometry.FlashDeconvRuntime
{
    /// <summary>
    /// Utility class to locate and optionally "preflight" (sanity check) the external FLASHDeconv executable.
    /// 
    /// WHY this exists:
    /// - Consumers (tests / runners / pipelines) often need a single, predictable place to resolve the binary.
    /// - Hides the search strategy so future location adjustments do not ripple through the codebase.
    /// - Provides a lightweight dependency probe (Preflight) that detects missing DLLs or environment issues
    ///   before a long analysis workflow starts (failing fast improves user feedback).
    /// 
    /// WHAT it does:
    /// - Resolve a candidate path using (explicit override ? environment variable ? conventional relative path).
    /// - Optionally run a short help invocation to ensure the process launches and dynamic dependencies load.
    /// 
    /// Design notes:
    /// - Returns null instead of throwing in FindExecutable: calling code can decide whether absence is fatal.
    /// - Preflight returns success if exit code is 0 or 1 because some tools use 1 for "printed help" rather
    ///   than a failure; accommodating that reduces false negatives.
    /// - Timeout is optional and can be disabled (null or <= 0) for scenarios where environment can be slow
    ///   to JIT / cache (e.g. first cold run in a container or network share).
    /// </summary>
    public static class FlashDeconvLocator
    {
        /// <summary>
        /// Attempts to find the FLASHDeconv executable.
        /// Resolution order (WHY this order):
        /// 1. Explicit override path (most specific / user intent wins).
        /// 2. Environment variable FLASHDECONV_PATH (CI or user profile configuration).
        /// 3. Conventional relative install: {AppContext.BaseDirectory}/External/OpenMS/FLASHDeconv.exe
        ///    (makes packaging deterministic; allows bundling inside published output).
        /// 
        /// Returns:
        /// - Full path string if found.
        /// - Null if no valid file exists (caller chooses how to handle).
        /// 
        /// NOTE: We do not validate executability (permissions); we only check existence to keep this fast.
        /// </summary>
        public static string? FindExecutable(string? overridePath = null)
        {
            // 1. Honor explicit direct path if supplied and exists.
            if (!string.IsNullOrWhiteSpace(overridePath) && File.Exists(overridePath))
                return overridePath;

            // 2. Check environment variable (enables user / CI configuration without code changes).
            var env = Environment.GetEnvironmentVariable("FLASHDECONV_PATH");
            if (!string.IsNullOrWhiteSpace(env) && File.Exists(env))
                return env;

            // 3. Fallback heuristic: packaged alongside application under External/OpenMS.
            var baseDir = AppContext.BaseDirectory;
            var candidate = Path.Combine(baseDir, "External", "OpenMS", "FLASHDeconv.exe");
            if (File.Exists(candidate))
                return candidate;

            // 4. Nothing found.
            return null;
        }

        /// <summary>
        /// Performs a "preflight" sanity check by invoking FLASHDeconv with a help argument (default -h).
        /// 
        /// WHY:
        /// - Detect missing native dependencies (DLL load failures appear immediately in stderr).
        /// - Confirm the binary is at least launchable before performing long runs.
        /// - Provide early diagnostics (stdout / stderr capture) to surface version or path issues.
        /// 
        /// Behavior details:
        /// - Treats exit code 0 OR 1 as success because some tools return 1 after printing help.
        /// - Captures both standard output & error even on success for optional version parsing.
        /// - If preflightTimeoutMs <= 0 or null ? wait indefinitely (caller explicitly opts-out of timeout).
        /// - On timeout: attempts to terminate process and returns false with a synthetic exit code -999.
        /// - On start failure returns false with exit = -998 providing a distinct sentinel.
        /// 
        /// Returns:
        /// - true if launch succeeded and exit code ? {0,1}.
        /// - false otherwise (stdout/stderr/exit still populated for diagnostics).
        /// </summary>
        /// <param name="exePath">Full path to FLASHDeconv.</param>
        /// <param name="stdout">Captured standard output (empty string on early failures).</param>
        /// <param name="stderr">Captured standard error (includes timeout message if triggered).</param>
        /// <param name="exit">Process exit code or sentinel (-998 start fail, -999 timeout).</param>
        /// <param name="preflightTimeoutMs">Max wait in ms. null or &lt;=0 disables timeout.</param>
        /// <param name="helpArg">Help argument; pass null/empty to run without arguments.</param>
        public static bool Preflight(
            string exePath,
            out string stdout,
            out string stderr,
            out int exit,
            int? preflightTimeoutMs = 30000,
            string? helpArg = "-h")
        {
            stdout = "";
            stderr = "";
            exit = -999; // default sentinel until replaced

            var psi = new ProcessStartInfo
            {
                FileName = exePath,
                Arguments = string.IsNullOrWhiteSpace(helpArg) ? "" : helpArg,
                UseShellExecute = false,              // WHY: Needed for redirection & non-shell context
                RedirectStandardError = true,         // WHY: Capture dependency load failures / version info
                RedirectStandardOutput = true,        // WHY: For version parsing & test validation
                CreateNoWindow = true                 // WHY: Prevent stray console windows in GUI/CI contexts
            };

            using var p = Process.Start(psi);
            if (p == null)
            {
                stderr = "Failed to start process.";
                exit = -998; // distinct sentinel: start failure
                return false;
            }

            // Determine wait strategy. We do not call ReadToEnd before WaitForExit to avoid deadlocks.
            if (preflightTimeoutMs is null || preflightTimeoutMs <= 0)
            {
                // Indefinite wait: explicit choice by caller to allow long cold starts.
                p.WaitForExit();
            }
            else if (!p.WaitForExit(preflightTimeoutMs.Value))
            {
                // Timeout: attempt best-effort kill to avoid orphan processes.
                try { p.Kill(true); } catch { /* ignore kill exceptions */ }
                stderr = $"Timeout ({preflightTimeoutMs} ms)";
                exit = -999;
                return false;
            }

            // Only after exit do we consume streams to completion.
            exit = p.ExitCode;
            stdout = p.StandardOutput.ReadToEnd();
            stderr = p.StandardError.ReadToEnd();

            // Accept both 0 (normal) and 1 (some tools treat help as 'informational non-zero').
            return exit == 0 || exit == 1;
        }
    }
}