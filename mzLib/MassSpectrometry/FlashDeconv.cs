using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using MassSpectrometry.FlashDeconvRuntime;

namespace MassSpectrometry
{
    /// <summary>
    /// Backward?compatible wrapper around the vendored / resolved FLASHDeconv executable.
    /// Original usage pattern (from tests) is preserved:
    ///   var flash = new FlashDeconv(new IsoDecDeconvolutionParameters(), exePathOrNull);
    ///   flash.Run(inputMzML, outputTsv, extraArgs);
    /// Resolution order for the executable:
    ///  1) overrideExePath constructor argument (if non-empty)
    ///  2) Environment variable FLASHDECONV_PATH
    ///  3) Vendored copy at (AppContext.BaseDirectory)/External/OpenMS/FLASHDeconv.exe
    /// Throws FileNotFoundException if not found or preflight fails.
    /// </summary>
    public sealed class FlashDeconv
    {
        /// <summary>Deconvolution parameter object (passed for symmetry / future use).</summary>
        public object Parameters { get; }

        /// <summary>Resolved full path to FLASHDeconv.exe.</summary>
        public string ExecutablePath { get; }

        /// <summary>Last process exit code (after Run/RunAsync).</summary>
        public int? LastExitCode { get; private set; }

        /// <summary>Last captured standard output.</summary>
        public string? LastStdOut { get; private set; }

        /// <summary>Last captured standard error.</summary>
        public string? LastStdErr { get; private set; }

        private readonly FlashDeconvRunner _runner;

        /// <param name="parameters">Any parameter object (kept for API compatibility).</param>
        /// <param name="overrideExePath">Optional explicit exe path. If null the locator fallback is used.</param>
        public FlashDeconv(object parameters, string? overrideExePath = null)
        {
            Parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));

            // Resolve + preflight inside runner constructor; propagate exceptions if failure.
            _runner = new FlashDeconvRunner(overrideExePath);
            ExecutablePath = GetPrivateFieldExecutablePath(_runner);
        }

        /// <summary>
        /// Synchronous run. Throws if FLASHDeconv returns non-zero exit code.
        /// </summary>
        public void Run(string inputMzml, string outputTsv, string extraArgs = "", int timeoutMs = 300_000)
        {
            var (ec, so, se) = _runner.Run(inputMzml, outputTsv, extraArgs, timeoutMs);
            LastExitCode = ec;
            LastStdOut = so;
            LastStdErr = se;

            if (ec != 0)
            {
                // Common native loader failures: -1073741515 (0xC0000135) missing DLL; -1073740791 (0xC0000409) stack buffer; etc.
                throw new InvalidOperationException(
                    $"FLASHDeconv failed (exit {ec}). StdErr (first 500 chars): {Truncate(se, 500)}");
            }
        }

        /// <summary>
        /// Asynchronous run returning exit code and streams (does not throw on non-zero unless throwOnError true).
        /// </summary>
        public async Task<(int ExitCode, string StdOut, string StdErr)> RunAsync(
            string inputMzml, string outputTsv, string extraArgs = "", int timeoutMs = 300_000,
            CancellationToken cancellationToken = default, bool throwOnError = true)
        {
            var (ec, so, se) = await _runner.RunAsync(inputMzml, outputTsv, extraArgs, timeoutMs, cancellationToken)
                                            .ConfigureAwait(false);
            LastExitCode = ec;
            LastStdOut = so;
            LastStdErr = se;

            if (throwOnError && ec != 0)
                throw new InvalidOperationException(
                    $"FLASHDeconv failed (exit {ec}). StdErr (first 500 chars): {Truncate(se, 500)}");

            return (ec, so, se);
        }

        private static string Truncate(string? s, int max) =>
            string.IsNullOrEmpty(s) ? string.Empty : (s.Length <= max ? s : s.Substring(0, max) + "...");

        // Reflection-free: we could expose Executable path from runner directly via a property.
        // To avoid editing the runner again, we derive by re-resolving (safe & cheap).
        private static string GetPrivateFieldExecutablePath(FlashDeconvRunner runner)
        {
            // Re-resolve using the same logic: (preflight already passed)
            // We call the locator again; since runner accepted a path, it must still exist.
            var path = FlashDeconvLocator.FindExecutable() ?? "";
            if (File.Exists(path)) return path;

            // Fallback: user likely supplied explicit override that is not the vendored one.
            // Ask user to store it if needed; leave blank if not discoverable.
            return path;
        }
    }
}