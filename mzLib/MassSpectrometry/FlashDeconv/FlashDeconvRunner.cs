using System;
using System.IO;
using System.Text;
using System.Threading;
using System.Diagnostics;
using System.Threading.Tasks;

namespace MassSpectrometry.FlashDeconvRuntime;

public sealed class FlashDeconvRunner
{
    private readonly string _exe;

    /// <summary>
    /// Creates a runner.
    /// preflight:
    ///   true  -> run a quick (-h / --help) load test (with optional timeout)
    ///   false -> skip entirely (no timeout risk)
    /// preflightTimeoutMs:
    ///   >0  -> timeout in ms
    ///   null or <=0 -> wait indefinitely for preflight (no timeout)
    /// Env overrides still respected:
    ///   FLASHDECONV_SKIP_PREFLIGHT (forces skip)
    ///   FLASHDECONV_PREFLIGHT_TIMEOUT_MS
    ///   FLASHDECONV_PREFLIGHT_ARG
    /// </summary>
    public FlashDeconvRunner(
        string? overridePath = null,
        bool preflight = true,
        int? preflightTimeoutMs = 30000,
        string? preflightArg = null)
    {
        _exe = FlashDeconvLocator.FindExecutable(overridePath)
               ?? throw new FileNotFoundException("FLASHDeconv executable not found (env FLASHDECONV_PATH or vendored runtime).");

        // Allow env var to force skip
        var skip = Environment.GetEnvironmentVariable("FLASHDECONV_SKIP_PREFLIGHT");
        if (string.Equals(skip, "1", StringComparison.OrdinalIgnoreCase) ||
            string.Equals(skip, "true", StringComparison.OrdinalIgnoreCase))
            preflight = false;

        if (!preflight)
            return;

        // Env overrides for timeout / arg
        if (Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_TIMEOUT_MS") is string tov &&
            int.TryParse(tov, out var envTimeout))
            preflightTimeoutMs = envTimeout;

        preflightArg ??= Environment.GetEnvironmentVariable("FLASHDECONV_PREFLIGHT_ARG") ?? "-h";

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
                .AppendLine($"  Timeout: {(preflightTimeoutMs is null or <=0 ? "None (infinite)" : preflightTimeoutMs + " ms")}")
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

    public (int ExitCode, string StdOut, string StdErr) Run(
        string inputMzml, string outputTsv, string extraArgs = "", int timeoutMs = 0)
        => RunAsync(inputMzml, outputTsv, extraArgs, timeoutMs, CancellationToken.None).GetAwaiter().GetResult();

    /// <summary>
    /// timeoutMs:
    ///   <=0 => no timeout (let the process run until natural completion).
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
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            CreateNoWindow = true
        };

        using var p = new Process { StartInfo = psi };
        var so = new StringBuilder();
        var se = new StringBuilder();

        p.OutputDataReceived += (_, e) => { if (e.Data != null) so.AppendLine(e.Data); };
        p.ErrorDataReceived += (_, e) => { if (e.Data != null) se.AppendLine(e.Data); };

        if (!p.Start())
            throw new InvalidOperationException("Failed to start FLASHDeconv process.");

        p.BeginOutputReadLine();
        p.BeginErrorReadLine();

        using var linkedCts = CancellationTokenSource.CreateLinkedTokenSource(ct);
        if (timeoutMs > 0)
            linkedCts.CancelAfter(timeoutMs);

        try
        {
            await Task.Run(() => p.WaitForExit(), linkedCts.Token);
        }
        catch (OperationCanceledException)
        {
            if (!p.HasExited)
            {
                try { p.Kill(true); } catch { }
                throw new TimeoutException($"FLASHDeconv run timed out after {timeoutMs} ms");
            }
        }

        return (p.ExitCode, so.ToString(), se.ToString());
    }

    private static string Trim(string s, int max) =>
        s.Length <= max ? s : s.Substring(0, max) + "...";
}