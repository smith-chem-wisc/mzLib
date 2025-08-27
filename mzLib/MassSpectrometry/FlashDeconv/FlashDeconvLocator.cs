using System;
using System.IO;
using System.Diagnostics;

namespace MassSpectrometry.FlashDeconvRuntime;

public static class FlashDeconvLocator
{
    public static string? FindExecutable(string? overridePath = null)
    {
        if (!string.IsNullOrWhiteSpace(overridePath) && File.Exists(overridePath))
            return overridePath;

        var env = Environment.GetEnvironmentVariable("FLASHDECONV_PATH");
        if (!string.IsNullOrWhiteSpace(env) && File.Exists(env))
            return env;

        var baseDir = AppContext.BaseDirectory;
        var candidate = Path.Combine(baseDir, "External", "OpenMS", "FLASHDeconv.exe");
        if (File.Exists(candidate))
            return candidate;

        return null;
    }

    /// <summary>
    /// Preflight: run help command to ensure dependencies load.
    /// preflightTimeoutMs:
    ///   null or <=0 => wait indefinitely (no timeout).
    /// </summary>
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
        exit = -999;

        var psi = new ProcessStartInfo
        {
            FileName = exePath,
            Arguments = string.IsNullOrWhiteSpace(helpArg) ? "" : helpArg,
            UseShellExecute = false,
            RedirectStandardError = true,
            RedirectStandardOutput = true,
            CreateNoWindow = true
        };

        using var p = Process.Start(psi);
        if (p == null)
        {
            stderr = "Failed to start process.";
            exit = -998;
            return false;
        }

        if (preflightTimeoutMs is null || preflightTimeoutMs <= 0)
        {
            p.WaitForExit(); // infinite wait
        }
        else if (!p.WaitForExit(preflightTimeoutMs.Value))
        {
            try { p.Kill(true); } catch { }
            stderr = $"Timeout ({preflightTimeoutMs} ms)";
            exit = -999;
            return false;
        }

        exit = p.ExitCode;
        stdout = p.StandardOutput.ReadToEnd();
        stderr = p.StandardError.ReadToEnd();

        return exit == 0 || exit == 1;
    }
}