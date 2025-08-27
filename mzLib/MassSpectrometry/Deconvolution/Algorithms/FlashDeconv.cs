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
    /// Minimal wrapper for the external FLASHDeconv executable.
    /// (Simplified: no dependency probing, no temp dirs, no advanced error handling.)
    /// </summary>
    public sealed class FlashDeconv : DeconvolutionAlgorithm
    {
        private static readonly string[] CandidateExeNames =
        {
            "FLASHDeconv.exe",
            "FlashDeconv.exe"
        };

        private readonly string _exePath;

        public string ExecutablePath => _exePath;

        public FlashDeconv(DeconvolutionParameters deconParameters, string? exePath = null)
            : base(deconParameters)
        {
            _exePath = ResolveExecutable(exePath);
            if (!File.Exists(_exePath))
            {
                throw new FileNotFoundException(
                    $"FLASHDeconv executable not found at '{_exePath}'. " +
                    "Provide full path or set FLASHDECONV_PATH.", _exePath);
            }
        }

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

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            throw new NotImplementedException("FLASHDeconv result parsing to IsotopicEnvelope not implemented.");
        }

        /// <summary>
        /// Run FLASHDeconv with minimal arguments.
        /// Caller is responsible for selecting appropriate extra arguments (may be blank).
        /// </summary>
        public void Run(string inputFile, string outputFile, string? extraArgs = null)
        {
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
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                CreateNoWindow = true,
                WorkingDirectory = Path.GetDirectoryName(_exePath) ?? Environment.CurrentDirectory
            };

            using var proc = new Process { StartInfo = psi };
            var stdout = new StringBuilder();
            var stderr = new StringBuilder();

            proc.OutputDataReceived += (_, e) => { if (e.Data != null) stdout.AppendLine(e.Data); };
            proc.ErrorDataReceived  += (_, e) => { if (e.Data != null) stderr.AppendLine(e.Data); };

            if (!proc.Start())
                throw new InvalidOperationException("Failed to start FLASHDeconv.");

            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();

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
