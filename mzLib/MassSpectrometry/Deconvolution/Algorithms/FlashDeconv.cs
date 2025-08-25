using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    internal class FlashDeconv : DeconvolutionAlgorithm
    {
        private readonly string exePath;
        internal FlashDeconv(DeconvolutionParameters deconParameters) : base(deconParameters)
        {
        }

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            throw new NotImplementedException();
        }
        /// <summary>
        /// Constructor.
        /// exePath: path to FlashDeconv.exe; if null, assumes it's in the same folder as the app.
        /// deconParameters: parameters for deconvolution algorithm.
        /// </summary>
        public FlashDeconv(DeconvolutionParameters deconParameters, string exePath = null)
            : base(deconParameters)
        {
            if (string.IsNullOrEmpty(exePath))
            {
                // Assume FlashDeconv.exe is in the same folder as the app
                this.exePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "FlashDeconv.exe");
            }
            else
            {
                this.exePath = exePath;
            }

            if (!File.Exists(this.exePath))
                throw new FileNotFoundException("FlashDeconv executable not found", this.exePath);
        }

        /// <summary>
        /// Run FlashDeconv with specified input and output files.
        /// </summary>
        /// <param name="inputFile">Path to input mzML file</param>
        /// <param name="outputFile">Path to output mzML file</param>
        /// <param name="extraArgs">Optional extra arguments for FlashDeconv</param>
        public void Run(string inputFile, string outputFile, string extraArgs = "")
        {
            if (!File.Exists(inputFile))
                throw new FileNotFoundException("Input file not found", inputFile);

            Process proc = new Process();
            proc.StartInfo.FileName = exePath;
            proc.StartInfo.Arguments = $"-in \"{inputFile}\" -out \"{outputFile}\" {extraArgs}";
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.RedirectStandardError = true;

            proc.Start();

            // Capture output and errors
            string stdout = proc.StandardOutput.ReadToEnd();
            string stderr = proc.StandardError.ReadToEnd();

            proc.WaitForExit();

            if (proc.ExitCode != 0)
            {
                throw new Exception($"FlashDeconv exited with code {proc.ExitCode}.\nError output:\n{stderr}");
            }

            Console.WriteLine(stdout);
        }
    }
}
