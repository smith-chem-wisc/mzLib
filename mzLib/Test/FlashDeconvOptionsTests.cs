using System;
using System.IO;
using NUnit.Framework;
using MassSpectrometry.FlashDeconvRuntime;      
using System.Diagnostics.CodeAnalysis;


namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class FlashDeconvOptionsTests
    {
        [Test]
        public void FlashDeconvolutionOptionsTest()
        {
            // Re?use the large example path from the existing test (adjust if different on your machine).
            const string inputFile = @"E:\Projects\LVS_TD_Yeast\05-26-17_B7A_yeast_td_fract7_rep1.mzML";

            if (!File.Exists(inputFile))
            {
                Assert.Ignore("Missing input file: " + inputFile);
            }

            // Prepare output folder under the test run directory
            string outDir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "FlashDeconvOptionsRun");
            Directory.CreateDirectory(outDir);

            string outTsv = Path.Combine(outDir, "deconv_features.tsv");
            string outMzml = Path.Combine(outDir, "deconv_spectra.mzML");

            // Build strongly-typed options
            var opts = new FlashDeconvOptions
            {
                InputMzMlPath = inputFile,
                OutputTsvPath = outTsv,
                OutMzmlDeconvolved = outMzml,
                WriteDetail = true,          // ask for detailed peak info
                MaxMsLevel = 3,
                PrecedingMs1Count = 3
            };

            // Narrow mass / charge space for a faster demo run (tune as appropriate)
            opts.Algorithm.MinMass = 500.0;
            opts.Algorithm.MaxMass = 20000.0;
            opts.Algorithm.MinCharge = 1;
            opts.Algorithm.MaxCharge = 60;
            // Tighten isotope cosine for MS1 only (leave others default)
            opts.Algorithm.MinIsotopeCosinePerMsLevel.Clear();
            opts.Algorithm.MinIsotopeCosinePerMsLevel.Add(0.9);

            // Example: set FeatureTracing mass error ppm to use algorithm MS1 tol (leave -1)
            opts.FeatureTracing.MassErrorPpm = -1.0;
            // Request area-based quant (default) explicitly
            opts.FeatureTracing.QuantMethod = FlashDeconvOptions.FeatureTracingGroup.QuantMethodEnum.area;

            // Disable progress to keep stdout minimal
            opts.Common.NoProgress = true;
            // Use more threads if desired
            opts.Common.Threads = Environment.ProcessorCount > 4 ? 4 : 1;

            // Create runner (skip preflight; long startup previously caused timeouts)
            FlashDeconvRunner runner;
            try
            {
                runner = new FlashDeconvRunner(preflight: false);
            }
            catch (FileNotFoundException fnf)
            {
                Assert.Ignore("FLASHDeconv executable not found: " + fnf.Message);
                return;
            }

            // Run (timeoutMs = 0 => no artificial timeout)
            var (exitCode, stdout, stderr) = opts.RunWith(runner, timeoutMs: 0);

            // Basic assertions
            Assert.That(exitCode, Is.EqualTo(0), $"FLASHDeconv failed.\nStdErr:\n{stderr}");
            Assert.That(File.Exists(outTsv), "Expected output feature TSV missing.");
            Assert.That(new FileInfo(outTsv).Length, Is.GreaterThan(0), "Output TSV is empty.");

            if (opts.OutMzmlDeconvolved is not null)
            {
                Assert.That(File.Exists(outMzml), "Expected deconvolved mzML missing.");
            }

            TestContext.WriteLine("FLASHDeconvOptionsTest completed successfully.");
            TestContext.WriteLine("Args: " + opts.ToCommandLineString());
        }
    }
}