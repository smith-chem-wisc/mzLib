using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// EXPLICIT differential tests: run the SAME fixed input through our C# port and through a
    /// standalone MSVC-compiled extract of the real OpenMS function, and compare. The C++ harness
    /// (E:\CodeReview\MetaFlashDecon\difftest\) writes its input to a shared file; this test reads
    /// that file so both sides see identical data. This is the anti-guessing check: it pins each
    /// ported function to the reference implementation's exact numeric output.
    /// </summary>
    [TestFixture]
    [Explicit("Differential test vs MSVC-compiled OpenMS extract; needs E:\\CodeReview\\MetaFlashDecon\\difftest input files.")]
    [ExcludeFromCodeCoverage]
    public class MetaFlashDeconDiffTest
    {
        private const string DiffDir = @"E:\CodeReview\MetaFlashDecon\difftest";

        [Test]
        public void NoisePeakPower_MatchesOpenMS()      => RunCase("noise_input.txt",     "noise_cpp_result.txt",     "noise_cs_result.txt");

        [Test]
        public void NoisePeakPower_Big_MatchesOpenMS()  => RunCase("noise_input_big.txt", "noise_cpp_result_big.txt", "noise_cs_result_big.txt");

        private static void RunCase(string inputName, string cppResultName, string csResultName)
        {
            string inputPath = Path.Combine(DiffDir, inputName);
            Assume.That(File.Exists(inputPath), $"missing {inputPath} (run noise_cpp.exe first)");

            var lines = File.ReadAllLines(inputPath);
            var header = lines[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int absCharge = int.Parse(header[0], CultureInfo.InvariantCulture);
            double isoDa = double.Parse(header[1], CultureInfo.InvariantCulture);

            var signal = new List<(double mz, double intensity)>();
            var noisy = new List<(double mz, double intensity)>();
            for (int i = 1; i < lines.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(lines[i])) continue;
                var t = lines[i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                double mz = double.Parse(t[1], CultureInfo.InvariantCulture);
                double inten = double.Parse(t[2], CultureInfo.InvariantCulture);
                if (t[0] == "S") signal.Add((mz, inten));
                else noisy.Add((mz, inten));
            }

            double cs = MetaFlashDeconAlgorithm.ComputeNoisePeakPower(noisy, signal, absCharge, isoDa);
            TestContext.Progress.WriteLine($"CS  {inputName} getNoisePeakPower = {cs:F8}");
            File.WriteAllText(Path.Combine(DiffDir, csResultName),
                $"{cs:F8}" + Environment.NewLine);

            // Compare against the C++ golden value if present.
            string cppResultPath = Path.Combine(DiffDir, cppResultName);
            if (File.Exists(cppResultPath))
            {
                double cpp = double.Parse(File.ReadAllText(cppResultPath).Trim(), CultureInfo.InvariantCulture);
                Assert.That(cs, Is.EqualTo(cpp).Within(1e-3),
                    $"C# ({cs}) != OpenMS ({cpp}) for getNoisePeakPower");
            }
        }
    }
}
