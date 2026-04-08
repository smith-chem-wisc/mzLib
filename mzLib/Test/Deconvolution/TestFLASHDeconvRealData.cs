using System;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.Deconvolution
{
    /// <summary>
    /// Critical tests against the three real averaged mzML files already in the test suite.
    /// These catch failure modes that synthetic spectra cannot: noise, missing peaks,
    /// overlapping envelopes, and apex-peak misidentification in real data.
    ///
    /// Known masses from StandardDeconvolutionTest.cs (Classic/IsoDec agree on these):
    ///   CytoC   ~ 12,367 Da (charges 9–15 in this dataset)
    ///   HGH     ~ 22,138 Da (charges 11–19)
    ///   Ubiquitin ~ 10,040 Da (charges 8–16)
    ///
    /// Tolerance: 3 Da (the paper states all tools find signals within 3 Da).
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvRealData
    {
        private static readonly string TestDataDir = FindTestDataDir();

        private static string FindTestDataDir()
        {
            var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
            while (dir != null)
            {
                string candidate = Path.Combine(dir.FullName, "Development",
                    "DeconvolutionDevelopment", "TestData");
                if (Directory.Exists(candidate) &&
                    File.Exists(Path.Combine(candidate, "Averaged_221110_CytoOnly.mzML")))
                    return candidate;

                candidate = Path.Combine(dir.FullName, "DeconvolutionDevelopment", "TestData");
                if (Directory.Exists(candidate) &&
                    File.Exists(Path.Combine(candidate, "Averaged_221110_CytoOnly.mzML")))
                    return candidate;

                dir = dir.Parent;
            }
            throw new DirectoryNotFoundException(
                "Could not find DeconvolutionDevelopment/TestData containing the averaged mzML files. " +
                "Searched upward from: " + TestContext.CurrentContext.TestDirectory);
        }

        private static FLASHDeconvolutionParameters Params(int maxZ = 60) =>
            new FLASHDeconvolutionParameters(
                minCharge: 1, maxCharge: maxZ,
                deconvolutionTolerancePpm: 10.0,
                minIsotopicPeakCount: 3,
                minCosineScore: 0.4,
                minMassRange: 50.0,
                maxMassRange: 100_000.0);

        private static MzSpectrum LoadFirstSpectrum(string path)
        {
            var file = MsDataFileReader.GetDataFile(path);
            file.LoadAllStaticData();
            return file.GetAllScansList()[0].MassSpectrum;
        }

        // ══════════════════════════════════════════════════════════════════════
        // 1. Real-data smoke tests — one per file
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void RealData_CytochromeC_FindsMassWithin3Da()
        {
            // Reference: real OpenMS FLASHDeconv output (Averaged_221110_CytoOnly_ms1.tsv).
            // Most abundant peak: 12351.35 Da, cosine=0.997, z=7-21.
            // Note: Classic/IsoDec report ~12367 Da for the same spectrum because they
            // land on a different observable isotope. FLASHDeconv's convention
            // (apex - DiffToMonoisotopic) gives the Averagine formula mono.
            const double knownMass = 12_351.35;
            var spectrum = LoadFirstSpectrum(Path.Combine(TestDataDir, "Averaged_221110_CytoOnly.mzML"));
            var results = Deconvoluter.Deconvolute(spectrum, Params()).ToList();

            Assert.That(results, Is.Not.Empty, "FLASHDeconv should produce at least one envelope on the CytoC spectrum");

            bool found = results.Any(e => Math.Abs(e.MonoisotopicMass - knownMass) <= 3.0);
            Assert.That(found, Is.True,
                $"Expected an envelope within 3 Da of {knownMass} Da (real FLASHDeconv reference). " +
                $"Closest: {results.Min(e => Math.Abs(e.MonoisotopicMass - knownMass)):F2} Da away. " +
                $"All masses (first 10): [{string.Join(", ", results.Take(10).Select(e => $"{e.MonoisotopicMass:F1}"))}]");
        }

        [Test]
        public void RealData_HumanGrowthHormone_FindsMassWithin3Da()
        {
            // Reference: real OpenMS FLASHDeconv output (Averaged_221110_HGHOnly_ms1.tsv).
            // Most abundant peak: 22111.13 Da, cosine=0.996, z=10-20.
            const double knownMass = 22_111.13;
            var spectrum = LoadFirstSpectrum(Path.Combine(TestDataDir, "Averaged_221110_HGHOnly.mzML"));
            var results = Deconvoluter.Deconvolute(spectrum, Params()).ToList();

            Assert.That(results, Is.Not.Empty, "FLASHDeconv should produce at least one envelope on the HGH spectrum");

            bool found = results.Any(e => Math.Abs(e.MonoisotopicMass - knownMass) <= 3.0);
            Assert.That(found, Is.True,
                $"Expected an envelope within 3 Da of {knownMass} Da (real FLASHDeconv reference). " +
                $"Closest: {results.Min(e => Math.Abs(e.MonoisotopicMass - knownMass)):F2} Da away. " +
                $"All masses (first 10): [{string.Join(", ", results.Take(10).Select(e => $"{e.MonoisotopicMass:F1}"))}]");
        }

        [Test]
        public void RealData_Ubiquitin_FindsMassWithin3Da()
        {
            // Reference: real OpenMS FLASHDeconv output (Averaged_221110_UbiqOnly_ms1.tsv).
            // Most abundant peak: 10025.34 Da, cosine=0.997, z=6-13.
            const double knownMass = 10_025.34;
            var spectrum = LoadFirstSpectrum(Path.Combine(TestDataDir, "Averaged_221110_UbiqOnly.mzML"));
            var results = Deconvoluter.Deconvolute(spectrum, Params()).ToList();

            Assert.That(results, Is.Not.Empty, "FLASHDeconv should produce at least one envelope on the Ubiquitin spectrum");

            bool found = results.Any(e => Math.Abs(e.MonoisotopicMass - knownMass) <= 3.0);
            Assert.That(found, Is.True,
                $"Expected an envelope within 3 Da of {knownMass} Da (real FLASHDeconv reference). " +
                $"Closest: {results.Min(e => Math.Abs(e.MonoisotopicMass - knownMass)):F2} Da away. " +
                $"All masses (first 10): [{string.Join(", ", results.Take(10).Select(e => $"{e.MonoisotopicMass:F1}"))}]");
        }

        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void ApexPeak_ZLookup_DoesNotThrowWhenApexRecruitedFromMultipleCharges()
        {
            // The algorithm finds apexEntry by matching apexPeak.mz against recruitedPeaks.
            // When the same raw peak is recruited at multiple charges, its mz appears
            // multiple times. The Where().First() must not throw even if tolerance
            // matching is slightly off due to floating-point.
            // This test uses a spectrum where the apex is in the center of the charge range
            // so every charge state recruits it.
            var p = Params();
            int avgIdx = p.AverageResidueModel.GetMostIntenseMassIndex(12_223.2);
            double apex = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx)[0];
            double diff = p.AverageResidueModel.GetDiffToMonoisotopic(avgIdx);
            double specMono = apex - diff;

            // Build spectrum with dense charge range so apex peak is recruited many times
            double[] masses = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx);
            double[] intens = p.AverageResidueModel.GetAllTheoreticalIntensities(avgIdx);
            var sortedAvg = masses.Zip(intens).OrderBy(x => x.First).ToArray();
            double[] sortedIntens = sortedAvg.Select(x => x.Second).ToArray();

            var mzList = new System.Collections.Generic.List<double>();
            var itList = new System.Collections.Generic.List<double>();
            foreach (int z in new[] { 8, 9, 10, 11, 12, 13, 14, 15 })  // wide charge range
                for (int n = 0; n < sortedAvg.Length; n++)
                {
                    double mz = (specMono + n * Chemistry.Constants.C13MinusC12).ToMz(z);
                    double it = 100_000.0 * sortedIntens[n];
                    if (it < 100.0) continue;
                    mzList.Add(mz);
                    itList.Add(it);
                }
            var pairs = mzList.Zip(itList).OrderBy(x => x.First).ToArray();
            var spectrum = new MzSpectrum(
                pairs.Select(x => x.First).ToArray(),
                pairs.Select(x => x.Second).ToArray(), false);

            // Must not throw, must produce results
            System.Collections.Generic.List<IsotopicEnvelope> results = null;
            Assert.DoesNotThrow(() => results = Deconvoluter.Deconvolute(spectrum, p).ToList(),
                "Deconvolute should not throw when apex peak is recruited across many charges");
            Assert.That(results, Is.Not.Empty,
                "Wide charge range spectrum should yield at least one envelope");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. All envelopes from real data respect invariants
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        [TestCase("Averaged_221110_CytoOnly.mzML")]
        [TestCase("Averaged_221110_HGHOnly.mzML")]
        [TestCase("Averaged_221110_UbiqOnly.mzML")]
        public void RealData_AllEnvelopes_RespectInvariants(string fileName)
        {
            // Every envelope produced on real data must satisfy the same invariants
            // as synthetic data: score in [0,1], mass in range, peaks >= min count,
            // positive charge in positive mode.
            var p = Params();
            var spectrum = LoadFirstSpectrum(Path.Combine(TestDataDir, fileName));
            var results = Deconvoluter.Deconvolute(spectrum, p).ToList();

            foreach (var env in results)
            {
                Assert.That(env.Score, Is.GreaterThanOrEqualTo(p.MinCosineScore)
                                         .And.LessThanOrEqualTo(1.0 + 1e-9),
                    $"Score {env.Score:F4} out of range in {fileName}");
                Assert.That(env.MonoisotopicMass, Is.GreaterThanOrEqualTo(p.MinMassRange)
                                                    .And.LessThanOrEqualTo(p.MaxMassRange),
                    $"Mass {env.MonoisotopicMass:F1} out of range in {fileName}");
                Assert.That(env.Peaks.Count, Is.GreaterThanOrEqualTo(p.MinIsotopicPeakCount),
                    $"Peak count {env.Peaks.Count} below minimum in {fileName}");
                Assert.That(env.Charge, Is.GreaterThan(0),
                    $"Charge {env.Charge} not positive in positive-mode {fileName}");
                Assert.That(env.TotalIntensity, Is.GreaterThan(0.0),
                    $"TotalIntensity {env.TotalIntensity} not positive in {fileName}");
            }
        }
    }
}