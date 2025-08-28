using Chemistry;
using Test.FileReadingTests;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;

#nullable enable

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestDeconvolution
    {
        #region Old Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]
        [TestCase(740.372202090153, 19, 108419280, 740.37)]
        [TestCase(1081.385183, 13, 35454636, 1081.385)]
        public void TestDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess,
            double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);
            int n = spectrumLines.Length;
            double[] mzs = new double[n];
            double[] intensities = new double[n];
            for (int i = 0; i < n; i++)
            {
                var p = spectrumLines[i].Split('\t');
                mzs[i] = Convert.ToDouble(p[0], CultureInfo.InvariantCulture);
                intensities[i] = Convert.ToDouble(p[1], CultureInfo.InvariantCulture);
            }

            var spectrum = new MzSpectrum(mzs, intensities, false);
            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617),
                "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz,
                selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);
            var fake = new FakeMsDataFile(Scans);
            var scan = fake.GetAllScansList()[0];

            var isolated = scan.GetIsolatedMassesAndCharges(spectrum, new ClassicDeconvolutionParameters(1, 60, 4, 3))
                .Select(m => m.MonoisotopicMass).ToList();
            Assert.That(isolated[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK",
            "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD",
            "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            var protein = new Protein(peptide, "Acc");
            var pep = new PeptideWithSetModifications(protein, new DigestionParams(), 1, protein.Length,
                CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double expectedMz = pep.MostAbundantMonoisotopicMass.ToMz(charge);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            var scan = reader.GetAllScansList()[0];
            var spectrum = scan.MassSpectrum;
            var range = new MzRange(spectrum.XArray.Min(), spectrum.XArray.Max());

            var dp = new ClassicDeconvolutionParameters(1, 60, 20, 3);
            var envs1 = Deconvoluter.Deconvolute(spectrum, dp, range).ToList();
            var envsCharge = envs1.Where(e => e.Charge == charge).ToList();
            Assert.That(envsCharge[0].MostAbundantObservedIsotopicMass / charge, Is.EqualTo(expectedMz).Within(0.1));

            var envs2 = Deconvoluter.Deconvolute(spectrum, dp, range).ToList();
            Assert.That(envs1.Select(e => e.MostAbundantObservedIsotopicMass),
                Is.EqualTo(envs2.Select(e => e.MostAbundantObservedIsotopicMass)).Within(.0005));
        }

        #endregion

        #region Classic Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]
        [TestCase(740.372202090153, 19, 108419280, 740.37)]
        [TestCase(1081.385183, 13, 35454636, 1081.385)]
        public void TestClassicDeconvolutionProteoformMultiChargeState(double selectedIonMz,
            int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");
            var lines = File.ReadAllLines(path);
            double[] mzs = new double[lines.Length];
            double[] ints = new double[lines.Length];
            for (int i = 0; i < lines.Length; i++)
            {
                var p = lines[i].Split('\t');
                mzs[i] = Convert.ToDouble(p[0], CultureInfo.InvariantCulture);
                ints[i] = Convert.ToDouble(p[1], CultureInfo.InvariantCulture);
            }
            var spectrum = new MzSpectrum(mzs, ints, false);
            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617),
                "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz,
                selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);
            var fake = new FakeMsDataFile(Scans);
            var scan = fake.GetAllScansList()[0];

            var dp = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            var iso1 = scan.GetIsolatedMassesAndCharges(scan, dp).Select(i => i.MonoisotopicMass).ToList();
            var iso2 = scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, dp).Select(i => i.MonoisotopicMass).ToList();
            Assert.That(iso1[0], Is.EqualTo(14037.926829).Within(.0005));
            Assert.That(iso2[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        #endregion

        #region IsoDec / Other Tests (unchanged from existing suite)

        // The remaining deconvolution-related tests from the original file are preserved below.

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]
        [TestCase(740.372202090153, 19, 108419280, 740.37)]
        [TestCase(1081.385183, 13, 35454636, 1081.385)]
        public void TestIsoDecDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");
            var lines = File.ReadAllLines(path);
            double[] mzs = new double[lines.Length];
            double[] ints = new double[lines.Length];
            for (int i = 0; i < lines.Length; i++)
            {
                var p = lines[i].Split('\t');
                mzs[i] = Convert.ToDouble(p[0], CultureInfo.InvariantCulture);
                ints[i] = Convert.ToDouble(p[1], CultureInfo.InvariantCulture);
            }
            var spectrum = new MzSpectrum(mzs, ints, false);
            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617),
                "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz,
                selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);
            var fake = new FakeMsDataFile(Scans);
            var scan = fake.GetAllScansList()[0];

            var dp = new IsoDecDeconvolutionParameters();
            var alg = new IsoDecAlgorithm(dp);
            var allMasses = alg.Deconvolute(scan.MassSpectrum,
                new MzRange((double)scan.MassSpectrum.FirstX!, (double)scan.MassSpectrum.LastX!)).ToList();

            var iso1 = scan.GetIsolatedMassesAndCharges(scan, dp).ToList();
            var iso2 = scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, dp).ToList();
            var mono1 = iso1.Select(m => m.MonoisotopicMass).ToList();
            var mono2 = iso2.Select(m => m.MonoisotopicMass).ToList();
            Assert.That(mono2.Count, Is.EqualTo(mono1.Count));

            double ppmwidth = (14037.926829 / 1e6) * 5;
            bool any1 = mono1.Any(m => m >= 14037.926829 - ppmwidth && m <= 14037.926829 + ppmwidth);
            bool any2 = mono2.Any(m => m >= 14037.926829 - ppmwidth && m <= 14037.926829 + ppmwidth);
            Assert.That(any1);
            Assert.That(any2);
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK", "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD", "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckIsoDecGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            var protein = new Protein(peptide, "Acc");
            var pep = new PeptideWithSetModifications(protein, new DigestionParams(), 1, protein.Length,
                CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double expectedMz = pep.MostAbundantMonoisotopicMass.ToMz(charge);

            string pth = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            var reader = MsDataFileReader.GetDataFile(pth).LoadAllStaticData();
            var scan = ((Mzml)reader).GetAllScansList()[0];
            var spec = scan.MassSpectrum;
            var range = new MzRange(spec.XArray.Min(), spec.XArray.Max());
            var dp = new IsoDecDeconvolutionParameters();

            var envs1 = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            var envCharge = envs1.Where(e => e.Charge == charge).ToList();
            Assert.That(envCharge[0].MostAbundantObservedIsotopicMass / charge, Is.EqualTo(expectedMz).Within(0.1));

            var envs2 = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            Assert.That(envs1.Select(e => e.MostAbundantObservedIsotopicMass),
                Is.EqualTo(envs2.Select(e => e.MostAbundantObservedIsotopicMass)));
        }

        [Test]
        [TestCase(373.85, -5, 1874.28)]
        [TestCase(936.13, -2, 1874.28)]
        [TestCase(473.05, -4, 1896.26)]
        [TestCase(631.07, -3, 1896.26)]
        [TestCase(947.121, -2, 1896.26)]
        public void TestNegativeModeIsoDecDeconvolution(double expectedMz, int expectedCharge, double expectedMonoMass)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var scan = MsDataFileReader.GetDataFile(filePath).GetAllScansList().First();
            var tol = new PpmTolerance(20);
            var dp = new IsoDecDeconvolutionParameters(Polarity.Negative);

            var envs = Deconvoluter.Deconvolute(scan, dp).ToList();
            var match = envs.FirstOrDefault(e => e.Peaks.Any(p => tol.Within(p.mz, expectedMz)));
            Assert.That(match, Is.Not.Null);
            Assert.That(match!.MonoisotopicMass, Is.EqualTo(expectedMonoMass).Within(0.01));
            Assert.That(match.Charge, Is.EqualTo(expectedCharge));
        }

        #endregion

        #region Example + Helper Tests

        [Test]
        public static void TestExampleNewDeconvolutionInDeconvoluter()
        {
            var dp = new ExampleNewDeconvolutionParametersTemplate(1, 60);
            var dataFile = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DataFiles", "GUACUG_NegativeMode_Sliced.mzML"));
            dataFile.InitiateDynamicConnection();
            var scan = dataFile.GetOneBasedScanFromDynamicConnection(726);
            var spec = scan.MassSpectrum;
            dataFile.CloseDynamicConnection();

            Assert.Throws<NotImplementedException>(() => _ = Deconvoluter.Deconvolute(spec, dp).ToList());
            Assert.Throws<NotImplementedException>(() => _ = Deconvoluter.Deconvolute(scan, dp).ToList());

            var bad = (DeconvolutionType)int.MaxValue;
            dp.GetType().GetProperty("DeconvolutionType")!.SetValue(dp, bad);
            Assert.Throws<MzLibException>(() => _ = Deconvoluter.Deconvolute(spec, dp).ToList());
            Assert.Throws<MzLibException>(() => _ = Deconvoluter.Deconvolute(scan, dp).ToList());
        }

        [Test]
        public static void Test_MsDataScan_GetIsolatedMassesAndCharges()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var df = MsDataFileReader.GetDataFile(filePath);
            var ms1 = df.GetOneBasedScan(1);
            var ms2 = df.GetOneBasedScan(2);
            var dp = new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative);

            var r1 = ms1.GetIsolatedMassesAndCharges(ms1.MassSpectrum, dp).ToList();
            Assert.That(r1.Count, Is.EqualTo(0));
            r1 = ms1.GetIsolatedMassesAndCharges(ms1, dp).ToList();
            Assert.That(r1.Count, Is.EqualTo(0));

            var r2 = ms2.GetIsolatedMassesAndCharges(ms1.MassSpectrum, dp).ToList();
            Assert.That(r2.Count, Is.EqualTo(1));
            r2 = ms2.GetIsolatedMassesAndCharges(ms1, dp).ToList();
            Assert.That(r2.Count, Is.EqualTo(1));
        }

        #endregion

        #region NeutralMassSpectrum Tests

        [Test]
        public void NeutralMassSpectrum_Deconvolute_AllInRange()
        {
            var x = new[] { 260.774188159546, 391.660998843979 };
            var y = new[] { 1000.0, 1.0 };
            var c = new[] { 1, 1 };
            var spec = new NeutralMassSpectrum(x, y, c, false);
            var dp = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var range = new MzRange(260.0, 400.0);

            var result = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            Assert.That(result.Count, Is.EqualTo(2));
            for (int i = 0; i < result.Count; i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(x[i]));
                Assert.That(result[i].TotalIntensity, Is.EqualTo(y[i]));
                Assert.That(result[i].Charge, Is.EqualTo(c[i]));
            }
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_AllInRange_Charged()
        {
            var x = new[] { 260.774188159546, 391.660998843979 };
            var y = new[] { 1000.0, 1.0 };
            var c = new[] { 3, 3 };
            var spec = new NeutralMassSpectrum(x, y, c, false);
            var dp = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var range = new MzRange(0, 200.0);

            var result = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            Assert.That(result.Count, Is.EqualTo(2));
            for (int i = 0; i < result.Count; i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(x[i]));
                Assert.That(result[i].Charge, Is.EqualTo(c[i]));
            }
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_SomeInRange()
        {
            var x = new[] { 260.774188159546, 391.660998843979 };
            var y = new[] { 1000.0, 1.0 };
            var c = new[] { 1, 1 };
            var spec = new NeutralMassSpectrum(x, y, c, false);
            var dp = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var range = new MzRange(260.0, 300.0);

            var result = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].MonoisotopicMass, Is.EqualTo(x[0]));
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_SomeInRange_Charged()
        {
            var x = new[] { 260.774188159546, 391.660998843979 };
            var y = new[] { 1000.0, 1.0 };
            var c = new[] { 1, 20 };
            var spec = new NeutralMassSpectrum(x, y, c, false);
            var dp = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var range = new MzRange(260.0, 300.0);

            var result = Deconvoluter.Deconvolute(spec, dp, range).ToList();
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[0].MonoisotopicMass, Is.EqualTo(x[0]));
        }

        #endregion

        #region Simplified FLASHDeconv Integration Test

        [Test]
        public void FlashDeconvolutionTest()
        {
            string inputFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "TopDown", "lvsYeastTopDownSnip.mzML");
            string expectedOutputTsv = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "TopDown", "lvsYeastSnipFlashDeconvExpected.tsv");

            if (!File.Exists(inputFile))
                Assert.Ignore("Missing input file: " + inputFile);
            if (!File.Exists(expectedOutputTsv))
                Assert.Ignore("Missing expected output TSV: " + expectedOutputTsv);

            string tempDir = Path.Combine(Path.GetTempPath(), "FlashDeconvTest_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(tempDir);
            string outputTsv = Path.Combine(tempDir, "deconv.tsv");

            try
            {
                var runner = new MassSpectrometry.FlashDeconvRuntime.FlashDeconvRunner(preflight: false);
                string extraArgs = (Environment.GetEnvironmentVariable("FLASHDECONV_ARGS") ?? "").Trim();
                var result = runner.Run(inputFile, outputTsv, extraArgs, timeoutMs: 0);

                Assert.That(result.ExitCode, Is.EqualTo(0), "FLASHDeconv failed: " + result.StdErr);
                Assert.That(File.Exists(outputTsv), "TSV output not created.");
                Assert.That(new FileInfo(outputTsv).Length, Is.GreaterThan(0), "TSV output empty.");

                CompareTsvIgnoringColumn(expectedOutputTsv, outputTsv,
                    ignoreColumnIndex: 1,        // Ignore FileName column (2nd column, zero-based index = 1)
                    numericAbsTol: 1e-3,
                    numericRelTol: 1e-3);
            }
            finally
            {
                try { Directory.Delete(tempDir, true); } catch { }
            }
        }

        private static void CompareTsvIgnoringColumn(string expectedPath, string actualPath,
            int ignoreColumnIndex, double numericAbsTol, double numericRelTol)
        {
            var expectedLines = File.ReadAllLines(expectedPath)
                .Where(l => !string.IsNullOrWhiteSpace(l)).ToList();
            var actualLines = File.ReadAllLines(actualPath)
                .Where(l => !string.IsNullOrWhiteSpace(l)).ToList();

            Assert.That(actualLines.Count, Is.EqualTo(expectedLines.Count),
                $"Line count mismatch. Expected {expectedLines.Count}, got {actualLines.Count}.");

            if (expectedLines.Count == 0)
                Assert.Fail("Expected TSV is empty.");

            AssertHeadersCompatibleSkip(expectedLines[0], actualLines[0], ignoreColumnIndex);

            for (int i = 1; i < expectedLines.Count; i++)
            {
                var expCols = expectedLines[i].Split('\t');
                var actCols = actualLines[i].Split('\t');

                int maxCols = Math.Max(expCols.Length, actCols.Length);
                if (expCols.Length != actCols.Length)
                {
                    Array.Resize(ref expCols, maxCols);
                    Array.Resize(ref actCols, maxCols);
                }

                for (int c = 0; c < maxCols; c++)
                {
                    if (c == ignoreColumnIndex) // skip FileName column
                        continue;

                    var e = (expCols[c] ?? "").Trim();
                    var a = (actCols[c] ?? "").Trim();

                    if (string.Equals(e, a, StringComparison.Ordinal))
                        continue;

                    if (double.TryParse(e, NumberStyles.Any, CultureInfo.InvariantCulture, out double ev) &&
                        double.TryParse(a, NumberStyles.Any, CultureInfo.InvariantCulture, out double av))
                    {
                        double absDiff = Math.Abs(ev - av);
                        double relDiff = ev != 0 ? absDiff / Math.Abs(ev) : absDiff;
                        if (absDiff <= numericAbsTol || relDiff <= numericRelTol)
                            continue;

                        Assert.Fail($"Numeric mismatch (line {i + 1}, col {c + 1}): expected {ev}, actual {av}, absDiff={absDiff}, relDiff={relDiff}");
                    }
                    else
                    {
                        if (!(string.IsNullOrEmpty(e) && string.IsNullOrEmpty(a)))
                        {
                            Assert.Fail($"Text mismatch (line {i + 1}, col {c + 1}): expected '{e}', actual '{a}'");
                        }
                    }
                }
            }
        }

        private static void AssertHeadersCompatibleSkip(string expectedHeader, string actualHeader, int ignoreColumnIndex)
        {
            if (expectedHeader == actualHeader)
                return;

            var expCols = expectedHeader.Split('\t');
            var actCols = actualHeader.Split('\t');
            int max = Math.Max(expCols.Length, actCols.Length);

            for (int c = 0; c < max; c++)
            {
                if (c == ignoreColumnIndex)
                    continue;

                string e = c < expCols.Length ? expCols[c] : "";
                string a = c < actCols.Length ? actCols[c] : "";

                if (!string.Equals(e, a, StringComparison.Ordinal))
                    Assert.Fail($"Header mismatch at column {c + 1}: expected '{e}', actual '{a}'");
            }
        }
        #endregion
    }
}

