using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Test.FileReadingTests;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestDeconvolution
    {
        #region Old Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)] //This is a lesser abundant charge state envelope at the low mz end
        [TestCase(740.372202090153, 19, 108419280, 740.37)] //This is the most abundant charge state envelope
        [TestCase(1081.385183, 13, 35454636,
            1081.385)] //This is a lesser abundant charge state envelope at the high mz end
        public void TestDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess,
            double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];

            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617),
                "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz,
                selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            List<IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(spectrum, 1, 60, 4, 3).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            Assert.That(monoIsotopicMasses[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK",
            "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD",
            "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein test1 = new Protein(peptide, "Accession");
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications pw = new PeptideWithSetModifications(test1, d, 1, test1.Length,
                CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double m = pw.MostAbundantMonoisotopicMass.ToMz(charge);

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            var reader = MsDataFileReader.GetDataFile(singleScan);
            reader.LoadAllStaticData();

            List<MsDataScan> singlescan = reader.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());
            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;

            //check assigned correctly
            List<IsotopicEnvelope> lie2 = singlespec.Deconvolute(singleRange, minAssumedChargeState,
                maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMass / charge, Is.EqualTo(m).Within(0.1));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = singlespec.Deconvolute(singleRange, minAssumedChargeState,
                maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
            Assert.That(lie2.Select(p => p.MostAbundantObservedIsotopicMass),
                Is.EqualTo(lie3.Select(p => p.MostAbundantObservedIsotopicMass)).Within(.0005));
        }

        #endregion

        #region Classic Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)] //This is a lesser abundant charge state envelope at the low mz end
        [TestCase(740.372202090153, 19, 108419280, 740.37)] //This is the most abundant charge state envelope
        [TestCase(1081.385183, 13, 35454636,
            1081.385)] //This is a lesser abundant charge state envelope at the high mz end
        public void TestClassicDeconvolutionProteoformMultiChargeState(double selectedIonMz,
            int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];

            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617),
                "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz,
                selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            // The ones marked 2 are for checking an overload method

            DeconvolutionParameters deconParameters = new ClassicDeconvolutionParameters(1, 60, 4, 3);

            List<IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(scan, deconParameters).ToList();
            List<IsotopicEnvelope> isolatedMasses2 =
                scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, deconParameters).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();
            List<double> monoIsotopicMasses2 = isolatedMasses2.Select(m => m.MonoisotopicMass).ToList();

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            Assert.That(monoIsotopicMasses[0], Is.EqualTo(14037.926829).Within(.0005));
            Assert.That(monoIsotopicMasses2[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK",
            "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD",
            "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckClassicGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein test1 = new Protein(peptide, "Accession");
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications pw = new PeptideWithSetModifications(test1, d, 1, test1.Length,
                CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double m = pw.MostAbundantMonoisotopicMass.ToMz(charge);

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            Mzml singleMZML = (Mzml)MsDataFileReader.GetDataFile(singleScan).LoadAllStaticData();

            List<MsDataScan> singlescan = singleMZML.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());
            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;

            DeconvolutionParameters deconParameters =
                new ClassicDeconvolutionParameters(minAssumedChargeState, maxAssumedChargeState,
                    deconvolutionTolerancePpm,
                    intensityRatioLimit);

            //check assigned correctly

            List<IsotopicEnvelope> lie2 = Deconvoluter.Deconvolute(singlespec, deconParameters, singleRange).ToList();
            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMass / charge, Is.EqualTo(m).Within(0.1));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = Deconvoluter.Deconvolute(singlespec, deconParameters, singleRange).ToList();
            Assert.That(lie2.Select(p => p.MostAbundantObservedIsotopicMass), Is.EqualTo(lie3.Select(p => p.MostAbundantObservedIsotopicMass)));
        }

        [Test]
        [TestCase(373.85, -5, 1874.28)] // GUAGUC -5
        [TestCase(467.57, -4, 1874.28)] // GUAGUC -4
        [TestCase(623.75, -3, 1874.28)] // GUAGUC -3
        [TestCase(936.13, -2, 1874.28)] // GUAGUC -2
        [TestCase(473.05, -4, 1896.26)] // GUAGUC +Na -H -4
        [TestCase(631.07, -3, 1896.26)] // GUAGUC +Na -H -3
        [TestCase(947.121, -2, 1896.26)] // GUAGUC +Na -H -2
        public void TestNegativeModeClassicDeconvolution_Averagine(double expectedMz, int expectedCharge, double expectedMonoMass)
        {
            // get scan
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var scan = MsDataFileReader.GetDataFile(filePath).GetAllScansList().First();
            var tolerance = new PpmTolerance(20);

            // set up deconvolution
            DeconvolutionParameters deconParams = new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative);

            List<IsotopicEnvelope> deconvolutionResults = Deconvoluter.Deconvolute(scan, deconParams).ToList();
            // ensure each expected result is found, with correct mz, charge, and monoisotopic mass
            var resultsWithPeakOfInterest = deconvolutionResults.FirstOrDefault(envelope =>
                envelope.Peaks.Any(peak => tolerance.Within(peak.mz, expectedMz)));
            if (resultsWithPeakOfInterest is null) Assert.Fail();

            Assert.That(expectedMonoMass, Is.EqualTo(resultsWithPeakOfInterest.MonoisotopicMass).Within(0.01));
            Assert.That(expectedCharge, Is.EqualTo(resultsWithPeakOfInterest.Charge));
        }

        [Test]
        [TestCase(373.85, -5, 1874.28)] // GUAGUC -5
        [TestCase(467.57, -4, 1874.28)] // GUAGUC -4
        [TestCase(623.75, -3, 1874.28)] // GUAGUC -3
        [TestCase(936.13, -2, 1874.28)] // GUAGUC -2
        [TestCase(473.05, -4, 1896.26)] // GUAGUC +Na -H -4
        [TestCase(631.07, -3, 1896.26)] // GUAGUC +Na -H -3
        [TestCase(947.121, -2, 1896.26)] // GUAGUC +Na -H -2
        public void TestNegativeModeClassicDeconvolution_Averagtide(double expectedMz, int expectedCharge, double expectedMonoMass)
        {
            // get scan
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var scan = MsDataFileReader.GetDataFile(filePath).GetAllScansList().First();
            var tolerance = new PpmTolerance(20);

            // set up deconvolution
            DeconvolutionParameters deconParams = new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative, new OxyriboAveragine());

            List<IsotopicEnvelope> deconvolutionResults = Deconvoluter.Deconvolute(scan, deconParams).ToList();
            // ensure each expected result is found, with correct mz, charge, and monoisotopic mass
            var resultsWithPeakOfInterest = deconvolutionResults.FirstOrDefault(envelope =>
                envelope.Peaks.Any(peak => tolerance.Within(peak.mz, expectedMz)));
            if (resultsWithPeakOfInterest is null) Assert.Fail();

            Assert.That(expectedMonoMass, Is.EqualTo(resultsWithPeakOfInterest.MonoisotopicMass).Within(0.01));
            Assert.That(expectedCharge, Is.EqualTo(resultsWithPeakOfInterest.Charge));
        }

        #endregion

        #region IsoDec Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]//This is a lesser abundant charge state envelope at the low mz end
        [TestCase(740.372202090153, 19, 108419280, 740.37)]//This is the most abundant charge state envelope
        [TestCase(1081.385183, 13, 35454636, 1081.385)]//This is a lesser abundant charge state envelope at the high mz end
        public void TestIsoDecDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];

            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617), "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            // The ones marked 2 are for checking an overload method

            DeconvolutionParameters deconParameters = new IsoDecDeconvolutionParameters();

            IsoDecAlgorithm alg  = new IsoDecAlgorithm(deconParameters);
            List<IsotopicEnvelope> allMasses = alg.Deconvolute(scan.MassSpectrum, new MzRange((double)scan.MassSpectrum.FirstX, (double)scan.MassSpectrum.LastX)).ToList();

            List<IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(scan, deconParameters).ToList();
            List<IsotopicEnvelope> isolatedMasses2 = scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, deconParameters).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();
            List<double> monoIsotopicMasses2 = isolatedMasses2.Select(m => m.MonoisotopicMass).ToList();
            Assert.That(monoIsotopicMasses2.Count, Is.EqualTo(monoIsotopicMasses.Count));

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            //Updated the tolerance on this test to be 5 ppm (which felt reasonable to me? --JGP)
            double ppmwidth = (14037.926829 / 1e6) * 5;
            bool isAnyEqual1 = monoIsotopicMasses.Any(m =>  m >= 14037.926829 - ppmwidth && m <= 14037.926826 + ppmwidth);
            bool isAnyEqual2 = monoIsotopicMasses2.Any(m => m >= 14037.926829 - ppmwidth && m <= 14037.926826 + ppmwidth);
            Assert.That(isAnyEqual1);
            Assert.That(isAnyEqual2);
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK", "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD", "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckIsoDecGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein test1 = new Protein(peptide, "Accession");
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications pw = new PeptideWithSetModifications(test1, d, 1, test1.Length, CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double m = pw.MostAbundantMonoisotopicMass.ToMz(charge);

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            Mzml singleMZML = (Mzml)MsDataFileReader.GetDataFile(singleScan).LoadAllStaticData();

            List<MsDataScan> singlescan = singleMZML.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());
            DeconvolutionParameters deconParameters = new IsoDecDeconvolutionParameters();

            //check assigned correctly
            List<IsotopicEnvelope> lie2 = Deconvoluter.Deconvolute(singlespec, deconParameters, singleRange).ToList();
            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMass / charge, Is.EqualTo(m).Within(0.1));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = Deconvoluter.Deconvolute(singlespec, deconParameters, singleRange).ToList();
            Assert.That(lie2.Select(p => p.MostAbundantObservedIsotopicMass), Is.EqualTo(lie3.Select(p => p.MostAbundantObservedIsotopicMass)));
        }

        [Test]
        [TestCase(373.85, -5, 1874.28)] // GUAGUC -5
        [TestCase(936.13, -2, 1874.28)] // GUAGUC -2
        [TestCase(473.05, -4, 1896.26)] // GUAGUC +Na -H -4
        [TestCase(631.07, -3, 1896.26)] // GUAGUC +Na -H -3
        [TestCase(947.121, -2, 1896.26)] // GUAGUC +Na -H -2
        public void TestNegativeModeIsoDecDeconvolution(double expectedMz, int expectedCharge, double expectedMonoMass)
        {
            // get scan
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var scan = MsDataFileReader.GetDataFile(filePath).GetAllScansList().First();
            var tolerance = new PpmTolerance(20);

            // set up deconvolution
            DeconvolutionParameters deconParams = new IsoDecDeconvolutionParameters(Polarity.Negative);

            List<IsotopicEnvelope> deconvolutionResults = Deconvoluter.Deconvolute(scan, deconParams).ToList();
            // ensure each expected result is found, with correct mz, charge, and monoisotopic mass
            /*var resultsWithPeakOfInterest = deconvolutionResults.FirstOrDefault(envelope => tolerance.Within(envelope.MonoisotopicMassString, expectedMonoMass));*/
            var resultsWithPeakOfInterest = deconvolutionResults.FirstOrDefault(envelope =>
                envelope.Peaks.Any(peak => tolerance.Within(peak.mz, expectedMz)));
            if (resultsWithPeakOfInterest is null) Assert.Fail();

            Assert.That(expectedMonoMass, Is.EqualTo(resultsWithPeakOfInterest.MonoisotopicMass).Within(0.01));
            Assert.That(expectedCharge, Is.EqualTo(resultsWithPeakOfInterest.Charge));
        }

        #endregion

        [Test]
        public static void TestExampleNewDeconvolutionInDeconvoluter()
        {
            DeconvolutionParameters deconParams = new ExampleNewDeconvolutionParametersTemplate(1, 60);
            var dataFile = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DataFiles", "GUACUG_NegativeMode_Sliced.mzML"));
            dataFile.InitiateDynamicConnection();
            var scan = dataFile.GetOneBasedScanFromDynamicConnection(726);
            var spectrum = scan.MassSpectrum;
            dataFile.CloseDynamicConnection();

            // test switch statements in Deconvoluter
            Assert.Throws<NotImplementedException>(() => _ = Deconvoluter.Deconvolute(spectrum, deconParams).ToList());
            Assert.Throws<NotImplementedException>(() => _ = Deconvoluter.Deconvolute(scan, deconParams).ToList());

            // test default exceptions in deconvoluter
            var badEnumValue = (DeconvolutionType)Int32.MaxValue;
            deconParams.GetType().GetProperty("DeconvolutionType")!.SetValue(deconParams, badEnumValue);
            Assert.Throws<MzLibException>(() => _ = Deconvoluter.Deconvolute(spectrum, deconParams).ToList());
            Assert.Throws<MzLibException>(() => _ = Deconvoluter.Deconvolute(scan, deconParams).ToList());
        }


        [Test]
        public static void Test_MsDataScan_GetIsolatedMassesAndCharges()
        {
            // get scan
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var dataFile = MsDataFileReader.GetDataFile(filePath);
            var precursorScan = dataFile.GetOneBasedScan(1);
            var fragmentationScan = dataFile.GetOneBasedScan(2);

            // set up deconvolution
            DeconvolutionParameters deconParams = new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative);

            // get isolated masses and charges on an MS1 scan. This means the isolation window is null.
            var ms1Result = precursorScan.GetIsolatedMassesAndCharges(precursorScan.MassSpectrum, deconParams).ToList();
            Assert.That(ms1Result.Count, Is.EqualTo(0));
            ms1Result = precursorScan.GetIsolatedMassesAndCharges(precursorScan, deconParams).ToList();
            Assert.That(ms1Result.Count, Is.EqualTo(0));

            // get isolated masses and charges on an MS2 scan. This should work correctly
            var ms2Result = fragmentationScan.GetIsolatedMassesAndCharges(precursorScan.MassSpectrum, deconParams)
                .ToList();
            Assert.That(ms2Result.Count, Is.EqualTo(1));
            ms2Result = fragmentationScan.GetIsolatedMassesAndCharges(precursorScan, deconParams).ToList();
            Assert.That(ms2Result.Count, Is.EqualTo(1));
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_AllInRange()
        {
            // Arrange
            var xArray = new[] { 260.774188159546, 391.660998843979 };
            var yArray = new[] { 1000.0, 1.0 };
            var charges = new[] { 1, 1 };
            var spectrum = new NeutralMassSpectrum(xArray, yArray, charges, false);
            var deconvolutionParameters = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var rangeToGetPeaksFrom = new MzRange(260.0, 400.0);

            // Act
            var result = Deconvoluter.Deconvolute(spectrum, deconvolutionParameters, rangeToGetPeaksFrom).ToList();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.InstanceOf<IEnumerable<IsotopicEnvelope>>());
            Assert.That(result.Count(), Is.EqualTo(2));

            for (int i = 0; i < result.Count(); i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(xArray[i]));
                Assert.That(result[i].TotalIntensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Peaks.Count, Is.EqualTo(1));
                Assert.That(result[i].Peaks.First().mz, Is.EqualTo(xArray[i].ToMz(charges[i])));
                Assert.That(result[i].Peaks.First().intensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Charge, Is.EqualTo(charges[i]));
            }
        }


        [Test]
        public void NeutralMassSpectrum_Deconvolute_AllInRange_Charged()
        {
            // Arrange
            var xArray = new[] { 260.774188159546, 391.660998843979 };
            var yArray = new[] { 1000.0, 1.0 };
            var charges = new[] { 3, 3 };
            var spectrum = new NeutralMassSpectrum(xArray, yArray, charges, false);
            var deconvolutionParameters = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var rangeToGetPeaksFrom = new MzRange(00, 200.0);

            // Act
            var result = Deconvoluter.Deconvolute(spectrum, deconvolutionParameters, rangeToGetPeaksFrom).ToList();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.InstanceOf<IEnumerable<IsotopicEnvelope>>());
            Assert.That(result.Count(), Is.EqualTo(2));

            for (int i = 0; i < result.Count(); i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(xArray[i]));
                Assert.That(result[i].TotalIntensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Peaks.Count, Is.EqualTo(1));
                Assert.That(result[i].Peaks.First().mz, Is.EqualTo(xArray[i].ToMz(charges[i])));
                Assert.That(result[i].Peaks.First().intensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Charge, Is.EqualTo(charges[i]));
            }
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_SomeInRange()
        {
            // Arrange
            var xArray = new[] { 260.774188159546, 391.660998843979 };
            var yArray = new[] { 1000.0, 1.0 };
            var charges = new[] { 1, 1 };
            var spectrum = new NeutralMassSpectrum(xArray, yArray, charges, false);
            var deconvolutionParameters = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var rangeToGetPeaksFrom = new MzRange(260.0, 300.0);

            // Act
            var result = Deconvoluter.Deconvolute(spectrum, deconvolutionParameters, rangeToGetPeaksFrom).ToList();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.InstanceOf<IEnumerable<IsotopicEnvelope>>());
            Assert.That(result.Count(), Is.EqualTo(1));

            for (int i = 0; i < result.Count(); i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(xArray[i]));
                Assert.That(result[i].TotalIntensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Peaks.Count, Is.EqualTo(1));
                Assert.That(result[i].Peaks.First().mz, Is.EqualTo(xArray[i].ToMz(charges[i])));
                Assert.That(result[i].Peaks.First().intensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Charge, Is.EqualTo(charges[i]));
            }
        }

        [Test]
        public void NeutralMassSpectrum_Deconvolute_SomeInRange_Charged()
        {
            // Arrange
            var xArray = new[] { 260.774188159546, 391.660998843979 };
            var yArray = new[] { 1000.0, 1.0 };
            var charges = new[] { 1, 20 };
            var spectrum = new NeutralMassSpectrum(xArray, yArray, charges, false);
            var deconvolutionParameters = new ClassicDeconvolutionParameters(1, 60, 20, 2);
            var rangeToGetPeaksFrom = new MzRange(260.0, 300.0);

            // Act
            var result = Deconvoluter.Deconvolute(spectrum, deconvolutionParameters, rangeToGetPeaksFrom).ToList();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.InstanceOf<IEnumerable<IsotopicEnvelope>>());
            Assert.That(result.Count(), Is.EqualTo(1));

            for (int i = 0; i < result.Count(); i++)
            {
                Assert.That(result[i].MonoisotopicMass, Is.EqualTo(xArray[i]));
                Assert.That(result[i].TotalIntensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Peaks.Count, Is.EqualTo(1));
                Assert.That(result[i].Peaks.First().mz, Is.EqualTo(xArray[i].ToMz(charges[i])));
                Assert.That(result[i].Peaks.First().intensity, Is.EqualTo(yArray[i]));
                Assert.That(result[i].Charge, Is.EqualTo(charges[i]));
            }
        }

    }
}
