using System;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using FlashLFQ;
using MzLibUtil;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MassSpectrometry;
using Test.FileReadingTests;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class IndexingEngineTests
    {

        private static string _fileToWrite = "indexingEngineTests.mzML";
        private static string _testMzMlFullFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, _fileToWrite);

        [OneTimeSetUp]
        public static void SetUp()
        {

            double intensity = 1e6;

            MsDataScan[] scans = new MsDataScan[9];
            double[] intensityMultipliers = { 1, 2, 3, 4, 5, 4, 3, 2, 1 };

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = new double[] { 500, 500.5, 501, 501.5, 502 };
                double[] intensities = Enumerable.Repeat(intensityMultipliers[s]*intensity, 5).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                _testMzMlFullFilePath, false);
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            File.Delete(_testMzMlFullFilePath);
        }

        [Test]
        public static void TestXicWithDoubleMzPeaks()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.SelectMany(v => new List<double> { v.ToMz(1), (v + 0.0001).ToMz(1) }).ToArray();
                double[] intensities = dist.Intensities.SelectMany(v => new List<double> { v * intensity * intensityMultipliers[s], v * intensity }).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            PeakIndexingEngine indexEngine = PeakIndexingEngine.InitializeIndexingEngine(file1);

            var xic = indexEngine.GetXic(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
            var shiftedXic = indexEngine.GetXic((dist.Masses.First() + 0.0001).ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);

            Assert.That(xic.Count, Is.EqualTo(10));
            Assert.That(shiftedXic.Count, Is.EqualTo(10));
            Assert.That(Math.Abs(xic.First().M - shiftedXic.First().M), Is.GreaterThan(0.00001));

            // Check that the xics don't overlap at all
            xic.AddRange(shiftedXic);
            Assert.That(xic.ToHashSet().Count(), Is.EqualTo(20)); // The intersection of the two should yield twenty different IndexedMzPeaks 
        }

        [Test]
        public static void TestXicStops()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 0, 0, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.Select(v =>  v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                // Turns out you can write a spectrum with a zero intensity peak and it's still recognized by FlashLFQ
                // So, this is a hack so that scans 3 and 4 don't have a peak at the mass We're looking for
                if(s == 2 || s == 3)
                {
                    mz = dist.Masses.Select(v => v.ToMz(2)).ToArray();
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);
            PeakIndexingEngine indexEngine = PeakIndexingEngine.InitializeIndexingEngine(file1);

            var xic = indexEngine.GetXic(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 7, new PpmTolerance(20), 1);
            //var shiftedXic = indexEngine.GetXic((dist.Masses.First() + 0.0001).ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);

            Assert.That(xic.Count, Is.EqualTo(6));
        }

        [Test]
        public static void TestXicIsBuiltFromDifferentBins()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                if (s % 2 == 0)
                {
                    mz = dist.Masses.Select(v => (v + 0.01).ToMz(1)).ToArray();
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            PeakIndexingEngine indexEngine = PeakIndexingEngine.InitializeIndexingEngine(file1);

            // due to the wide tolerance and the way the data was constructed, the XIC and shiftedXIC should match
            var xic = indexEngine.GetXic(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
            var shiftedXic = indexEngine.GetXic((dist.Masses.First() + 0.01).ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);

            Assert.That(xic.Count, Is.EqualTo(10));
            Assert.That(shiftedXic.Count, Is.EqualTo(10));
            Assert.That(Math.Abs(xic.First().M - shiftedXic.First().M), Is.LessThan(0.00001));

            // Check that the xics are the same
            xic.AddRange(shiftedXic);
            Assert.That(xic.ToHashSet().Count(), Is.EqualTo(10)); // The intersection of the two should yield ten unique peaks
        }

        [Test]
        public static void TestInitializeIndexingEngine()
        {
            MsDataScan[] scans = new MsDataScan[10];
            for (int s = 0; s < scans.Length; s++)
            {
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 0, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
            Assert.IsNotNull(indexingEngine);
        }

        [Test]
        public static void TestGetIndexedPeakWithoutChargeState()
        {
            MsDataFile testFile = MsDataFileReader.GetDataFile(_testMzMlFullFilePath);
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testFile);
            var peak = indexingEngine.GetIndexedPeak(500, 5, new PpmTolerance(20));
            Assert.IsNotNull(peak);
            Assert.That(peak.M, Is.EqualTo(500).Within(0.001));
            Assert.That(peak.ZeroBasedScanIndex, Is.EqualTo(5));
        }

        [Test]
        public static void TestGetXicWithRetentionTime()
        {
            MsDataFile testFile = MsDataFileReader.GetDataFile(_testMzMlFullFilePath);
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testFile);
            var xic = indexingEngine.GetXic(500.0, 1.0, new PpmTolerance(20), 1);
            Assert.IsNotNull(xic);
            Assert.That(xic.Count, Is.EqualTo(9));
        }

        [Test]
        public static void TestGetXicWithStartIndex()
        {
            MsDataFile testFile = MsDataFileReader.GetDataFile(_testMzMlFullFilePath);
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testFile);
            var xic = indexingEngine.GetXic(500.0, 0, new PpmTolerance(20), 1);
            Assert.IsNotNull(xic);
            Assert.That(xic.Count, Is.EqualTo(9));
        }

        [Test]
        public static void TestMissingXic()
        {
            MsDataFile testFile = MsDataFileReader.GetDataFile(_testMzMlFullFilePath);
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testFile);
            var xic = indexingEngine.GetXic(400.0, 0, new PpmTolerance(20), 1);
            Assert.IsNotNull(xic);
            Assert.IsEmpty(xic);
        }

        [Test]
        public static void TestNotIndexedException()
         {
            PeakIndexingEngine indexingEngine = new();
            Assert.Throws<MzLibException>(() => indexingEngine.GetXic(500.0, 0, new PpmTolerance(20), 1));
        }

        [Test]
        public static void TestMassIndexingEngine()
        {
            string peptide = "PEPTIDE";
            double intensity = 1e6;
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);

            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.Select(v => v.ToMz(2)).Concat(dist.Masses.Select(v => v.ToMz(1))).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).Concat(dist.Intensities.Select(v => v * intensity * intensityMultipliers[s])).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Test the mass indexing function
            var massIndexingEngine = new MassIndexingEngine();
            Assert.IsTrue(massIndexingEngine.IndexPeaks(scans, deconParameters));

            //Test GetXIC with indexed masses
            var xic1 = massIndexingEngine.GetXic(cf.MonoisotopicMass, 5, new PpmTolerance(20), 2, 1);
            Assert.That(xic1.Count, Is.EqualTo(10));

            //Test GetXIC with different charge states
            var xic2 = massIndexingEngine.GetXic(cf.MonoisotopicMass, 5, new PpmTolerance(20), 2, 2);
            Assert.That(xic2.Count, Is.EqualTo(10));

            //Test GetXIC with different starting scan and they should return the same list of peaks
            var xic3 = massIndexingEngine.GetXic(cf.MonoisotopicMass, 1, new PpmTolerance(20), 2, 1);
            for (int i = 0; i < xic1.Count; i++)
            {
                Assert.That(Object.ReferenceEquals(xic1[i], xic3[i]));
            }

            //Get XIC with a mass that does not belong to any bins, should return an empty list
            var xic4 = massIndexingEngine.GetXic(5000.0, 5, new PpmTolerance(20), 2, 1);
            Assert.That(xic4.IsNullOrEmpty());
        }

        [Test]
        public static void TestMassIndexingExceptions()
        {
            var massIndexingEngine = new MassIndexingEngine();
            Assert.That(massIndexingEngine.IndexPeaks(new MsDataScan[] { }, new ClassicDeconvolutionParameters(1, 20, 4, 3)), Is.EqualTo(false));
            try
            {
                massIndexingEngine.GetXic(500.0,  5, new PpmTolerance(20), 2, 1);
            } catch (MzLibException e)
            {
                Assert.That(e.Message, Is.EqualTo("Error: Attempt to retrieve XIC before peak indexing was performed"));
            }
        }
    }
}
