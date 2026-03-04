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

namespace Test.FlashLFQ
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

            MsDataScan[] scanArrayAvailableToAllUnitTests = new MsDataScan[9];
            double[] intensityMultipliers = { 1, 2, 3, 4, 5, 4, 3, 2, 1 };

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < scanArrayAvailableToAllUnitTests.Length; s++)
            {
                double[] mz = new double[] { 500, 500.5, 501, 501.5, 502 };
                double[] intensities = Enumerable.Repeat(intensityMultipliers[s]*intensity, 5).ToArray();

                // add the scan
                scanArrayAvailableToAllUnitTests[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scanArrayAvailableToAllUnitTests),
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

            var xic = indexEngine.GetXicByScanIndex(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
            var shiftedXic = indexEngine.GetXicByScanIndex((dist.Masses.First() + 0.0001).ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);

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

            var xic = indexEngine.GetXicByScanIndex(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 7, new PpmTolerance(20), 1);
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
            var xic = indexEngine.GetXicByScanIndex(dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
            var shiftedXic = indexEngine.GetXicByScanIndex((dist.Masses.First() + 0.01).ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);

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
            var xic = indexingEngine.GetXicByScanIndex(500.0, 0, new PpmTolerance(20), 1);
            Assert.IsNotNull(xic);
            Assert.That(xic.Count, Is.EqualTo(9));
        }

        [Test]
        public static void TestMissingXic()
        {
            MsDataFile testFile = MsDataFileReader.GetDataFile(_testMzMlFullFilePath);
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testFile);
            var xic = indexingEngine.GetXicByScanIndex(400.0, 0, new PpmTolerance(20), 1);
            Assert.IsNotNull(xic);
            Assert.IsEmpty(xic);
        }

        [Test]
        public static void TestNotIndexedException()
         {
            PeakIndexingEngine indexingEngine = new();
            Assert.Throws<MzLibException>(() => indexingEngine.GetXicByScanIndex(500.0, 0, new PpmTolerance(20), 1));
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
            int minCharge = 1; // Set a minimum charge state for indexing

            // Create mzSpectra
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.Select(v => v.ToMz(minCharge + 1)).Concat(dist.Masses.Select(v => v.ToMz(minCharge))).ToArray();
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
            var xic1 = massIndexingEngine.GetXicByScanIndex(cf.MonoisotopicMass, 5, new PpmTolerance(20), 2, 1);
            Assert.That(xic1.Count, Is.EqualTo(10));

            //Test GetXIC with different charge states
            var xic2 = massIndexingEngine.GetXicByScanIndex(cf.MonoisotopicMass, 5, new PpmTolerance(20), 2, 2);
            Assert.That(xic2.Count, Is.EqualTo(10));

            //Test GetXIC with different starting scan and they should return the same list of peaks
            var xic3 = massIndexingEngine.GetXicByScanIndex(cf.MonoisotopicMass, 1, new PpmTolerance(20), 2, 1);
            for (int i = 0; i < xic1.Count; i++)
            {
                Assert.That(Object.ReferenceEquals(xic1[i], xic3[i]));
            }

            //Get XIC with a mass that does not belong to any bins, should return an empty list
            var xic4 = massIndexingEngine.GetXicByScanIndex(5000.0, 5, new PpmTolerance(20), 2, 1);
            Assert.That(xic4.IsNullOrEmpty());

            //Test the IndexPeaks method with a scan where the monoisotopic mass is less than the minimum
            massIndexingEngine = new MassIndexingEngine();
            double minMass = cf.MonoisotopicMass + 10; // Set a minimum mass that is greater than the peptide's monoisotopic mass

            //no peaks indexed
            Assert.Throws<MzLibException>(() => massIndexingEngine.IndexPeaks(scans, deconParameters, null, minMass, minCharge));


            //Test the IndexPeaks method with a scan where the scan charges are less than the minimum charge
            massIndexingEngine = new MassIndexingEngine();

            //no peaks indexed
            Assert.Throws<MzLibException>(() => massIndexingEngine.IndexPeaks(scans, deconParameters, null, cf.MonoisotopicMass, 10));
        }
        [Test]
        public static void TestIndexPeaksWithAPeptideThatIsTooLarge()
        {
            string peptide = "PEPTIDE";
            double intensity = 1e6;
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);

            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            List<ChemicalFormula> chemicalFormulaList = new List<ChemicalFormula>();
            for (int times = 0; times < 40; times++)
            {
                chemicalFormulaList.Add(cf);
            }
            ChemicalFormula cfTooLarge = ChemicalFormula.Combine(chemicalFormulaList);
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cfTooLarge, 0.125, 1e-8);
            int minCharge = 1; // Set a minimum charge state for indexing

            // Create mzSpectra
            for (int s = 0; s < scans.Length; s++)
            {
                double[] mz = dist.Masses.Select(v => v.ToMz(minCharge + 1)).Concat(dist.Masses.Select(v => v.ToMz(minCharge))).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).Concat(dist.Intensities.Select(v => v * intensity * intensityMultipliers[s])).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
            //Test the mass indexing function
            var massIndexingEngine = new MassIndexingEngine();
            // Expected behavior: IndexPeaks should throw MzLibException when no peaks can be indexed due to all scans having monoisotopic mass larger than the maximum mass.
            Assert.Throws<MzLibException>(() => massIndexingEngine.IndexPeaks(scans, deconParameters, null, 0, minCharge));
        }
        [Test]
        public static void TestMassIndexingEngineWithNearlyIsobaricPeptides()
        {
            int totalScans = 10;
            MsDataScan[] scans = new MsDataScan[totalScans];
            string lowerMassPeptide = "PEPTIDEFP";
            string higherMassPeptide = "PEPTIDELM"; // A peptide with a mass that is 0.01 Da higher than the lower mass peptide
            double intensity = 1e6;
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);

            double[] lowerMassIntensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            double[] higherMassIntensityMultipliers = { 1, 10, 20, 10, 1};
            ChemicalFormula chemicalFormulaLowerMassPeptide = new Proteomics.AminoAcidPolymer.Peptide(lowerMassPeptide).GetChemicalFormula();
            ChemicalFormula chemicalFormulaHigherMassPeptide = new Proteomics.AminoAcidPolymer.Peptide(higherMassPeptide).GetChemicalFormula();
            IsotopicDistribution lowerMassIsotopicDistribution = IsotopicDistribution.GetDistribution(chemicalFormulaLowerMassPeptide, 0.125, 1e-8);
            IsotopicDistribution higherMassIsotopicDistribution = IsotopicDistribution.GetDistribution(chemicalFormulaHigherMassPeptide, 0.125, 1e-8);
            Assert.That(1043.48114179643, Is.EqualTo(chemicalFormulaLowerMassPeptide.MonoisotopicMass).Within(0.001));
            Assert.That(1043.48451309975, Is.EqualTo(chemicalFormulaHigherMassPeptide.MonoisotopicMass).Within(0.001));

            double ppmDifference = (chemicalFormulaHigherMassPeptide.MonoisotopicMass - chemicalFormulaLowerMassPeptide.MonoisotopicMass) / chemicalFormulaLowerMassPeptide.MonoisotopicMass * 1e6;
            Assert.That(ppmDifference, Is.EqualTo(3.23).Within(0.01)); // 0.01 Da difference at 1043.48114179643 Da corresponds to approximately 3.23 ppm

            // Create mzSpectra
            for (int s = 0; s < totalScans; s++)
            {
                List<double> mzValues = lowerMassIsotopicDistribution.Masses.Select(v => v.ToMz(2)).Concat(lowerMassIsotopicDistribution.Masses.Select(v => v.ToMz(1))).ToList();
                List<double> intensities = lowerMassIsotopicDistribution.Intensities.Select(v => v * intensity * lowerMassIntensityMultipliers[s]).Concat(lowerMassIsotopicDistribution.Intensities.Select(v => v * intensity * lowerMassIntensityMultipliers[s])).ToList();

                if(s < 5)
                {
                    // For the first half of the scans, add the higher mass peptide
                    mzValues.AddRange(higherMassIsotopicDistribution.Masses.Select(v => v.ToMz(2)).Concat(higherMassIsotopicDistribution.Masses.Select(v => v.ToMz(1))));
                    intensities.AddRange(higherMassIsotopicDistribution.Intensities.Select(v => v * intensity * higherMassIntensityMultipliers[s]).Concat(higherMassIsotopicDistribution.Intensities.Select(v => v * intensity * higherMassIntensityMultipliers[s])));
                }


                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mzValues.ToArray(), intensities.ToArray(), false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(100, 2000), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Test the mass indexing function
            var massIndexingEngine = new MassIndexingEngine();
            Assert.IsTrue(massIndexingEngine.IndexPeaks(scans, deconParameters));

            //Test GetXIC with indexed masses
            var lowerMassXic1 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaLowerMassPeptide.MonoisotopicMass, 5, new PpmTolerance(20), 2, 1);
            Assert.That(lowerMassXic1.Count, Is.EqualTo(10));

            //Test GetXIC with different charge states
            var lowerMassXic2 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaLowerMassPeptide.MonoisotopicMass, 5, new PpmTolerance(20), 2, 2);
            Assert.That(lowerMassXic2.Count, Is.EqualTo(10));

            //Test GetXIC with different starting scan and they should return the same list of peaks
            var lowerMassXic3 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaLowerMassPeptide.MonoisotopicMass, 1, new PpmTolerance(20), 2, 1);
            for (int i = 0; i < lowerMassXic1.Count; i++)
            {
                Assert.That(lowerMassXic1[i], Is.SameAs(lowerMassXic3[i]));
            }

            //Test GetXIC with indexed masses
            //even thought the higher mass peptide is only present in the first 5 scans, it should still return 10 scans for the XIC because of the wide tolerance
            var higherMassXic1 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 5, new PpmTolerance(20), 2, 1);
            Assert.That(higherMassXic1.Count, Is.EqualTo(10));

            //Test GetXIC with different charge states
            var higherMassXic2 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 5, new PpmTolerance(20), 2, 2);
            Assert.That(higherMassXic2.Count, Is.EqualTo(10));

            //Test GetXIC with different starting scan and they should return the same list of peaks
            var higherMassXic3 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 1, new PpmTolerance(20), 2, 1);
            for (int i = 0; i < higherMassXic1.Count; i++)
            {
                Assert.That(higherMassXic1[i], Is.SameAs(higherMassXic3[i]));
            }

            //lower the tolerance here to separate the peaks. Now we should get 5 scans in the xic for the higher mass peptide
            //Test GetXIC with indexed masses
            higherMassXic1 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 5, new PpmTolerance(1), 2, 1);
            //this should be 5 but there is currently an error in deconvolution. when that gets fixed this should change to 5.
            Assert.That(higherMassXic1.Count, Is.EqualTo(4));

            //Test GetXIC with different charge states
            higherMassXic2 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 5, new PpmTolerance(1), 2, 2);
            //this should be 5 but there is currently an error in deconvolution. when that gets fixed this should change to 5.
            Assert.That(higherMassXic2.Count, Is.EqualTo(4));

            //Test GetXIC with different starting scan and they should return the same list of peaks
            higherMassXic3 = massIndexingEngine.GetXicByScanIndex(chemicalFormulaHigherMassPeptide.MonoisotopicMass, 1, new PpmTolerance(1), 2, 1);
            for (int i = 0; i < higherMassXic1.Count; i++)
            {
                Assert.That(higherMassXic1[i], Is.SameAs(higherMassXic3[i]));
            }

            //Get XIC with a mass that does not belong to any bins, should return an empty list
            var emptyXic4 = massIndexingEngine.GetXicByScanIndex(5000.0, 5, new PpmTolerance(20), 2, 1);
            Assert.That(emptyXic4.IsNullOrEmpty());

            //Test for proper handling of a null scan.
            MsDataScan nullScan = null;
            MsDataScan[] scansWithNull = scans.ToList().Append(nullScan).ToArray();
            massIndexingEngine = new MassIndexingEngine();
            Assert.That(massIndexingEngine.IndexPeaks(scansWithNull, deconParameters), Is.EqualTo(true));
            Assert.That(massIndexingEngine.GetXicByScanIndex(chemicalFormulaLowerMassPeptide.MonoisotopicMass, 5, new PpmTolerance(20), 2, 1).Count, Is.EqualTo(10));
        }
        [Test]
        public static void TestMassIndexingExceptions()
        {
            var massIndexingEngine = new MassIndexingEngine();
            Assert.That(massIndexingEngine.IndexPeaks(new MsDataScan[] { }, new ClassicDeconvolutionParameters(1, 20, 4, 3)), Is.EqualTo(false));
            try
            {
                massIndexingEngine.GetXicByScanIndex(500.0,  5, new PpmTolerance(20), 2, 1);
            } catch (MzLibException e)
            {
                Assert.That(e.Message, Is.EqualTo("Error: Attempt to retrieve XIC before peak indexing was performed"));
            }
        }
        [Test]
        public static void TestMassIndexingEngineNullable()
        {
            List<MsDataScan> scans = new List<MsDataScan>();
            for (int i = 0; i < 10; i++)
            {
                MsDataScan scan = null;
                scans.Add(scan);
            }
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var massIndexingEngine = MassIndexingEngine.InitializeMassIndexingEngine(scans.ToArray(), deconParameters);
            Assert.That(massIndexingEngine, Is.Null);
        }
    }
}
