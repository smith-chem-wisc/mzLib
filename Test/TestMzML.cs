using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml.Serialization;

namespace Test
{
    [TestFixture]
    public sealed class TestMzML
    {
        [Test]
        public static void AnotherMzMLtest()
        {
            MsDataScan[] scans = new MsDataScan[4];

            double[] intensities1 = new double[] { 1 };
            double[] mz1 = new double[] { 50 };
            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            scans[0] = new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1");

            //ms2
            double[] intensities2 = new double[] { 1 };
            double[] mz2 = new double[] { 30 };
            MzSpectrum massSpec2 = new MzSpectrum(mz2, intensities2, false);
            scans[1] = new MsDataScan(massSpec2, 2, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec2.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            double[] intensities3 = new double[] { 1 };
            double[] mz3 = new double[] { 50 };
            MzSpectrum massSpec3 = new MzSpectrum(mz3, intensities3, false);
            scans[2] = new MsDataScan(massSpec3, 3, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "3");

            //ms2
            double[] intensities4 = new double[] { 1 };
            double[] mz4 = new double[] { 30 };
            MzSpectrum massSpec4 = new MzSpectrum(mz4, intensities4, false);
            scans[3] = new MsDataScan(massSpec4, 4, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec2.SumOfAllY, null, null, "4", 50, null, null, 50, 1, DissociationType.CID, 3, null);

            FakeMsDataFile f = new FakeMsDataFile(scans);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"));

            var scanWithPrecursor = ok.GetAllScansList().Last(b => b.MsnOrder != 1);

            Assert.AreEqual(3, scanWithPrecursor.OneBasedPrecursorScanNumber);
        }

        [Test]
        public void LoadBadMzml()
        {
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "asdfasdfasdfasdfasdf.mzML")); // just to be sure
            Assert.Throws<FileNotFoundException>(() => Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "asdfasdfasdfasdfasdf.mzML")));
        }

        [Test]
        public static void TestPeakTrimmingWithOneWindow()
        {
            Random rand = new Random();
            int numPeaks = 200;
            double minRatio = 0.01;
            int numWindows = 1;

            var testFilteringParams = new FilteringParams(numPeaks, minRatio, numWindows, true, true);
            List<(double mz, double intensity)> myPeaks = new List<(double mz, double intensity)>();

            for (int mz = 400; mz < 1600; mz++)
            {
                myPeaks.Add((mz, rand.Next(1000, 1000000)));
            }

            double myMaxIntensity = myPeaks.Max(p => p.intensity);
            var myPeaksOrderedByIntensity = myPeaks.OrderByDescending(p => p.intensity).ToList();
            myPeaksOrderedByIntensity = myPeaksOrderedByIntensity.Take(numPeaks).ToList();
            myPeaksOrderedByIntensity = myPeaksOrderedByIntensity.Where(p => (p.intensity / myMaxIntensity) > minRatio).ToList();
            double sumOfAllIntensities = myPeaksOrderedByIntensity.Sum(p => p.intensity);

            double[] intensities1 = myPeaks.Select(p => p.intensity).ToArray();
            double[] mz1 = myPeaks.Select(p => p.mz).ToArray();

            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), testFilteringParams);

            int expNumPeaks = ok.GetAllScansList().First().MassSpectrum.XArray.Length;
            double expMinRatio = ok.GetAllScansList().First().MassSpectrum.YArray.Min(p => p / ok.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY).Value;
            List<(double mz, double intensity)> myExpPeaks = new List<(double mz, double intensity)>();

            for (int i = 0; i < ok.GetAllScansList().First().MassSpectrum.YArray.Length; i++)
            {
                myExpPeaks.Add((ok.GetAllScansList().First().MassSpectrum.XArray[i], ok.GetAllScansList().First().MassSpectrum.YArray[i]));
            }

            Assert.That(Math.Round(myMaxIntensity, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            Assert.That(Math.Round(sumOfAllIntensities, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.SumOfAllY, 0));
            Assert.That(myPeaksOrderedByIntensity.Count == ok.GetAllScansList().First().MassSpectrum.XArray.Length);
            Assert.That(expMinRatio >= minRatio);
            Assert.That(!myExpPeaks.Except(myPeaksOrderedByIntensity).Any());
            Assert.That(!myPeaksOrderedByIntensity.Except(myExpPeaks).Any());
        }

        [Test]
        public static void TestPeakTrimmingWithThreeWindows()
        {
            Random rand = new Random();
            int numPeaksPerWindow = 200;
            double minRatio = 0.01;
            int numWindows = 3;

            var testFilteringParams = new FilteringParams(numPeaksPerWindow, minRatio, numWindows, true, true);
            List<(double mz, double intensity)> myPeaks = new List<(double mz, double intensity)>();
            List<(double mz, double intensity)> myPeaksWindow1 = new List<(double mz, double intensity)>();
            List<(double mz, double intensity)> myPeaksWindow2 = new List<(double mz, double intensity)>();
            List<(double mz, double intensity)> myPeaksWindow3 = new List<(double mz, double intensity)>();

            int peakCounter = 0;
            for (int mz = 400; mz < 1599; mz++)
            {
                peakCounter++;

                int intensity = rand.Next(1000, 10000);
                myPeaks.Add((mz, intensity));
                if (peakCounter <= 399)
                {
                    myPeaksWindow1.Add((mz, intensity));
                }
                else if (peakCounter <= 799)
                {
                    myPeaksWindow2.Add((mz, intensity));
                }
                else if (peakCounter <= 1200)
                {
                    myPeaksWindow3.Add((mz, intensity));
                }
            }

            double[] intensities1 = myPeaks.Select(p => p.intensity).ToArray();
            double[] mz1 = myPeaks.Select(p => p.mz).ToArray();

            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), testFilteringParams);

            int expNumPeaks = ok.GetAllScansList().First().MassSpectrum.XArray.Length;
            double expMinRatio = ok.GetAllScansList().First().MassSpectrum.YArray.Min(p => p / ok.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY).Value;
            List<(double mz, double intensity)> myExpPeaks = new List<(double mz, double intensity)>();

            for (int i = 0; i < ok.GetAllScansList().First().MassSpectrum.YArray.Length; i++)
            {
                myExpPeaks.Add((ok.GetAllScansList().First().MassSpectrum.XArray[i], ok.GetAllScansList().First().MassSpectrum.YArray[i]));
            }

            double myMaxIntensity = myPeaks.Max(p => p.intensity);

            myPeaksWindow1 = myPeaksWindow1.OrderByDescending(p => p.intensity).Take(numPeaksPerWindow).Where(p => (p.intensity / myMaxIntensity) > minRatio).ToList();

            myPeaksWindow2 = myPeaksWindow2.OrderByDescending(p => p.intensity).Take(numPeaksPerWindow).Where(p => (p.intensity / myMaxIntensity) > minRatio).ToList();

            myPeaksWindow3 = myPeaksWindow3.OrderByDescending(p => p.intensity).Take(numPeaksPerWindow).Where(p => (p.intensity / myMaxIntensity) > minRatio).ToList();

            var allWindowPeaksCombined = myPeaksWindow1.Concat(myPeaksWindow2).Concat(myPeaksWindow3).ToList();

            double sumOfAllIntensities = allWindowPeaksCombined.Sum(p => p.intensity);

            Assert.That(Math.Round(myMaxIntensity, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            Assert.That(allWindowPeaksCombined.Count == ok.GetAllScansList().First().MassSpectrum.XArray.Length);
            Assert.That(expMinRatio >= minRatio);
            Assert.That(Math.Round(sumOfAllIntensities, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.SumOfAllY, 0));
            Assert.That(!myExpPeaks.Except(allWindowPeaksCombined).Any());
            Assert.That(!allWindowPeaksCombined.Except(myExpPeaks).Any());
        }

        [Test]
        public static void TestPeakTrimmingWithTooManyWindows()
        {
            Random rand = new Random();
            int numPeaks = 200;
            double minRatio = 0.01;
            int numWindows = 10;

            var testFilteringParams = new FilteringParams(numPeaks, minRatio, numWindows, true, true);
            // only 1 peak but 10 windows
            List<(double mz, double intensity)> myPeaks = new List<(double mz, double intensity)>
            {
                {(400, rand.Next(1000, 1000000)) }
            };

            double[] intensities1 = myPeaks.Select(p => p.intensity).ToArray();
            double[] mz1 = myPeaks.Select(p => p.mz).ToArray();

            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), testFilteringParams);

            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.SumOfAllY, 0));
            Assert.That(1 == ok.GetAllScansList().First().MassSpectrum.XArray.Length);
            Assert.That(Math.Round(myPeaks[0].mz, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.XArray[0], 0));
            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(ok.GetAllScansList().First().MassSpectrum.YArray[0], 0));
        }

        [Test]
        public static void WriteEmptyScan()
        {
            double[] intensities1 = new double[] { };
            double[] mz1 = new double[] { };
            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzmlWithEmptyScan.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzmlWithEmptyScan.mzML"));

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzmlWithEmptyScan2.mzML"), false);

            var testFilteringParams = new FilteringParams(200, 0.01, 5, true, true);
            ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzmlWithEmptyScan2.mzML"), testFilteringParams);
        }

        [Test]
        public static void DifferentAnalyzersTest()
        {
            MsDataScan[] scans = new MsDataScan[2];

            double[] intensities1 = new double[] { 1 };
            double[] mz1 = new double[] { 50 };
            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            scans[0] = new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1");

            double[] intensities2 = new double[] { 1 };
            double[] mz2 = new double[] { 30 };
            MzSpectrum massSpec2 = new MzSpectrum(mz2, intensities2, false);
            scans[1] = new MsDataScan(massSpec2, 2, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.IonTrap3D, massSpec2.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            FakeMsDataFile f = new FakeMsDataFile(scans);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "asdfefsf.mzML"), false);

            Mzml ok = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "asdfefsf.mzML"));

            Assert.AreEqual(MZAnalyzerType.Orbitrap, ok.GetAllScansList().First().MzAnalyzer);
            Assert.AreEqual(MZAnalyzerType.IonTrap3D, ok.GetAllScansList().Last().MzAnalyzer);
        }

        [Test]
        public static void Mzid111Test()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111
            {
                DataCollection = new mzIdentML111.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML111.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML111.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML111.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML111.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML111.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML111.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML111.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML111.Generated.CVParamType[]
                {
                    new mzIdentML111.Generated.CVParamType
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] =
                new mzIdentML111.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation =
                new mzIdentML111.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] =
                new mzIdentML111.Generated.IonTypeType
                {
                    FragmentArray = new mzIdentML111.Generated.FragmentArrayType[1]
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] =
                new mzIdentML111.Generated.FragmentArrayType
                {
                    values = new float[3] { 200, 300, 400 }
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef =
                new mzIdentML111.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] =
                new mzIdentML111.Generated.PeptideEvidenceRefType
                {
                    peptideEvidence_ref = "PE_1"
                };
            _mzid.DataCollection.Inputs = new mzIdentML111.Generated.InputsType
            {
                SpectraData = new mzIdentML111.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML111.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML111.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML111.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML111.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML111.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML111.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML111.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML111.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML111.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML111.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML111.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML111.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML111.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML111.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML111.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML111.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML111.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public static void Mzid120Test()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120()
            {
                DataCollection = new mzIdentML120.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML120.Generated.AnalysisDataType()
            {
                SpectrumIdentificationList = new mzIdentML120.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML120.Generated.SpectrumIdentificationListType()
            {
                SpectrumIdentificationResult = new mzIdentML120.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML120.Generated.SpectrumIdentificationResultType()
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML120.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML120.Generated.SpectrumIdentificationItemType()
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
                {
                    new mzIdentML120.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML120.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML120.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML120.Generated.IonTypeType()
            {
                FragmentArray = new mzIdentML120.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML120.Generated.FragmentArrayType()
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML120.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML120.Generated.PeptideEvidenceRefType()
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML120.Generated.InputsType()
            {
                SpectraData = new mzIdentML120.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML120.Generated.SpectraDataType()
            {
                FileFormat = new mzIdentML120.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML120.Generated.CVParamType()
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML120.Generated.SequenceCollectionType()
            {
                PeptideEvidence = new mzIdentML120.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML120.Generated.PeptideEvidenceType()
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML120.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML120.Generated.PeptideType()
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML120.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML120.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML120.Generated.DBSequenceType()
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML120.Generated.ModificationType()
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML120.Generated.CVParamType()
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML120.Generated.AnalysisProtocolCollectionType()
            {
                SpectrumIdentificationProtocol = new mzIdentML120.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML120.Generated.SpectrumIdentificationProtocolType()
            {
                ParentTolerance = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML120.Generated.CVParamType()
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML120.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML120.Generated.CVParamType()
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements(@"elements.dat");
        }

        [Test]
        public void LoadMzmlTest()
        {
            Assert.Throws<AggregateException>(() =>
            {
                Mzml.LoadAllStaticData(@"tiny.pwiz.1.1.mzML");
            }, "Reading profile mode mzmls not supported");
        }

        [Test]
        public void LoadMzmlFromConvertedMGFTest()
        {
            Mzml a = Mzml.LoadAllStaticData(@"tester.mzML");

            var ya = a.GetOneBasedScan(1).MassSpectrum;
            Assert.AreEqual(192, ya.Size);
            var ya2 = a.GetOneBasedScan(3).MassSpectrum;
            Assert.AreEqual(165, ya2.Size);
            var ya3 = a.GetOneBasedScan(5).MassSpectrum;
            Assert.AreEqual(551, ya3.Size);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, "CreateFileFromConvertedMGF.mzML", false);

            Mzml b = Mzml.LoadAllStaticData(@"CreateFileFromConvertedMGF.mzML");

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(b, "CreateFileFromConvertedMGF2.mzML", false);
        }

        [Test]
        public void WriteMzmlTest()
        {
            var peptide = new Peptide("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR");
            OldSchoolChemicalFormulaModification carbamidomethylationOfCMod = new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula("H3C2NO"), "carbamidomethylation of C", ModificationSites.C);
            peptide.AddModification(carbamidomethylationOfCMod);

            MzSpectrum MS1 = CreateSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);

            MzSpectrum MS2 = CreateMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

            MsDataScan[] Scans = new MsDataScan[2];

            Scans[0] = new MsDataScan(MS1, 1, 1, true, Polarity.Positive, 1.0, new MzRange(300, 2000), " first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, 1, null, "scan=1");
            //            scans[1] = new MsDataScanZR(massSpec2, 2, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.IonTrap3D, massSpec2.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            Scans[1] = new MsDataScan(MS2, 2, 2, true, Polarity.Positive, 2.0, new MzRange(100, 1500), " second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, 1, null, "scan=2", 1134.26091302033, 3, 0.141146966879759, 1134.3, 1, DissociationType.Unknown, 1, 1134.26091302033);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var oldFirstValue = myMsDataFile.GetOneBasedScan(1).MassSpectrum.FirstX;

            var secondScan = myMsDataFile.GetOneBasedScan(2);
            Assert.AreEqual(1, secondScan.IsolationRange.Maximum - secondScan.IsolationRange.Minimum);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, "argh.mzML", false);

            Mzml okay = Mzml.LoadAllStaticData(@"argh.mzML");
            okay.GetOneBasedScan(2);

            Assert.AreEqual(1, okay.GetClosestOneBasedSpectrumNumber(1));
            Assert.AreEqual(2, okay.GetClosestOneBasedSpectrumNumber(2));

            var newFirstValue = okay.GetOneBasedScan(1).MassSpectrum.FirstX;
            Assert.AreEqual(oldFirstValue.Value, newFirstValue.Value, 1e-9);

            var secondScan2 = okay.GetOneBasedScan(2);

            Assert.AreEqual(1, secondScan2.IsolationRange.Maximum - secondScan2.IsolationRange.Minimum);

            secondScan2.MassSpectrum.ReplaceXbyApplyingFunction((a) => 44);
            Assert.AreEqual(44, secondScan2.MassSpectrum.LastX);
        }

        [Test]
        public void MzidTest()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110
            {
                DataCollection = new mzIdentML110.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML110.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML110.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML110.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML110.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML110.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML110.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML110.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML110.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML110.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML110.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].passThreshold = true;

            _mzid.DataCollection.Inputs = new mzIdentML110.Generated.InputsType
            {
                SpectraData = new mzIdentML110.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML110.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML110.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML110.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML110.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                start = 2,
                end = 34,
                isDecoy = false,
                peptide_ref = "P_1",
                dBSequence_ref = "DB_1",
                id = "PE_1"
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML110.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML110.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML110.Generated.ModificationType[1],
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML110.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML110.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML110.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML110.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML110.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML110.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid110Test()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110
            {
                DataCollection = new mzIdentML110.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML110.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML110.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML110.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML110.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML110.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML110.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML110.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML110.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML110.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML110.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML110.Generated.InputsType
            {
                SpectraData = new mzIdentML110.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML110.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML110.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML110.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML110.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML110.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML110.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML110.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML110.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML110.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML110.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML110.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML110.Generated.SpectrumIdentificationProtocolType()
            {
                ParentTolerance = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML110.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid111Test_()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111
            {
                DataCollection = new mzIdentML111.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML111.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML111.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML111.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML111.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML111.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML111.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML111.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML111.Generated.CVParamType[]
                {
                    new mzIdentML111.Generated.CVParamType
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] =
                new mzIdentML111.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation =
                new mzIdentML111.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] =
                new mzIdentML111.Generated.IonTypeType
                {
                    FragmentArray = new mzIdentML111.Generated.FragmentArrayType[1]
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] =
                new mzIdentML111.Generated.FragmentArrayType
                {
                    values = new float[3] { 200, 300, 400 }
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef =
                new mzIdentML111.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] =
                new mzIdentML111.Generated.PeptideEvidenceRefType
                {
                    peptideEvidence_ref = "PE_1"
                };
            _mzid.DataCollection.Inputs = new mzIdentML111.Generated.InputsType
            {
                SpectraData = new mzIdentML111.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML111.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML111.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML111.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML111.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML111.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML111.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML111.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML111.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML111.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML111.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML111.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML111.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML111.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML111.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML111.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML111.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML111.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid120Test_()
        {
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120
            {
                DataCollection = new mzIdentML120.Generated.DataCollectionType()
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML120.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML120.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML120.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML120.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML120.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML120.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML120.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
                {
                    new mzIdentML120.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML120.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML120.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML120.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML120.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML120.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML120.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML120.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML120.Generated.InputsType
            {
                SpectraData = new mzIdentML120.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML120.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML120.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML120.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML120.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML120.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML120.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML120.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML120.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML120.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML120.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML120.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML120.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML120.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML120.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML120.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML120.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML120.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML120.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML120.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };
            TextWriter writer = new StreamWriter("myIdentifications.mzid");
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void MzmlFindPrecursorReferenceScan()
        {
            //some ms2 scans dont have properly assigned precursor scans
            //this unit test is intended to test the fallback option in MzML\Mzml.cs @ lines 459-47
            //constructs three ms1 scans and three ms2 scans
            //if ms2 scan precusor reference is null, assumes most recent ms1 scan is correct reference

            MsDataScan[] scans = new MsDataScan[6];
            MsDataScan[] scans1 = new MsDataScan[4];

            double[] intensities0 = new double[] { 1 };
            double[] mz0 = new double[] { 50 };
            MzSpectrum massSpec0 = new MzSpectrum(mz0, intensities0, false);
            scans[0] = new MsDataScan(massSpec0, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec0.SumOfAllY, null, null, "1");
            scans1[0] = scans[0];

            double[] intensities1 = new double[] { 1 };
            double[] mz1 = new double[] { 50 };
            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            scans[1] = new MsDataScan(massSpec0, 2, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1");

            double[] intensities2 = new double[] { 1 };
            double[] mz2 = new double[] { 50 };
            MzSpectrum massSpec2 = new MzSpectrum(mz2, intensities2, false);
            scans[2] = new MsDataScan(massSpec2, 3, 1, true, Polarity.Positive, 1, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec2.SumOfAllY, null, null, "1");

            //ms2
            double[] intensities3 = new double[] { 1 };
            double[] mz3 = new double[] { 30 };
            MzSpectrum massSpec3 = new MzSpectrum(mz3, intensities3, false);
            scans[3] = new MsDataScan(massSpec3, 4, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec3.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 3, null);
            scans1[1] = new MsDataScan(massSpec3, 2, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec3.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            //ms2
            double[] intensities4 = new double[] { 1 };
            double[] mz4 = new double[] { 30 };
            MzSpectrum massSpec4 = new MzSpectrum(mz4, intensities4, false);
            scans[4] = new MsDataScan(massSpec4, 5, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec4.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 3, null);
            scans1[2] = new MsDataScan(massSpec4, 3, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec4.SumOfAllY, null, null, "2", 50, null, null, 50, 1, DissociationType.CID, 1, null);

            //ms2
            double[] intensities5 = new double[] { 1 };
            double[] mz5 = new double[] { 30 };
            MzSpectrum massSpec5 = new MzSpectrum(mz5, intensities5, false);
            scans[5] = new MsDataScan(massSpec5, 6, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec5.SumOfAllY, null, null, "4", 50, null, null, 50, 1, DissociationType.CID, null, null);
            scans1[3] = new MsDataScan(massSpec5, 4, 2, true, Polarity.Positive, 2, new MzRange(1, 100), "f", MZAnalyzerType.Orbitrap, massSpec5.SumOfAllY, null, null, "4", 50, null, null, 50, 1, DissociationType.CID, null, null);

            FakeMsDataFile fakeFile = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fakeFile, Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"), false);
            Mzml fakeMzml = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"));

            FakeMsDataFile fakeFile1 = new FakeMsDataFile(scans1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fakeFile1, Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"), false);
            Mzml fakeMzml1 = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"));

            Assert.AreEqual(3, fakeMzml.GetAllScansList().ElementAt(5).OneBasedPrecursorScanNumber);
            Assert.AreEqual(1, fakeMzml1.GetAllScansList().ElementAt(3).OneBasedPrecursorScanNumber);
        }

        private MzSpectrum CreateMS2spectrum(IEnumerable<Fragment> fragments, int v1, int v2)
        {
            List<double> allMasses = new List<double>();
            List<double> allIntensities = new List<double>();
            foreach (ChemicalFormulaFragment f in fragments)
            {
                var spec = CreateSpectrum(f.ThisChemicalFormula, v1, v2, 2);
                for (int i = 0; i < spec.Size; i++)
                {
                    allMasses.Add(spec.XArray[i]);
                    allIntensities.Add(spec.YArray[i]);
                }
            }
            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensities.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);
            return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        private MzSpectrum CreateSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
        {
            IsotopicDistribution isodist = IsotopicDistribution.GetDistribution(f, 0.1);

            return new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);
            //massSpectrum1 = massSpectrum1.FilterByNumberOfMostIntense(5);

            //var chargeToLookAt = minCharge;
            //var correctedSpectrum = massSpectrum1.NewSpectrumApplyFunctionToX(s => s.ToMz(chargeToLookAt));

            //List<double> allMasses = new List<double>();
            //List<double> allIntensitiess = new List<double>();

            //while (correctedSpectrum.FirstX > lowerBound)
            //{
            //    foreach (var thisPeak in correctedSpectrum)
            //    {
            //        if (thisPeak.Mz > lowerBound && thisPeak.Mz < upperBound)
            //        {
            //            allMasses.Add(thisPeak.Mz);
            //            allIntensitiess.Add(thisPeak.Intensity);
            //        }
            //    }
            //    chargeToLookAt += 1;
            //    correctedSpectrum = massSpectrum1.NewSpectrumApplyFunctionToX(s => s.ToMz(chargeToLookAt));
            //}

            //var allMassesArray = allMasses.ToArray();
            //var allIntensitiessArray = allIntensitiess.ToArray();

            //Array.Sort(allMassesArray, allIntensitiessArray);

            //return new MzmlMzSpectrum(allMassesArray, allIntensitiessArray, false);
        }
    }
}