using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Xml.Serialization;
using Chemistry;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics.AminoAcidPolymer;
using Readers;
using Stopwatch = System.Diagnostics.Stopwatch;
using Easy.Common.Extensions;
using MathNet.Numerics.Statistics;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzML
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

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

            f.ExportAsMzML(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"), false);

            var reader =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"));
            reader.LoadAllStaticData();
            var scanWithPrecursor = reader.GetAllScansList().Last(b => b.MsnOrder != 1);
            Assert.AreEqual(3, scanWithPrecursor.OneBasedPrecursorScanNumber);
        }

        [Test]
        public void TextScanNumberNotReportedInMzml()
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "noScanNumber.mzML");
            FilteringParams filter = new(200, 0.01, 1, null, false, false, true);
            bool lookForPrecursorScans = false;
            Assert.Throws<AggregateException>(() => MsDataFileReader.GetDataFile(origDataFile).LoadAllStaticData(filter, 1), "Precursor scan number not define for scan=1");
        }

        [Test]
        public static void ReadMzMlInNewEra()
        {
            Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "BinGenerationTest.mzML");
            FilteringParams filter = new FilteringParams(200, 0.01, 1, null, false, false, true);
            var reader = MsDataFileReader.GetDataFile(origDataFile);
            reader.LoadAllStaticData(filter, 1);
            MyMsDataFiles[origDataFile] = reader;
            var scans = MyMsDataFiles[origDataFile].GetAllScansList();

            Assert.AreEqual(6, scans[0].MassSpectrum.XArray.Count());
            Assert.AreEqual(20, scans[1].MassSpectrum.XArray.Count());
        }

        [Test]
        public void LoadBadMzml()
        {
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "asdfasdfasdfasdfasdf.mzML")); // just to be sure

            Assert.Throws<FileNotFoundException>(() =>
            {
                var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                    "asdfasdfasdfasdfasdf.mzML"));
                reader.LoadAllStaticData();
            });
        }

        [Test]
        public static void TestPeakTrimmingWithOneWindow()
        {
            int numPeaks = 200;
            double minRatio = 0.01;

            var testFilteringParams = new FilteringParams(numPeaks, minRatio, null, 1, false, true, false);
            List<(double mz, double intensity)> myPeaks = new List<(double mz, double intensity)>();

            for (int mz = 400; mz < 1600; mz++)
            {
                myPeaks.Add((mz, 10d * mz));
            }

            double myMaxIntensity = myPeaks.Max(p => p.intensity);
            var myPeaksOrderedByIntensity = myPeaks.OrderByDescending(p => p.intensity).ToList();
            myPeaksOrderedByIntensity = myPeaksOrderedByIntensity.Take(numPeaks).ToList();
            myPeaksOrderedByIntensity = myPeaksOrderedByIntensity.Where(p => p.intensity / myMaxIntensity > minRatio).ToList();
            double sumOfAllIntensities = myPeaksOrderedByIntensity.Sum(p => p.intensity);

            double[] intensities1 = myPeaks.Select(p => p.intensity).ToArray();
            double[] mz1 = myPeaks.Select(p => p.mz).ToArray();

            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            f.ExportAsMzML(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"), false);

            var reader =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"));
            reader.LoadAllStaticData(testFilteringParams);

            int expNumPeaks = reader.GetAllScansList().First().MassSpectrum.XArray.Length;
            double expMinRatio = reader.GetAllScansList()
                .First().MassSpectrum.YArray
                .Min(p => p / reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY).Value;
            List<(double mz, double intensity)> myExpPeaks = new List<(double mz, double intensity)>();

            for (int i = 0; i < reader.GetAllScansList().First().MassSpectrum.YArray.Length; i++)
            {
                myExpPeaks.Add((reader.GetAllScansList().First().MassSpectrum.XArray[i], reader.GetAllScansList().First().MassSpectrum.YArray[i]));
            }

            Assert.That(Math.Round(myMaxIntensity, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            Assert.That(Math.Round(sumOfAllIntensities, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.SumOfAllY, 0));
            Assert.That(myPeaksOrderedByIntensity.Count == reader.GetAllScansList().First().MassSpectrum.XArray.Length);
            Assert.That(expMinRatio >= minRatio);
            Assert.That(!myExpPeaks.Except(myPeaksOrderedByIntensity).Any());
            Assert.That(!myPeaksOrderedByIntensity.Except(myExpPeaks).Any());
        }

        [Test]
        [TestCase(null, null, null, null, null, true, false, 1000)] // no filtering return all peaks
        [TestCase(200, null, null, null, null, true, false, 200)] // top 200 peaks only
        [TestCase(null, 0.01, null, null, null, true, false, 990)] // peaks above intensity ratio
        [TestCase(200, null, null, 2, null, true, false, 400)] // top 200 peaks in each of two windows
        [TestCase(200, null, null, 20, null, true, false, 1000)] // top 200 peaks in each of twenty windows exceeds actual peak count so return all peaks
        [TestCase(200, null, 500, null, null, true, false, 400)] // top 200 peaks in each 500 Da Window
        [TestCase(200, null, 100, null, null, true, false, 1000)] // top 200 peaks in each 100 Da Window exceeds actual peak count so return all peaks
        [TestCase(null, null, null, null, true, true, false, 1000)] // no filtering return all peaks max intensity of 100
        [TestCase(200, null, 500, 20, null, true, false, 400)] // nominal width takes precedent over number of windows
        public static void TestPeakTrimmingVarietyPack(int? numberOfPeaksToKeepPerWindow, double? minimumAllowedIntensityRatioToBasePeak, int? nominalWindowWidthDaltons, int? numberOfWindows, bool normalize, bool applyTrimminToMs1, bool applyTrimmingToMsMs, int expectedPeakCount)
        {
            var testFilteringParams = new FilteringParams(numberOfPeaksToKeepPerWindow, minimumAllowedIntensityRatioToBasePeak, nominalWindowWidthDaltons, numberOfWindows, normalize, applyTrimminToMs1, applyTrimmingToMsMs);
            List<(double mz, double intensity)> myPeaks = new List<(double mz, double intensity)>();

            for (int mz = 1; mz <= 1000; mz++)
            {
                myPeaks.Add((mz, mz));
            }

            double myMaxIntensity = myPeaks.Max(p => p.intensity);

            double[] intensities1 = myPeaks.Select(p => p.intensity).ToArray();
            double[] mz1 = myPeaks.Select(p => p.mz).ToArray();

            MzSpectrum massSpec1 = new MzSpectrum(mz1, intensities1, false);
            MsDataScan[] scans = new MsDataScan[]{
                new MsDataScan(massSpec1, 1, 1, true, Polarity.Positive, 1, new MzRange(1, 1000), "f", MZAnalyzerType.Orbitrap, massSpec1.SumOfAllY, null, null, "1")
            };
            FakeMsDataFile f = new FakeMsDataFile(scans);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(f, Path.Combine(TestContext.CurrentContext.TestDirectory, "variety_mzml.mzML"), false);

            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "variety_mzml.mzML"));
            reader.LoadAllStaticData(testFilteringParams);
            int expNumPeaks = reader.GetAllScansList().First().MassSpectrum.XArray.Length;
            double expMinRatio = reader.GetAllScansList().First().MassSpectrum.YArray.Min(p => p / reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY).Value;
            List<(double mz, double intensity)> myExpPeaks = new List<(double mz, double intensity)>();

            for (int i = 0; i < reader.GetAllScansList().First().MassSpectrum.YArray.Length; i++)
            {
                myExpPeaks.Add((reader.GetAllScansList().First().MassSpectrum.XArray[i], reader.GetAllScansList().First().MassSpectrum.YArray[i]));
            }

            Assert.AreEqual(expectedPeakCount, myExpPeaks.Count());

            if (normalize)
            {
                Assert.AreEqual(50.0d, Math.Round(reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            }
            else
            {
                Assert.That(Math.Round(myMaxIntensity, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            }

            if (minimumAllowedIntensityRatioToBasePeak != null && minimumAllowedIntensityRatioToBasePeak > 0)
            {
                Assert.That(expMinRatio >= minimumAllowedIntensityRatioToBasePeak.Value);
            }
            else
            {
                Assert.That(expMinRatio >= 0);
            }
        }

        [Test]
        public static void TestPeakTrimmingWithTooManyWindows()
        {
            Random rand = new Random();
            int numPeaks = 200;
            double minRatio = 0.01;
            int numWindows = 10;

            var testFilteringParams = new FilteringParams(numPeaks, minRatio, null, numWindows, false, true, true);
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

            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "mzml.mzML"));
            reader.LoadAllStaticData(testFilteringParams);

            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.YofPeakWithHighestY.Value, 0));
            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.SumOfAllY, 0));
            Assert.That(1 == reader.GetAllScansList().First().MassSpectrum.XArray.Length);
            Assert.That(Math.Round(myPeaks[0].mz, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.XArray[0], 0));
            Assert.That(Math.Round(myPeaks[0].intensity, 0) == Math.Round(reader.GetAllScansList().First().MassSpectrum.YArray[0], 0));
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

            var reader =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "mzmlWithEmptyScan.mzML"));
            reader.LoadAllStaticData();

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, Path.Combine(TestContext.CurrentContext.TestDirectory, "mzmlWithEmptyScan2.mzML"), false);

            var testFilteringParams = new FilteringParams(200, 0.01, null, 5, false, true, true);
            var reader2 = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                "mzmlWithEmptyScan2.mzML"));
            reader2.LoadAllStaticData(testFilteringParams);
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

            var reader =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "asdfefsf.mzML"));
            reader.LoadAllStaticData();

            Assert.AreEqual(MZAnalyzerType.Orbitrap, reader.GetAllScansList().First().MzAnalyzer);
            Assert.AreEqual(MZAnalyzerType.IonTrap3D, reader.GetAllScansList().Last().MzAnalyzer);
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

            UsefulProteomicsDatabases.Loaders.LoadElements();
        }

        [Test]
        public void LoadMzmlTest()
        {
            Assert.Throws<AggregateException>(() =>
            {
                var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "DataFiles", "tiny.pwiz.1.1.mzML"));
                reader.LoadAllStaticData();
            }, "Reading profile mode mzmls not supported");
        }

        [Test]
        public void EliminateZeroIntensityPeaksFromMzmlOnFileLoad()
        {
            //create an mzML file with zero intensity peaks in both MS1 and MS2 scans
            MzSpectrum spectrum1 = new MzSpectrum(new double[] { 1, 2, 3 }, new double[] { 0, 1, 2 }, false);
            MzSpectrum spectrum2 = new MzSpectrum(new double[] { 2, 3, 4 }, new double[] { 0, 2, 4 }, false);
            MsDataScan[] scans = new MsDataScan[2];
            scans[0] = new MsDataScan(spectrum1, 1, 1, true, Polarity.Positive, 1.0, new MzRange(300, 2000), "scan filter", MZAnalyzerType.Unknown, spectrum1.SumOfAllY, null, null, null);
            scans[1] = new MsDataScan(spectrum2, 2, 2, true, Polarity.Positive, 1, new MzRange(300, 2000), "scan filter", MZAnalyzerType.Unknown, spectrum2.SumOfAllY, 1, new double[,] { }, "nativeId", 2, 1, 1, 1, 1, DissociationType.Unknown, 1, 2, null);

            var testMsDataFile = new GenericMsDataFile(scans, new SourceFile("no nativeID format", "mzML format",
                    null, null, null));

            //check that scans have zero intensities
            Assert.IsTrue(testMsDataFile.GetAllScansList()[0].MassSpectrum.YArray.Contains(0));
            Assert.IsTrue(testMsDataFile.GetAllScansList()[1].MassSpectrum.YArray.Contains(0));

            //write the mzML file with zero intensity scans
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(testMsDataFile, "mzMLWithZeros.mzML", false);

            var reader = MsDataFileReader.GetDataFile("mzMLWithZeros.mzML");
            reader.LoadAllStaticData();
            //read the mzML file. zero intensity peaks should be eliminated during read

            //insure that read scans contain no zero intensity peaks
            Assert.IsFalse(reader.GetAllScansList()[0].MassSpectrum.YArray.Contains(0));
            Assert.IsFalse(reader.GetAllScansList()[1].MassSpectrum.YArray.Contains(0));

            File.Delete("mzMLWithZeros.mzML");
        }

        [Test]
        public void LoadMzmlFromConvertedMGFTest()
        {
            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "tester.mzML"));
            reader.LoadAllStaticData();

            var ya = reader.GetOneBasedScan(1).MassSpectrum;
            Assert.AreEqual(192, ya.Size);
            var ya2 = reader.GetOneBasedScan(3).MassSpectrum;
            Assert.AreEqual(165, ya2.Size);
            var ya3 = reader.GetOneBasedScan(5).MassSpectrum;
            Assert.AreEqual(551, ya3.Size);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, "CreateFileFromConvertedMGF.mzML", false);

            var reader2 = MsDataFileReader.GetDataFile(@"CreateFileFromConvertedMGF.mzML");
            reader2.LoadAllStaticData();

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader2, "CreateFileFromConvertedMGF2.mzML", false);
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

            var reader = MsDataFileReader.GetDataFile(@"argh.mzML");
            reader.LoadAllStaticData();
            reader.GetOneBasedScan(2);

            Assert.AreEqual(1, reader.GetClosestOneBasedSpectrumNumber(1));
            Assert.AreEqual(2, reader.GetClosestOneBasedSpectrumNumber(2));

            var newFirstValue = reader.GetOneBasedScan(1).MassSpectrum.FirstX;
            Assert.AreEqual(oldFirstValue.Value, newFirstValue.Value, 1e-9);

            var secondScan2 = reader.GetOneBasedScan(2);

            Assert.AreEqual(1, secondScan2.IsolationRange.Maximum - secondScan2.IsolationRange.Minimum);

            secondScan2.MassSpectrum.ReplaceXbyApplyingFunction((a) => 44);
            Assert.AreEqual(44, secondScan2.MassSpectrum.LastX);
        }

        [Test]
        public void LocalTestForMzmlEquivalence()
        {
            MsDataFile j14_Averaged = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\InstabilityInvestigation\02-17-20_jurkat_td_rep2_fract4-calib-averaged_1-14.mzML");
            MsDataFile j16_Averaged = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\InstabilityInvestigation\02-17-20_jurkat_td_rep2_fract4-calib-averaged_1-16.mzML");

            j14_Averaged.LoadAllStaticData();
            j16_Averaged.LoadAllStaticData();

            List<MsDataScan> variableScans = new();
            List<double> mzDiscrepancies = new();

            for(int i = 0; i < j14_Averaged.Scans.Length; i++)
            {
                if (!j14_Averaged.Scans[i].MassSpectrum.Equals(j16_Averaged.Scans[i].MassSpectrum))
                {
                    variableScans.Add(j14_Averaged.Scans[i]);
                    variableScans.Add(j16_Averaged.Scans[i]);
                    if (j14_Averaged.Scans[i].MassSpectrum.Size != j16_Averaged.Scans[i].MassSpectrum.Size)
                    {
                        Console.WriteLine("Spectra sized differ by: " + (j14_Averaged.Scans[i].MassSpectrum.Size - j16_Averaged.Scans[i].MassSpectrum.Size));
                        continue;
                    }
                        
                    for (int j = 0; j < j14_Averaged.Scans[i].MassSpectrum.Size; j++)
                    {
                        if (j14_Averaged.Scans[i].MassSpectrum.XArray[j] != j16_Averaged.Scans[i].MassSpectrum.XArray[j])
                            mzDiscrepancies.Add(Math.Abs(j14_Averaged.Scans[i].MassSpectrum.XArray[j] - j16_Averaged.Scans[i].MassSpectrum.XArray[j]));
                    }
                }
            }

            Console.WriteLine("Number of mismatched spectra: " + variableScans.Count / 2);
            Console.WriteLine("Total number of spectra: " + j14_Averaged.Scans.Length);
            Console.WriteLine("Median peak difference: " + mzDiscrepancies.Median());
            Console.WriteLine("Max peak difference: " + mzDiscrepancies.Max());
            Console.WriteLine("Total number of peak differences: " + mzDiscrepancies.Count);
        }

        [Test]
        public void LocalTestForMzmlEquivalenceCalibratedFiles()
        {
            MsDataFile j14_Averaged = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\InstabilityInvestigation\02-17-20_jurkat_td_rep2_fract4-calib_1-14.mzML");
            MsDataFile j16_Averaged = MsDataFileReader.GetDataFile(@"D:\MetaMorpheusVignette\InstabilityInvestigation\02-17-20_jurkat_td_rep2_fract4-calib_1-16.mzML");

            j14_Averaged.LoadAllStaticData();
            j16_Averaged.LoadAllStaticData();

            List<MsDataScan> variableScans = new();
            List<double> mzDiscrepancies = new();

            for (int i = 0; i < j14_Averaged.Scans.Length; i++)
            {
                if (!j14_Averaged.Scans[i].MassSpectrum.Equals(j16_Averaged.Scans[i].MassSpectrum))
                {
                    variableScans.Add(j14_Averaged.Scans[i]);
                    variableScans.Add(j16_Averaged.Scans[i]);
                    if (j14_Averaged.Scans[i].MassSpectrum.Size != j16_Averaged.Scans[i].MassSpectrum.Size)
                    {
                        Console.WriteLine("Spectra sized differ by: " + (j14_Averaged.Scans[i].MassSpectrum.Size - j16_Averaged.Scans[i].MassSpectrum.Size));
                        continue;
                    }

                    for (int j = 0; j < j14_Averaged.Scans[i].MassSpectrum.Size; j++)
                    {
                        if (j14_Averaged.Scans[i].MassSpectrum.XArray[j] != j16_Averaged.Scans[i].MassSpectrum.XArray[j])
                            mzDiscrepancies.Add(Math.Abs(j14_Averaged.Scans[i].MassSpectrum.XArray[j] - j16_Averaged.Scans[i].MassSpectrum.XArray[j]));
                    }
                }
            }

            Console.WriteLine("Number of mismatched spectra: " + variableScans.Count / 2);
            Console.WriteLine("Total number of spectra: " + j14_Averaged.Scans.Length);
            Console.WriteLine("Median peak difference: " + mzDiscrepancies.Median());
            Console.WriteLine("Max peak difference: " + mzDiscrepancies.Max());
            Console.WriteLine("Total number of peak differences: " + mzDiscrepancies.Count);
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
            var fakeMzml =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "what.mzML"));
            fakeMzml.LoadAllStaticData();

            FakeMsDataFile fakeFile1 = new FakeMsDataFile(scans1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fakeFile1, Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"), false);
            var fakeMzml1 =
                MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "what1.mzML"));
            fakeMzml1.LoadAllStaticData();

            Assert.AreEqual(3, fakeMzml.GetAllScansList().ElementAt(5).OneBasedPrecursorScanNumber);
            Assert.AreEqual(1, fakeMzml1.GetAllScansList().ElementAt(3).OneBasedPrecursorScanNumber);
        }

        [Test]
        [TestCase("tester.mzml")]
        [TestCase("SmallCalibratibleYeast.mzml")]
        [TestCase("small.raw", true)]
        [TestCase("small.raw", false)]
        [TestCase("testFileWMS2.raw", true)]
        [TestCase("testFileWMS2.raw", false)]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw", true)]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw", false)]
        public static void TestDynamicMzml(string fileName, bool writeIndexed = false)
        {
            if (Path.GetExtension(fileName).ToUpper() == ".RAW")
            {
                string append = writeIndexed ? "_indexed" : "_unindexed";
                string rawPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);
                var thermoReader = MsDataFileReader.GetDataFile(rawPath);
                thermoReader.LoadAllStaticData();

                string mzmlFilename = Path.GetFileNameWithoutExtension(fileName) + append + ".mzML";
                fileName = mzmlFilename;
                string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

                Task task = Task.Run(() => MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(thermoReader, mzmlPath, writeIndexed));
                while (!task.IsCompleted) { }
            }


            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData();
            reader.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in reader.GetAllScansList())
            {
                MsDataScan dynamicScan = reader.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);

                if (!double.IsNaN(dynamicScan.RetentionTime) || !double.IsNaN(staticScan.RetentionTime))
                {
                    Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                }

                Assert.That(dynamicScan.Polarity == staticScan.Polarity);

                if (!double.IsNaN(staticScan.ScanWindowRange.Minimum) || !double.IsNaN(staticScan.ScanWindowRange.Maximum)
                    || !double.IsNaN(dynamicScan.ScanWindowRange.Minimum) || !double.IsNaN(dynamicScan.ScanWindowRange.Maximum))
                {
                    Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                    Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                }

                Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.TotalIonCurrent == staticScan.TotalIonCurrent);
                Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);

                if (dynamicScan.IsolationWidth != null || staticScan.IsolationWidth != null)
                {
                    if (!double.IsNaN(dynamicScan.IsolationWidth.Value) || !double.IsNaN(staticScan.IsolationWidth.Value))
                    {
                        Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                    }
                }

                Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);
                }
            }
            reader.CloseDynamicConnection();
        }

        [Test]
        public static void TestDynamicMzmlWithPeakFiltering()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            FilteringParams filteringParams = new FilteringParams(200, 0.01, 1, null, false, true, true);

            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(filteringParams);
            reader.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in reader.GetAllScansList())
            {
                if (staticScan.OneBasedScanNumber > 10)
                {
                    // only the first 10 scans are checked because checking the entire file takes 7min on AppVeyor
                    break;
                }

                MsDataScan dynamicScan = reader.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber, filteringParams);

                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);

                if (!double.IsNaN(dynamicScan.RetentionTime) || !double.IsNaN(staticScan.RetentionTime))
                {
                    Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                }

                Assert.That(dynamicScan.Polarity == staticScan.Polarity);

                if (!double.IsNaN(staticScan.ScanWindowRange.Minimum) || !double.IsNaN(staticScan.ScanWindowRange.Maximum)
                    || !double.IsNaN(dynamicScan.ScanWindowRange.Minimum) || !double.IsNaN(dynamicScan.ScanWindowRange.Maximum))
                {
                    Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                    Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                }

                Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.TotalIonCurrent == staticScan.TotalIonCurrent);
                Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);

                if (dynamicScan.IsolationWidth != null || staticScan.IsolationWidth != null)
                {
                    if (!double.IsNaN(dynamicScan.IsolationWidth.Value) || !double.IsNaN(staticScan.IsolationWidth.Value))
                    {
                        Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                    }
                }

                Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);
                }
            }
        }

        [Test]
        /// This test reads an .mzML file where the intensities and m/z values are encoded in compressed 64-bit, and also an identical
        /// file that is encoded in mostly uncompressed 32-bit. Scan #73 in the latter file has its intensities encoded in compressed 64-bit.
        /// The intensities and m/z values for both of these files should be identical (sans rounding issues).
        /// Additionally, scan 73 has its retention time in seconds, rather than minutes.
        public static void TestDynamicMzmlWithMixedBits()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            string filePath2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast_mixed_bits.mzml");

            var reader1 = MsDataFileReader.GetDataFile(filePath);
            reader1.LoadAllStaticData();
            var reader2 = MsDataFileReader.GetDataFile(filePath2);
            reader2.LoadAllStaticData();

            reader1.InitiateDynamicConnection();
            reader2.InitiateDynamicConnection();

            foreach (MsDataScan staticScan in reader1.GetAllScansList())
            {
                // only the scans 72-74 are checked because checking the entire file takes a long time on AppVeyor
                if (staticScan.OneBasedScanNumber < 72)
                {
                    continue;
                }
                if (staticScan.OneBasedScanNumber > 74)
                {
                    break;
                }

                MsDataScan staticScan2 = reader2.GetOneBasedScan(staticScan.OneBasedScanNumber);
                MsDataScan dynamicScan1 = reader1.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);
                MsDataScan dynamicScan2 = reader2.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.That(staticScan2.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan1.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan2.OneBasedScanNumber == staticScan.OneBasedScanNumber);

                Assert.That(staticScan2.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan1.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan2.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);

                Assert.That(staticScan2.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);
                Assert.That(dynamicScan1.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);
                Assert.That(dynamicScan2.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double staticMz2 = staticScan2.MassSpectrum.XArray[i];
                    double staticIntensity2 = staticScan2.MassSpectrum.YArray[i];

                    Assert.That(staticMz2 == staticMz);
                    Assert.That(staticIntensity2 == staticIntensity);

                    double dynamicMz = dynamicScan1.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan1.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);

                    double dynamicMz2 = dynamicScan2.MassSpectrum.XArray[i];
                    double dynamicIntensity2 = dynamicScan2.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz2 == staticMz);
                    Assert.That(dynamicIntensity2 == staticIntensity);
                }
            }
        }

        [Test]
        public static void TestEthcdReading()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.mzML");
            var reader = MsDataFileReader.GetDataFile(filePath);
            reader.LoadAllStaticData(null, 1);
            var unknownScan = reader.GetOneBasedScan(3);
            Assert.That(unknownScan.DissociationType == DissociationType.Unknown);
            var hcdScan = reader.GetOneBasedScan(5);
            Assert.That(hcdScan.DissociationType == DissociationType.HCD);
            var ethcdScan = reader.GetOneBasedScan(6);
            Assert.That(ethcdScan.DissociationType == DissociationType.EThcD);
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
        }
    }
}