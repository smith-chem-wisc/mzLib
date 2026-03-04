// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestDataFile.cs) is part of MassSpectrometry.Tests.
//
// MassSpectrometry.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry.Tests. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMsDataFile
    {
        private MzSpectrum _mzSpectrumA;
        private FakeMsDataFile myMsDataFile;
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzSpectrum(mz, intensities, false);

            var peptide = new Peptide("KQEEQMETEQQNKDEGK");

            MzSpectrum MS1 = CreateSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
            MzSpectrum MS2 = CreateMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

            MsDataScan[] Scans = new MsDataScan[2];
            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, 1, null, "scan=1");

            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, 1, null, "scan=2", 693.9892, 3, .3872, 693.99, 1, DissociationType.Unknown, 1, 693.6550);

            myMsDataFile = new FakeMsDataFile(Scans);
        }

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
        public void SpectrumCount()
        {
            Assert.AreEqual(15, _mzSpectrumA.Size);
        }

        [Test]
        public void SpectrumFirstMZ()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumLastMZ()
        {
            Assert.AreEqual(723.35345, _mzSpectrumA.LastX);
        }

        [Test]
        public void DataFileTest()
        {
            MsDataScan theSpectrum = new MsDataScan(_mzSpectrumA, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 1000), "fake scan filter", MZAnalyzerType.Unknown, _mzSpectrumA.SumOfAllY, 1, null, "scan=1");

            MsDataScan[] theList = new MsDataScan[1];

            theList[0] = theSpectrum;

            FakeMsDataFile thefile = new FakeMsDataFile(theList);

            var theOneBasedScan = thefile.GetOneBasedScan(1);

            Assert.AreEqual("Scan #1", theOneBasedScan.ToString());

            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);
            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);

            Assert.AreEqual(1, thefile.NumSpectra);
            Assert.AreEqual(1, thefile.NumSpectra);

            Assert.IsTrue(thefile.GetOneBasedScan(1).IsCentroid);

            foreach (var ok in thefile)
            {
                Assert.AreEqual(300, ok.ScanWindowRange.Minimum, 1e-9);
                Assert.AreEqual(1000, ok.ScanWindowRange.Maximum, 1e-9);
            }

            int ok1 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(0, 2))
                ok1 += 1;

            Assert.AreEqual(1, ok1);

            int ok2 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(2, 4))
                ok2 += 1;

            Assert.AreEqual(0, ok2);

            int ok3 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(-4, -2))
                ok3 += 1;

            Assert.AreEqual(0, ok3);
        }

        [Test]
        public void TestFunctionsOfMsDataScan()
        {
            MsDataScan theSpectrum = new MsDataScan(_mzSpectrumA, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 1000), "fake scan filter", MZAnalyzerType.Unknown, _mzSpectrumA.SumOfAllY, 1, null, "scan=1");
            List<IsotopicEnvelope> isolatedMassesAndCharges = theSpectrum.GetIsolatedMassesAndCharges(_mzSpectrumA, 1, 10, 10, 1).ToList();
            Assert.AreEqual(0, isolatedMassesAndCharges.Count); //Isolation range is null, so we get an empty set

            Assert.Throws<MzLibException>(() => theSpectrum.RefineSelectedMzAndIntensity(_mzSpectrumA)); //no isolation Mz throws error 

            theSpectrum.SetOneBasedPrecursorScanNumber(6);
            Assert.AreEqual(6, theSpectrum.OneBasedPrecursorScanNumber);

            theSpectrum.SetNativeID("bubba");
            Assert.AreEqual("bubba", theSpectrum.NativeId);

            theSpectrum.SetIsolationMz(42);
            Assert.AreEqual(42, theSpectrum.IsolationMz);

            theSpectrum.SetMsnOrder(2);
            Assert.AreEqual(2, theSpectrum.MsnOrder);

            theSpectrum.SetIsolationRange(400, 2000);
            Assert.AreEqual(400, theSpectrum.IsolationRange.Minimum);
            Assert.AreEqual(2000, theSpectrum.IsolationRange.Maximum);
        }

        [Test]
        public void MoreMsDataFilesTests()
        {
            GenericMsDataFile fakeDataFile = new GenericMsDataFile(new MsDataScan[1],
                new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null));
            Assert.AreEqual(1, fakeDataFile.NumSpectra);
            Assert.AreEqual("scan number only nativeID format", fakeDataFile.SourceFile.NativeIdFormat);
            Assert.AreEqual("mzML format", fakeDataFile.SourceFile.MassSpectrometerFileFormat);
            Assert.IsNull(fakeDataFile.SourceFile.CheckSum);
            Assert.AreEqual("SHA-1", fakeDataFile.SourceFile.FileChecksumType);
            Assert.IsNull(fakeDataFile.SourceFile.Id);
        }

        [Test]
        public void TestAMoreRealFile()
        {
            var theScan = myMsDataFile.GetOneBasedScan(2);
            Assert.AreEqual(1, theScan.IsolationRange.Width);
            Assert.AreEqual(DissociationType.Unknown, theScan.DissociationType);
            Assert.AreEqual(693.99, theScan.IsolationMz);
            Assert.AreEqual(1, theScan.IsolationRange.Maximum - theScan.IsolationRange.Minimum);
            Assert.AreEqual(1, theScan.OneBasedPrecursorScanNumber);
            Assert.AreEqual(3, theScan.SelectedIonChargeStateGuess.Value);

            var precursorScan = myMsDataFile.GetOneBasedScan(theScan.OneBasedPrecursorScanNumber.Value);
            theScan.RefineSelectedMzAndIntensity(precursorScan.MassSpectrum);
            Assert.AreEqual(.32872, (double)theScan.SelectedIonIntensity, 0.01);
            Assert.AreEqual(693.9892, (double)theScan.SelectedIonMZ, 0.01);
            Assert.AreEqual(693.655, (double)theScan.SelectedIonMonoisotopicGuessMz, 0.001);

            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);

            Assert.AreEqual(myMsDataFile.GetOneBasedScan(2).MassSpectrum.CopyTo2DArray()[0, 0], myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);

            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);

            theScan.ComputeMonoisotopicPeakIntensity(precursorScan.MassSpectrum);

            theScan.TransformMzs(b => 0, b => 0);

            Assert.AreEqual("Scan #2", myMsDataFile.GetOneBasedScan(2).ToString());

            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);
            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);
            Assert.AreEqual(0, theScan.SelectedIonMZ);

            List<MsDataScan> a = myMsDataFile.GetAllScansList();
            foreach (var b in a)
                Assert.IsFalse(b.IsCentroid);
            foreach (var b in myMsDataFile)
                Assert.AreEqual(Polarity.Positive, b.Polarity);
        }

        [Test]
        public static void TestXicExtraction()
        {
            string dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            var reader = MsDataFileReader.GetDataFile(dataFilePath);
            reader.LoadAllStaticData();

            var peptide = new PeptideWithSetModifications("KAPAGGAADAAAK", new Dictionary<string, Modification>());
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
        [Test]
        public static void SkipEmptyScansWhenReadingMzml()
        {
            //this original mzML has four scans including 1 MS1 and three MS2s. The second MS2 scan does
            //not have mz or intensity values. We skip this scan when reading.
            string dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", @"badScan7192.mzml");

            var reader = MsDataFileReader.GetDataFile(dataFilePath);
            reader.LoadAllStaticData();

            MsDataScan[] ms1Scans = reader.GetMS1Scans().ToArray();
            MsDataScan[] allScans = reader.GetAllScansList().ToArray();

            Assert.AreEqual(1, ms1Scans.Length);
            Assert.AreEqual(3, allScans.Length);
            List<int> expectedScanNumbers = new() { 1, 2, 4 };
            CollectionAssert.AreEquivalent(expectedScanNumbers, allScans.Select(s => s.OneBasedScanNumber).ToList());
            Assert.AreEqual(1, allScans[1].OneBasedPrecursorScanNumber);

            //even with the deleted scan, the second MS2Scan should still refer to the correct MS1
            Assert.AreEqual(1, allScans[2].OneBasedPrecursorScanNumber);
        }
        private MzSpectrum CreateSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
        {
            IsotopicDistribution isodist = IsotopicDistribution.GetDistribution(f, 0.1, 0.001);
            MzSpectrum notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);

            notActuallyMzS.ReplaceXbyApplyingFunction(b => b.Mz.ToMz(1));

            List<double> allMasses = new List<double>();
            List<double> allIntensitiess = new List<double>();

            while (notActuallyMzS.FirstX > lowerBound)
            {
                for (int i = 0; i < notActuallyMzS.Size; i++)
                {
                    if (notActuallyMzS.XArray[i] > lowerBound && notActuallyMzS.XArray[i] < upperBound)
                    {
                        allMasses.Add(notActuallyMzS.XArray[i]);
                        allIntensitiess.Add(notActuallyMzS.YArray[i]);
                    }
                }
                minCharge += 1;
                notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);
                notActuallyMzS.ReplaceXbyApplyingFunction(s => s.Mz.ToMz(minCharge));
            }

            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensitiess.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);

            return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        [Test]
        public void TestMsDataFileWithoutLoadingStaticData()
        {
            string dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml");
            var dataFile = MsDataFileReader.GetDataFile(dataFilePath);

            Assert.That(!dataFile.CheckIfScansLoaded());
            var scansArray = dataFile.GetMsDataScans();
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(scansArray.Length == 142);

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var scanslist = dataFile.GetAllScansList();
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(scanslist.Count == 142);

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var ms1ScanList = dataFile.GetMS1Scans().ToList();
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(ms1ScanList.All(p => p.MsnOrder == 1));
            Assert.That(scanslist.Count == 142);

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var scan = dataFile.GetOneBasedScan(1);
            Assert.That(dataFile.CheckIfScansLoaded());

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var scanslistInTimeRange = dataFile.GetMsScansInTimeRange(0, 100).ToList();
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(scanslistInTimeRange.Count == 142);

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var scanslistInIndexRange = dataFile.GetMsScansInIndexRange(1, 10).ToList();
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(scanslistInIndexRange.Count == 10);

            dataFile = MsDataFileReader.GetDataFile(dataFilePath);
            Assert.That(!dataFile.CheckIfScansLoaded());
            var index = dataFile.GetClosestOneBasedSpectrumNumber(5);
            Assert.That(dataFile.CheckIfScansLoaded());
            Assert.That(scanslist.Count == 142);
        }

        [Test]
        public static void NegativeModeSetsCorrectCharge_FromFile()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
                "GUACUG_NegativeMode_Sliced.mzML");
            var scans = MsDataFileReader.GetDataFile(filePath).GetAllScansList();
            var ms1 = scans.FirstOrDefault(s => s.MsnOrder == 1);

            Assert.That(ms1, Is.Not.Null, "No MS1 scan found in the file.");
            Assert.That(ms1.Polarity, Is.EqualTo(Polarity.Negative), "MS1 scan polarity is not negative.");

            var ms2 = scans.FirstOrDefault(s => s.MsnOrder == 2);
            Assert.That(ms2, Is.Not.Null, "No MS2 scan found in the file.");
            Assert.That(ms2.Polarity, Is.EqualTo(Polarity.Negative), "MS2 scan polarity is not negative.");
            Assert.That(ms2.SelectedIonChargeStateGuess.HasValue, Is.True, filePath + " does not have charge state guess for MS2 scan.");
            Assert.That(ms2.SelectedIonChargeStateGuess!.Value, Is.EqualTo(-3), "MS2 scan charge state guess is not -3.");
        }

        [TestCase(Polarity.Positive, 3, 3, TestName = "PositiveMode_PositiveCharge")]
        [TestCase(Polarity.Negative, -3, -3, TestName = "NegativeMode_NegativeCharge")]
        [TestCase(Polarity.Negative, 3, -3, TestName = "NegativeMode_PositiveCharge")]
        [TestCase(Polarity.Positive, -3, 3, TestName = "PositiveMode_NegativeCharge")]
        public static void NegativeModeSetsCorrectCharge_FromConstructedScan(Polarity polarity, int inputCharge, int expectedCharge)
        {
            var scan = new MsDataScan(new([], [], false), 1, 2, true, polarity, 2, new(100, 1000), "", MZAnalyzerType.Orbitrap, 100, 2, null, "", 120, inputCharge, 20, 120, inputCharge, DissociationType.CID, 1, 120, "30", "");

            Assert.That(scan.Polarity, Is.EqualTo(polarity), "Newly created scan has incorrect polarity.");
            Assert.That(scan.SelectedIonChargeStateGuess, Is.Not.Null, "Charge state guess for newly created scan.");
            Assert.That(scan.SelectedIonChargeStateGuess!.Value, Is.EqualTo(expectedCharge), $"Newly created scan charge state guess is not {expectedCharge}.");
        }
    }
}