using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    internal class TestMsAlign
    {
        public static Dictionary<string, MsAlign> MsAlignTestFiles { get; set; }

        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            MsAlignTestFiles = new Dictionary<string, MsAlign>
            {
                {
                    "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign",
                    new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign"))
                },
                {
                    "Ms1Align_IsoDec1.0.0_ms1.msalign",
                    new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms1Align_IsoDec1.0.0_ms1.msalign"))
                },
                {
                    "Ms1Align_TopFDv1.6.2_ms1.msalign",
                    new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms1Align_TopFDv1.6.2_ms1.msalign"))
                },
                {
                    "Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign",
                    new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign"))
                },
                {
                    "Ms2Align_IsoDec1.0.0_ms2.msalign",
                    new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms2Align_IsoDec1.0.0_ms2.msalign"))
                },
                {
                    "Ms2Align_TopFDv1.6.2_ms2.msalign",
                    new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        @"FileReadingTests\ExternalFileTypes", "Ms2Align_TopFDv1.6.2_ms2.msalign"))
                }
            };

            foreach (var msAlign in MsAlignTestFiles.Values)
                msAlign.LoadAllStaticData();
        }

        [Test]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 7, 1)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 12, 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 22, 2)]
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 42, 2)]
        public void TestMsAlignResultsLoadsAndCountCorrect(string path, int recordCount, int expectedMsAlignNumber)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", path);

            MsAlign file = expectedMsAlignNumber == 1
                ? new Ms1Align(filePath)
                : new Ms2Align(filePath);

            file.LoadAllStaticData();
            Assert.That(file.Count(), Is.EqualTo(recordCount));
            Assert.That(file.DefaultMsnOrder, Is.EqualTo(expectedMsAlignNumber));
        }

        [Test]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 7, 1)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 12, 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 22, 2)]
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 42, 2)]
        public void TestMsAlignResultsLoadsAndCountCorrect_GenericMsDataFileReader(string path, int recordCount,
            int expectedMsAlignNumber)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", path);

            MsDataFile file = MsDataFileReader.GetDataFile(filePath);
            file.LoadAllStaticData();
            Assert.That(file.Count(), Is.EqualTo(recordCount));

            Assert.That(file is MsAlign msa);
            if (expectedMsAlignNumber == 1)
                Assert.That(file is Ms1Align);
            else
                Assert.That(file is Ms2Align);

            Assert.That(((MsAlign)file).DefaultMsnOrder, Is.EqualTo(expectedMsAlignNumber));
        }

        [Test]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 7, 1)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 12, 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 22, 2)]
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 42, 2)]
        public void TestMsAlignResultsLoadsAndCountCorrect_GenericFileReader(string path, int recordCount,
            int expectedMsAlignNumber)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", path);

            MsDataFile file = FileReader.ReadFile<MsDataFileToResultFileAdapter>(filePath);
            Assert.That(file.Count(), Is.EqualTo(recordCount));
            Assert.That(file.Scans.Length, Is.EqualTo(recordCount));
            Assert.That(file.NumSpectra, Is.EqualTo(recordCount));
            Assert.That(file.SourceFile.MassSpectrometerFileFormat, Contains.Substring("align"));
            Assert.That(file.SourceFile.MassSpectrometerFileFormat,
                Contains.Substring(expectedMsAlignNumber.ToString()));
        }

        [Test]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 1, 1, 0.0035, 190, null, null, null, null, null, null,
            null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 23, 1, 0.4035, 172, null, null, null, null, null,
            null, null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 26, 1, 0.4412, 181, null, null, null, null, null,
            null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 1, 1, 0.0205, 3, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1, 7.4313, 16, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 1, 1, 0.0021, 13, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1, 0.9861, 14, null, null, null, null, null, null, null)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 2, 2, 0.87 / 60.0, 13, 1, 344.721008, 2, 85166.05, DissociationType.HCD, 344.721008 - 1.5, 344.721008 + 1.5)] 
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 5, 2, 3.73 / 60.0, 3, 1, 350.708160, 2, 9323.61, DissociationType.HCD, 350.708160 - 1.5, 350.708160 + 1.5)] 
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 580, 2, 970.65 / 60.0, 2, 579, 0.00000, 0, 0.00, DissociationType.HCD, -2, 2)] 
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 1669, 2, 2726.74 / 60.0, 0, 1665, 1239.30464, 8, 1329780.59, DissociationType.HCD, 1239.30464 - 2, 1239.30464 + 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1366, 2, 832.915128799 / 60.0, 24, 1365, 1005.95751953125, 5, 619272.56, DissociationType.HCD, 1005.1555541753769, 1006.3555542230606)] // First scan
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1398, 2, 862.421989215 / 60, 22, 1395, 973.6828002929688, 6, 122901736.0, DissociationType.HCD, 973.0828002691269, 974.2828003168106)] // Last scan
        public void TestMsAlign_DynamnicConnectionAndHeaderComponents(string path, int oneBasedScanNumber, int msnOrder,
            double retentionTime,
            int peakCount, int? oneBasePrecursorScanNumber, double? precursorMz, int? precursorCharge,
            double? precursorIntensity, DissociationType? dissociationType,
            double? precursorIsolationMzStart, double? precursorIsolationMzEnd)
        {
            var file = MsAlignTestFiles[path];
            file.InitiateDynamicConnection();

            var scanToTest = file.GetOneBasedScanFromDynamicConnection(oneBasedScanNumber);

            Assert.That(scanToTest.OneBasedScanNumber, Is.EqualTo(oneBasedScanNumber));
            Assert.That(scanToTest.MsnOrder, Is.EqualTo(msnOrder));
            Assert.That(scanToTest.RetentionTime, Is.EqualTo(retentionTime).Within(0.001));
            Assert.That(scanToTest.MassSpectrum.Size, Is.EqualTo(peakCount));
            Assert.That(scanToTest.OneBasedPrecursorScanNumber, Is.EqualTo(oneBasePrecursorScanNumber));
            Assert.That(scanToTest.SelectedIonMZ, Is.EqualTo(precursorMz).Within(0.001));
            Assert.That(scanToTest.SelectedIonChargeStateGuess, Is.EqualTo(precursorCharge));
            Assert.That(scanToTest.SelectedIonIntensity, Is.EqualTo(precursorIntensity).Within(0.001));
            Assert.That(scanToTest.DissociationType, Is.EqualTo(dissociationType));

            if (msnOrder is 1)
            {
                Assert.That(scanToTest.IsolationRange, Is.EqualTo(precursorIsolationMzStart));
                Assert.That(scanToTest.IsolationRange, Is.EqualTo(precursorIsolationMzEnd));
            }
            else
            {
                Assert.That(scanToTest.IsolationRange.Minimum, Is.EqualTo(precursorIsolationMzStart));
                Assert.That(scanToTest.IsolationRange.Maximum, Is.EqualTo(precursorIsolationMzEnd));
            }

            file.CloseDynamicConnection();
        }

        [Test]
        public void GetOneBasedScanFromDynamicConnection_ExistingScanNumber_ReturnsMsDataScan()
        {
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            var scan = msAlign.GetOneBasedScanFromDynamicConnection(1);

            Assert.That(scan, Is.Not.Null);
            Assert.That(1, Is.EqualTo(scan.OneBasedScanNumber));
        }

        [Test]
        public void GetOneBasedScanFromDynamicConnection_NonExistingScanNumber_ThrowsException()
        {
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            Assert.Throws<MzLibException>(() => msAlign.GetOneBasedScanFromDynamicConnection(100));
        }

        [Test]
        public void CloseDynamicConnection_ConnectionOpen_ClosesConnection()
        {
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            msAlign.CloseDynamicConnection();

            Assert.Throws<MzLibException>(() => msAlign.GetOneBasedScanFromDynamicConnection(14560790));
        }

        [Test]
        public void InitiateDynamicConnection_FileDoesNotExist_FileNotFoundExceptionThrown()
        {
            var msAlign = new Ms2Align("NonExistentFile.mzML");

            Assert.Throws<FileNotFoundException>(() => msAlign.InitiateDynamicConnection());
        }

        
    }
}
