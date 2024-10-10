using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
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
        public void TestMsAlign_DynamnicConnectionAndHeaderComponents(string path, int oneBasedScanNumber, int msnOrder,
            double retentionTime,
            int peakCount, int? oneBasePrecursorScanNumber, double? precursorMz, int? precursorCharge,
            double? precursorIntensity, DissociationType? dissociationType,
            double? precursorIsolationMzStart, double? precursorIsolationMzEnd)
        {
            var file = MsAlignTestFiles[path];
            file.InitiateDynamicConnection();

            var scanToTest = file.GetOneBasedScanFromDynamicConnection(oneBasedScanNumber);

            Assert.That(oneBasedScanNumber, Is.EqualTo(scanToTest.OneBasedScanNumber));
            Assert.That(msnOrder, Is.EqualTo(scanToTest.MsnOrder));
            Assert.That(retentionTime, Is.EqualTo(scanToTest.RetentionTime).Within(0.001));
            Assert.That(peakCount, Is.EqualTo(scanToTest.MassSpectrum.Size));
            Assert.That(oneBasePrecursorScanNumber, Is.EqualTo(scanToTest.OneBasedPrecursorScanNumber));
            Assert.That(precursorMz, Is.EqualTo(scanToTest.SelectedIonMZ).Within(0.001));
            Assert.That(precursorCharge, Is.EqualTo(scanToTest.SelectedIonChargeStateGuess));
            Assert.That(precursorIntensity, Is.EqualTo(scanToTest.SelectedIonIntensity).Within(0.001));
            Assert.That(dissociationType, Is.EqualTo(scanToTest.DissociationType));

            if (msnOrder is 1)
            {
                Assert.That(precursorIsolationMzStart, Is.EqualTo(scanToTest.IsolationRange));
                Assert.That(precursorIsolationMzEnd, Is.EqualTo(scanToTest.IsolationRange));
            }
            else
            {
                Assert.That(precursorIsolationMzStart, Is.EqualTo(scanToTest.IsolationRange.Minimum));
                Assert.That(precursorIsolationMzEnd, Is.EqualTo(scanToTest.IsolationRange.Maximum));
            }

            file.CloseDynamicConnection();
        }

        [Test]
        public void GetOneBasedScanFromDynamicConnection_ExistingScanNumber_ReturnsMsDataScan()
        {
            // Arrange
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            // Act
            var scan = msAlign.GetOneBasedScanFromDynamicConnection(1);

            // Assert
            Assert.That(scan, Is.Not.Null);
            Assert.That(1, Is.EqualTo(scan.OneBasedScanNumber));
        }

        [Test]
        public void GetOneBasedScanFromDynamicConnection_NonExistingScanNumber_ThrowsException()
        {
            // Arrange
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            // Act and Assert
            Assert.Throws<MzLibException>(() => msAlign.GetOneBasedScanFromDynamicConnection(100));
        }

        [Test]
        public void CloseDynamicConnection_ConnectionOpen_ClosesConnection()
        {
            // Arrange
            var msAlign = MsAlignTestFiles.First().Value;
            msAlign.InitiateDynamicConnection();

            // Act
            msAlign.CloseDynamicConnection();

            // Assert
            Assert.Throws<MzLibException>(() => msAlign.GetOneBasedScanFromDynamicConnection(14560790));
        }

        [Test]
        public void InitiateDynamicConnection_FileDoesNotExist_FileNotFoundExceptionThrown()
        {
            // Arrange
            var msAlign = new Ms2Align("NonExistentFile.mzML");

            // Act & Assert
            Assert.Throws<FileNotFoundException>(() => msAlign.InitiateDynamicConnection());
        }

        
    }
}
