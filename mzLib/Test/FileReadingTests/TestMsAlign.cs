using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Proteomics;
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
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 1, 1, 0.0035, 190, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 23, 1, 0.4035, 172, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 26, 1, 0.4412, 181, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 1, 1, 0.0205, 3, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1, 7.4313, 16, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 1, 1, 0.0021, 13, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1, 0.9861, 14, null, null, null, null, null, null, null)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 2, 2, 0.87 / 60.0, 13, 1, 344.721008, 2, 85166.05, DissociationType.HCD, 344.721008 - 1.5, 344.721008 + 1.5)] 
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 5, 2, 3.73 / 60.0, 3, 1, 350.708160, 2, 9323.61, DissociationType.HCD, 350.708160 - 1.5, 350.708160 + 1.5)] 
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 580, 2, 970.65 / 60.0, 2, 579, 0.00000, 0, 0.00, DissociationType.HCD, -2, 2)] 
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 1669, 2, 2726.74 / 60.0, 0, 1665, 1239.30464, 8, 1329780.59, DissociationType.HCD, 1239.30464 - 2, 1239.30464 + 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1366, 2, 832.915128799 / 60.0, 24, 1365, 1005.95751953125, 5, 619272.56, DissociationType.HCD, 1005.1555541753769, 1006.3555542230606)] 
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1398, 2, 862.421989215 / 60, 22, 1395, 973.6828002929688, 6, 122901736.0, DissociationType.HCD, 973.0828002691269, 974.2828003168106)] 
        public void TestMsAlign_HeaderComponents_DynamicConnection(string path, int oneBasedScanNumber, int msnOrder,
            double retentionTime,
            int peakCount, int? oneBasePrecursorScanNumber, double? precursorMz, int? precursorCharge,
            double? precursorIntensity, DissociationType? dissociationType,
            double? precursorIsolationMzStart, double? precursorIsolationMzEnd)
        {
            var file = MsAlignTestFiles[path];
            Assert.That(file.GetSourceFile().Id, Is.EqualTo(msnOrder.ToString()));
            Assert.That(file.GetSourceFile().FileName, Is.EqualTo(path.GetPeriodTolerantFilenameWithoutExtension().Replace("_ms1", "").Replace("_ms2", "")));
            Assert.That(file.SourceFile.NativeIdFormat, Is.EqualTo("no nativeID format"));
            Assert.That(file.SourceFile.MassSpectrometerFileFormat, Is.EqualTo($"ms{msnOrder}.msalign format"));

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
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 1, 1, 0.0035, 190, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 23, 1, 0.4035, 172, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 26, 1, 0.4412, 181, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 1, 1, 0.0205, 3, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_TopFDv1.6.2_ms1.msalign", 265, 1, 7.4313, 16, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 1, 1, 0.0021, 13, null, null, null, null, null, null, null)]
        [TestCase(@"Ms1Align_IsoDec1.0.0_ms1.msalign", 99, 1, 0.9861, 14, null, null, null, null, null, null, null)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 2, 2, 0.87 / 60.0, 13, 1, 344.721008, 2, 85166.05, DissociationType.HCD, 344.721008 - 1.5, 344.721008 + 1.5)]
        [TestCase(@"Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", 5, 2, 3.73 / 60.0, 3, 1, 350.708160, 2, 9323.61, DissociationType.HCD, 350.708160 - 1.5, 350.708160 + 1.5)]
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 580, 2, 970.65 / 60.0, 2, 579, 0.00000, 0, 0.00, DissociationType.HCD, -2, 2)]
        [TestCase(@"Ms2Align_TopFDv1.6.2_ms2.msalign", 1669, 2, 2726.74 / 60.0, 0, 1665, 1239.30464, 8, 1329780.59, DissociationType.HCD, 1239.30464 - 2, 1239.30464 + 2)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1366, 2, 832.915128799 / 60.0, 24, 1365, 1005.95751953125, 5, 619272.56, DissociationType.HCD, 1005.1555541753769, 1006.3555542230606)]
        [TestCase(@"Ms2Align_IsoDec1.0.0_ms2.msalign", 1398, 2, 862.421989215 / 60, 22, 1395, 973.6828002929688, 6, 122901736.0, DissociationType.HCD, 973.0828002691269, 974.2828003168106)]
        public void TestMsAlign_HeaderComponents_StaticLoading(string path, int oneBasedScanNumber, int msnOrder,
           double retentionTime,
           int peakCount, int? oneBasePrecursorScanNumber, double? precursorMz, int? precursorCharge,
           double? precursorIntensity, DissociationType? dissociationType,
           double? precursorIsolationMzStart, double? precursorIsolationMzEnd)
        {
            var file = MsAlignTestFiles[path];
            Assert.That(file.GetSourceFile().Id, Is.EqualTo(msnOrder.ToString()));
            Assert.That(file.GetSourceFile().FileName, Is.EqualTo(path.GetPeriodTolerantFilenameWithoutExtension().Replace("_ms1", "").Replace("_ms2", "")));
            Assert.That(file.GetSourceFile().NativeIdFormat, Is.EqualTo("no nativeID format"));
            Assert.That(file.GetSourceFile().MassSpectrometerFileFormat, Is.EqualTo($"ms{msnOrder}.msalign format"));

            var scanToTest = file.GetOneBasedScan(oneBasedScanNumber);

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
        public void GetOneBasedScanFromDynamicConnection_ConnectionNotInitialized_ThrowsException()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign");

            var msAlign = new Ms1Align(path);

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

        [Test]
        public static void TestAlternativeConstructors_Ms1Align()
        {
            var msAlign = MsAlignTestFiles.First().Value;

            var align = new Ms1Align(msAlign.Scans, msAlign.SourceFile);
            Assert.That(align.Scans.Length, Is.EqualTo(msAlign.Scans.Length));
            for (var index = 0; index < msAlign.Scans.Length; index++)
            {
                var ogScan = msAlign.Scans[index];
                var otherScan = align.Scans[index];
                Assert.That(otherScan, Is.EqualTo(ogScan));
            }

            align = new Ms1Align(msAlign.NumSpectra, msAlign.SourceFile);
            Assert.That(align.Scans.Length, Is.EqualTo(msAlign.NumSpectra));
        }

        [Test]
        public static void TestAlternativeConstructors_Ms2Align()
        {
            var msAlign = MsAlignTestFiles.Last().Value;

            var align = new Ms2Align(msAlign.Scans, msAlign.SourceFile);
            Assert.That(align.Scans.Length, Is.EqualTo(msAlign.Scans.Length));
            for (var index = 0; index < msAlign.Scans.Length; index++)
            {
                var ogScan = msAlign.Scans[index];
                var otherScan = align.Scans[index];
                Assert.That(otherScan, Is.EqualTo(ogScan));
            }

            align = new Ms2Align(msAlign.NumSpectra, msAlign.SourceFile);
            Assert.That(align.Scans.Length, Is.EqualTo(msAlign.NumSpectra));
        }

        [Test]
        public void ParseHeader_HeaderNotPresent()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign");

            var msAlign = new Ms1Align(path);
            msAlign.LoadAllStaticData();

            Assert.That(msAlign.FileName, Is.Null);
            Assert.That(msAlign.FaimsData, Is.Null);
            Assert.That(msAlign.FaimsVoltage, Is.Null);
            Assert.That(msAlign.DissociationType, Is.Null);
            Assert.That(msAlign.Ms1ScanCount, Is.Null);
            Assert.That(msAlign.Ms2ScanCount, Is.Null);
            Assert.That(msAlign.SpectralDataType, Is.Null);
            Assert.That(msAlign.MaxAssumedChargeState, Is.Null);
            Assert.That(msAlign.MaxAssumedMonoisotopicMass, Is.Null);
            Assert.That(msAlign.PeakErrorTolerance, Is.Null);
            Assert.That(msAlign.Ms1SnRRatio, Is.Null);
            Assert.That(msAlign.Ms2SnRRatio, Is.Null);
            Assert.That(msAlign.MaxThreadsToUse, Is.Null);
            Assert.That(msAlign.PrecursorWindowSize, Is.Null);
            Assert.That(msAlign.UseEnvCnnModel, Is.Null);
            Assert.That(msAlign.MissMs1Spectra, Is.Null);
            Assert.That(msAlign.SoftwareVersion, Is.Null);
            Assert.That(msAlign.Software, Is.EqualTo(Software.Unspecified));
            Assert.That(msAlign.UseMsDeconvScore, Is.Null);
        }

        [Test]
        public void ParseHeader_HeaderPresent()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", "Ms1Align_TopFDv1.6.2_ms1.msalign");

            Ms1Align msAlign = new Ms1Align(path);
            msAlign.LoadAllStaticData();

            Assert.That(msAlign.FileName, Is.Null);
            Assert.That(msAlign.FaimsData, Is.Null);
            Assert.That(msAlign.FaimsVoltage, Is.Null);
            Assert.That(msAlign.DissociationType, Is.EqualTo(DissociationType.Autodetect));
            Assert.That(msAlign.Ms1ScanCount, Is.EqualTo(2906));
            Assert.That(msAlign.Ms2ScanCount, Is.EqualTo(2201));
            Assert.That(msAlign.SpectralDataType, Is.EqualTo("Centroid"));
            Assert.That(msAlign.MaxAssumedChargeState, Is.EqualTo(60));
            Assert.That(msAlign.MaxAssumedMonoisotopicMass, Is.EqualTo(70000));
            Assert.That(msAlign.PeakErrorTolerance, Is.EqualTo("0.02 m/z"));
            Assert.That(msAlign.Ms1SnRRatio, Is.EqualTo(3));
            Assert.That(msAlign.Ms2SnRRatio, Is.EqualTo(1));
            Assert.That(msAlign.MaxThreadsToUse, Is.EqualTo(4));
            Assert.That(msAlign.PrecursorWindowSize, Is.EqualTo(4));
            Assert.That(msAlign.UseEnvCnnModel, Is.EqualTo(false));
            Assert.That(msAlign.MissMs1Spectra, Is.EqualTo(false));
            Assert.That(msAlign.SoftwareVersion, Is.EqualTo("1.6.2"));
            Assert.That(msAlign.Software, Is.EqualTo(Software.TopFD));
            Assert.That(msAlign.UseMsDeconvScore, Is.Null);
        }

        [Test]
        public void ParseHeader_Ms2Align_IsoDec1_0_0()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", "Ms2Align_IsoDec1.0.0_ms2.msalign");

            var msAlign = new Ms2Align(path);
            msAlign.LoadAllStaticData();

            Assert.That(msAlign.FileName, Is.EqualTo("20170105_L_MaD_ColC4_Ecoli20161108-10_01.msalign"));
            Assert.That(msAlign.FaimsData, Is.EqualTo(false));
            Assert.That(msAlign.FaimsVoltage, Is.Null);
            Assert.That(msAlign.Ms1ScanCount, Is.EqualTo(4430));
            Assert.That(msAlign.Ms2ScanCount, Is.EqualTo(6111));
            Assert.That(msAlign.SpectralDataType, Is.EqualTo("centroid"));
            Assert.That(msAlign.PeakErrorTolerance, Is.EqualTo("0.01"));
            Assert.That(msAlign.Ms1SnRRatio, Is.EqualTo(3));
            Assert.That(msAlign.Ms2SnRRatio, Is.EqualTo(1));
            Assert.That(msAlign.MaxThreadsToUse, Is.EqualTo(1));
            Assert.That(msAlign.PrecursorWindowSize, Is.EqualTo(3));
            Assert.That(msAlign.DissociationType, Is.EqualTo(DissociationType.Autodetect));
            Assert.That(msAlign.UseMsDeconvScore, Is.EqualTo(false));
            Assert.That(msAlign.MissMs1Spectra, Is.EqualTo(true));
            Assert.That(msAlign.SoftwareVersion, Is.EqualTo("1.0.0"));
            Assert.That(msAlign.Software, Is.EqualTo(Software.IsoDec));
        }


        [Test]
        public void NeutralMassSpectrum_FromMsAlign()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", "Ms2Align_IsoDec1.0.0_ms2.msalign");

            var msAlign = new Ms2Align(path);
            msAlign.LoadAllStaticData();

            int firstScanNumber = 1366;
            int lastScanNumber = 1398;
            int peakCount = 24;
            int ms2ScansPresent = 22;

            string scanIonText =
                "5020.7255859375\t432147.2\t5\t1\r\n5022.744140625\t73592.586\t3\t1\r\n5021.72802734375\t279469.47\t4\t1\r\n4648.5078125\t18255.004\t3\t1\r\n1673.0023193359375\t2019.978\t2\t1\r\n5004.7060546875\t58158.74\t4\t1\r\n4520.5087890625\t8129.425\t3\t1\r\n5005.7216796875\t14691.099\t3\t1\r\n559.2601318359375\t820.32465\t1\t1\r\n4804.61083984375\t7133.3687\t4\t1\r\n4977.7041015625\t36881.18\t4\t1\r\n4121.30615234375\t5057.6055\t3\t1\r\n4463.48291015625\t6635.46\t3\t1\r\n4577.5322265625\t6578.043\t3\t1\r\n4771.61474609375\t6805.8135\t3\t1\r\n4978.716796875\t8181.794\t3\t1\r\n3674.077880859375\t7278.764\t3\t1\r\n4963.6884765625\t22462.277\t4\t1\r\n4994.73974609375\t6223.2183\t4\t1\r\n2578.2978515625\t1839.8989\t2\t1\r\n4492.3974609375\t6426.9375\t3\t1\r\n4634.5498046875\t5513.514\t3\t1\r\n4963.689453125\t3275.681\t3\t1\r\n1675.367431640625\t2546.7227\t1\t1";
            (double Mono, double Intensity, int Charge, double Mz)[] ions = scanIonText.Split("\r\n").Select(p =>
            {
                var items = p.Split("\t");
                return (double.Parse(items[0]), double.Parse(items[1]), int.Parse(items[2]), double.Parse(items[0]).ToMz(int.Parse(items[2])));
            }).OrderBy(p => p.Item1).ToArray();

            Assert.That(msAlign.Scans.Length, Is.EqualTo(ms2ScansPresent));

            // Check that isolation range was read correctly
            var firstScan = msAlign.GetOneBasedScan(firstScanNumber);
            Assert.That(firstScan.OneBasedPrecursorScanNumber, Is.EqualTo(1365));
            Assert.That(firstScan.IsolationRange.Minimum, Is.EqualTo(1005.1555541753769));
            Assert.That(firstScan.IsolationRange.Maximum, Is.EqualTo(1006.3555542230606));

            // test the neutral mass spectrum
            Assert.That(firstScan.MassSpectrum, Is.TypeOf<NeutralMassSpectrum>());
            NeutralMassSpectrum firstSpectrum = (NeutralMassSpectrum)firstScan.MassSpectrum;
            Assert.That(firstSpectrum.Size, Is.EqualTo(peakCount));
            Assert.That(firstSpectrum.XArray.Length, Is.EqualTo(peakCount));
            Assert.That(firstSpectrum.YArray.Length, Is.EqualTo(peakCount));
            Assert.That(firstSpectrum.Charges.Length, Is.EqualTo(peakCount));
            Assert.That(firstSpectrum.FirstX, Is.EqualTo(ions.Min(p => p.Mz)));
            Assert.That(firstSpectrum.LastX, Is.EqualTo(ions.Max(p => p.Mz)));

            for (int i = 0; i < firstSpectrum.Size; i++)
            {
                double mass = firstSpectrum.XArray[i];
                double intensity = firstSpectrum.YArray[i];
                int charge = firstSpectrum.Charges[i];

                double expectedMass = ions[i].Mono;
                double expectedIntensity = ions[i].Intensity;
                int expectedCharge = ions[i].Charge;

                Assert.That(mass, Is.EqualTo(expectedMass));
                Assert.That(intensity, Is.EqualTo(expectedIntensity));
                Assert.That(charge, Is.EqualTo(expectedCharge));
            }
        }

        [Test]
        public void CombineMsAlign_ShouldCombineMs1AndMs2Align()
        {
            Ms1Align ms1Align = (Ms1Align)MsAlignTestFiles["Ms1Align_IsoDec1.0.0_ms1.msalign"];
            Ms2Align ms2Align = (Ms2Align)MsAlignTestFiles["Ms2Align_IsoDec1.0.0_ms2.msalign"];
            var expectedScansCount = ms1Align.Scans.Length + ms2Align.Scans.Length;
            
            var combinedMsAlign = MsAlign.CombineMsAlign(ms1Align, ms2Align);
            
            Assert.That(expectedScansCount, Is.EqualTo(combinedMsAlign.Scans.Length));
            
            int previousScanNumber = -1;
            foreach(var scan in combinedMsAlign)
            {
                Assert.That(scan.OneBasedScanNumber, Is.GreaterThanOrEqualTo(previousScanNumber));
                previousScanNumber = scan.OneBasedScanNumber;
            }
        }
    }
}
