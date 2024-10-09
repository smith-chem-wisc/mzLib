using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
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
                { "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign")) },
                { "Ms1Align_IsoDec1.0.0_ms1.msalign", new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms1Align_IsoDec1.0.0_ms1.msalign")) },
                { "Ms1Align_TopFDv1.6.2_ms1.msalign", new Ms1Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms1Align_TopFDv1.6.2_ms1.msalign")) },
                { "Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign", new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign")) },
                { "Ms2Align_IsoDec1.0.0_ms2.msalign", new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms2Align_IsoDec1.0.0_ms2.msalign")) },
                { "Ms2Align_TopFDv1.6.2_ms2.msalign", new Ms2Align(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes", "Ms2Align_TopFDv1.6.2_ms2.msalign")) }
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
        public void TestMsAlignResultsLoadsAndCountCorrect_GenericMsDataFileReader(string path, int recordCount, int expectedMsAlignNumber)
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
        public void TestMsAlignResultsLoadsAndCountCorrect_GenericFileReader(string path, int recordCount, int expectedMsAlignNumber)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes", path);

            MsDataFile file = FileReader.ReadFile<MsDataFileToResultFileAdapter>(filePath);
            Assert.That(file.Count(), Is.EqualTo(recordCount)); 
            Assert.That(file.Scans.Length, Is.EqualTo(recordCount)); 
            Assert.That(file.NumSpectra, Is.EqualTo(recordCount)); 
            Assert.That(file.SourceFile.MassSpectrometerFileFormat, Contains.Substring("align")); 
            Assert.That(file.SourceFile.MassSpectrometerFileFormat, Contains.Substring(expectedMsAlignNumber.ToString())); 
        }

        [Test]
        [TestCase(@"Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign", 1, 1, 0.0035, 189, null, null, null, null, null, null, null, null)]
        public void TestMsAlign_AllSpectrumHeaderComponents(string path, int oneBasedScanNumber, int msnOrder, double retentionTime,
            int peakCount, int? oneBasePrecursorScanNumber, double? precursorMz, int? precursorCharge, double? precursorMass, 
            double? precursorIntensity, DissociationType? dissociationType,
           double? precursorIsolationMzStart, double? precursorIsolationMzEnd)
        {
            var file = MsAlignTestFiles[path];

            var scanToTest = file.GetOneBasedScan(oneBasedScanNumber);

        }
    }
}
