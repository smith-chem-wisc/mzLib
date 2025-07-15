using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Readers;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    internal class TestQuantifiedPeak
    {
        internal static string TestDirectory;
        internal static string TestFilePath;

        [OneTimeSetUp]
        public void SetUp()
        {
            TestDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ReadingWritingTests");
            TestFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv");
            Directory.CreateDirectory(TestDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(TestDirectory, true);
        }

        [Test]
        public static void TestFileLoadsAndCountCorrect()
        {
            QuantifiedPeakFile file = new QuantifiedPeakFile(TestFilePath);
            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));

            file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));
        }

        [Test]
        public static void TestFileFirstAndLastAreCorrect()
        {
            QuantifiedPeakFile file = new QuantifiedPeakFile(TestFilePath);
            var first = file.First();

            Assert.That(first.FileName, Is.EqualTo("20100721_Velos1_TaGe_SA_A549_06-calib-averaged"));
            Assert.That(first.BaseSequence, Is.EqualTo("DICNDVLSLLEK"));
            Assert.That(first.FullSequence, Is.EqualTo("DIC[Common Fixed:Carbamidomethyl on C]NDVLSLLEK"));
            Assert.That(first.ProteinGroup, Is.EqualTo("P63104"));
            Assert.That(first.PeptideMonoisotopicMass, Is.EqualTo(1417.712281191));
            Assert.That(first.MS2RetentionTime, Is.EqualTo(188.04623));
            Assert.That(first.PrecursorCharge, Is.EqualTo(3));
            Assert.That(first.TheoreticalMZ, Is.EqualTo(473.578036863879));
            Assert.That(first.PeakIntensity, Is.EqualTo(61339740.974156484));
            Assert.That(first.PeakRTStart, Is.EqualTo(187.71813));
            Assert.That(first.PeakRTApex, Is.EqualTo(188.27129333333335));
            Assert.That(first.PeakRTEnd, Is.EqualTo(195.134625));
            Assert.That(first.PeakMz, Is.EqualTo(709.8649279277056));
            Assert.That(first.PeakCharge, Is.EqualTo(2));
            Assert.That(first.NumChargeStatesObserved, Is.EqualTo(3));
            Assert.That(first.PeakDetectionType, Is.EqualTo("MSMS"));
            Assert.That(first.MBRScore, Is.EqualTo(0));
            Assert.That(first.PSMsMapped, Is.EqualTo(2));
            Assert.That(first.BaseSequencesMapped, Is.EqualTo(1));
            Assert.That(first.FullSequencesMapped, Is.EqualTo(1));
            Assert.That(first.PeakSplitValleyRT, Is.EqualTo(0));
            Assert.That(first.PeakApexMassError, Is.EqualTo(2.1314131880687888));

            var last = file.Last();
            Assert.That(last.FileName, Is.EqualTo("20101230_Velos1_TaGe_SA_Jurkat6-calib-averaged"));
            Assert.That(last.BaseSequence, Is.EqualTo("QDLEAQIRGLREEVEK"));
            Assert.That(last.FullSequence, Is.EqualTo("QDLEAQIRGLREEVEK"));
            Assert.That(last.ProteinGroup, Is.EqualTo("A1A5D9"));
            Assert.That(last.PeptideMonoisotopicMass, Is.EqualTo(1912.001403494));
            Assert.That(last.MS2RetentionTime, Is.EqualTo(71.64922));
            Assert.That(last.PrecursorCharge, Is.EqualTo(4));
            Assert.That(last.TheoreticalMZ, Is.EqualTo(479.007627340379));
            Assert.That(last.PeakIntensity, Is.EqualTo(0));
            Assert.That(last.PeakRTStart, Is.Null);
            Assert.That(last.PeakRTApex, Is.Null);
            Assert.That(last.PeakRTEnd, Is.Null);
            Assert.That(last.PeakMz, Is.Null);
            Assert.That(last.PeakCharge, Is.Null);
            Assert.That(last.NumChargeStatesObserved, Is.EqualTo(0));
            Assert.That(last.PeakDetectionType, Is.EqualTo("MSMS"));
            Assert.That(last.MBRScore, Is.EqualTo(0));
            Assert.That(last.PSMsMapped, Is.EqualTo(1));
            Assert.That(last.BaseSequencesMapped, Is.EqualTo(1));
            Assert.That(last.FullSequencesMapped, Is.EqualTo(1));
            Assert.That(last.PeakSplitValleyRT, Is.EqualTo(0));
            Assert.That(last.PeakApexMassError, Is.EqualTo(double.NaN));
        }

        [Test]
        public static void TestFileReadWrite_WithoutExtensionInPath()
        {
            var file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            var testOutputPath = Path.Combine(TestDirectory, "TestOutput");

            file.WriteResults(testOutputPath);
            var newPath = testOutputPath + file.FileType.GetFileExtension();
            Assert.That(File.Exists(newPath));

            var writtenFile = new QuantifiedPeakFile(newPath);
            Assert.That(file.Count(), Is.EqualTo(writtenFile.Count()));

            for (int i = 0; i < file.Count(); i++)
            {
                var originalPeak = file.Results[i];
                var writtenPeak = writtenFile.Results[i];
                Assert.That(originalPeak.FileName, Is.EqualTo(writtenPeak.FileName));
                Assert.That(originalPeak.BaseSequence, Is.EqualTo(writtenPeak.BaseSequence));
                Assert.That(originalPeak.FullSequence, Is.EqualTo(writtenPeak.FullSequence));
                Assert.That(originalPeak.ProteinGroup, Is.EqualTo(writtenPeak.ProteinGroup));
                Assert.That(originalPeak.PeptideMonoisotopicMass, Is.EqualTo(writtenPeak.PeptideMonoisotopicMass).Within(0.0000001));
                Assert.That(originalPeak.MS2RetentionTime, Is.EqualTo(writtenPeak.MS2RetentionTime).Within(0.0000001));
                Assert.That(originalPeak.PrecursorCharge, Is.EqualTo(writtenPeak.PrecursorCharge));
                Assert.That(originalPeak.TheoreticalMZ, Is.EqualTo(writtenPeak.TheoreticalMZ).Within(0.0000001));
                Assert.That(originalPeak.PeakIntensity, Is.EqualTo(writtenPeak.PeakIntensity).Within(0.0000001));
                Assert.That(originalPeak.PeakRTStart, Is.EqualTo(writtenPeak.PeakRTStart).Within(0.0000001));
                Assert.That(originalPeak.PeakRTApex, Is.EqualTo(writtenPeak.PeakRTApex).Within(0.0000001));
                Assert.That(originalPeak.PeakRTEnd, Is.EqualTo(writtenPeak.PeakRTEnd).Within(0.0000001));
                Assert.That(originalPeak.PeakMz, Is.EqualTo(writtenPeak.PeakMz).Within(0.0000001));
                Assert.That(originalPeak.PeakCharge, Is.EqualTo(writtenPeak.PeakCharge));
                Assert.That(originalPeak.NumChargeStatesObserved, Is.EqualTo(writtenPeak.NumChargeStatesObserved));
                Assert.That(originalPeak.PeakDetectionType, Is.EqualTo(writtenPeak.PeakDetectionType));
                Assert.That(originalPeak.MBRScore, Is.EqualTo(writtenPeak.MBRScore).Within(0.0000001));
                Assert.That(originalPeak.PSMsMapped, Is.EqualTo(writtenPeak.PSMsMapped));
                Assert.That(originalPeak.BaseSequencesMapped, Is.EqualTo(writtenPeak.BaseSequencesMapped));
                Assert.That(originalPeak.FullSequencesMapped, Is.EqualTo(writtenPeak.FullSequencesMapped));
                Assert.That(originalPeak.PeakSplitValleyRT, Is.EqualTo(writtenPeak.PeakSplitValleyRT));
                Assert.That(originalPeak.PeakApexMassError, Is.EqualTo(writtenPeak.PeakApexMassError).Within(0.0000001));
            }
        }

        [Test]
        public static void TestFileReadWrite_WithExtensionInPath()
        {
            var file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            var testOutputPath = Path.Combine(TestDirectory, "TestOutput_QuantifiedPeaks.tsv");

            file.WriteResults(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            var writtenFile = new QuantifiedPeakFile(testOutputPath);
            Assert.That(file.Count(), Is.EqualTo(writtenFile.Count()));

            for (int i = 0; i < file.Count(); i++)
            {
                var originalPeak = file.Results[i];
                var writtenPeak = writtenFile.Results[i];
                Assert.That(originalPeak.FileName, Is.EqualTo(writtenPeak.FileName));
                Assert.That(originalPeak.BaseSequence, Is.EqualTo(writtenPeak.BaseSequence));
                Assert.That(originalPeak.FullSequence, Is.EqualTo(writtenPeak.FullSequence));
                Assert.That(originalPeak.ProteinGroup, Is.EqualTo(writtenPeak.ProteinGroup));
                Assert.That(originalPeak.PeptideMonoisotopicMass, Is.EqualTo(writtenPeak.PeptideMonoisotopicMass).Within(0.0000001));
                Assert.That(originalPeak.MS2RetentionTime, Is.EqualTo(writtenPeak.MS2RetentionTime).Within(0.0000001));
                Assert.That(originalPeak.PrecursorCharge, Is.EqualTo(writtenPeak.PrecursorCharge));
                Assert.That(originalPeak.TheoreticalMZ, Is.EqualTo(writtenPeak.TheoreticalMZ).Within(0.0000001));
                Assert.That(originalPeak.PeakIntensity, Is.EqualTo(writtenPeak.PeakIntensity).Within(0.0000001));
                Assert.That(originalPeak.PeakRTStart, Is.EqualTo(writtenPeak.PeakRTStart).Within(0.0000001));
                Assert.That(originalPeak.PeakRTApex, Is.EqualTo(writtenPeak.PeakRTApex).Within(0.0000001));
                Assert.That(originalPeak.PeakRTEnd, Is.EqualTo(writtenPeak.PeakRTEnd).Within(0.0000001));
                Assert.That(originalPeak.PeakMz, Is.EqualTo(writtenPeak.PeakMz).Within(0.0000001));
                Assert.That(originalPeak.PeakCharge, Is.EqualTo(writtenPeak.PeakCharge));
                Assert.That(originalPeak.NumChargeStatesObserved, Is.EqualTo(writtenPeak.NumChargeStatesObserved));
                Assert.That(originalPeak.PeakDetectionType, Is.EqualTo(writtenPeak.PeakDetectionType));
                Assert.That(originalPeak.MBRScore, Is.EqualTo(writtenPeak.MBRScore).Within(0.0000001));
                Assert.That(originalPeak.PSMsMapped, Is.EqualTo(writtenPeak.PSMsMapped));
                Assert.That(originalPeak.BaseSequencesMapped, Is.EqualTo(writtenPeak.BaseSequencesMapped));
                Assert.That(originalPeak.FullSequencesMapped, Is.EqualTo(writtenPeak.FullSequencesMapped));
                Assert.That(originalPeak.PeakSplitValleyRT, Is.EqualTo(writtenPeak.PeakSplitValleyRT));
                Assert.That(originalPeak.PeakApexMassError, Is.EqualTo(writtenPeak.PeakApexMassError).Within(0.0000001));
            }
        }
    }
}
