using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSupportedFileExtensions
    {
        private static IEnumerable<TestCaseData> SupportedFileTypeTestCases()
        {
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvOpenMs3.0.0_ms2.feature", SupportedFileType.Ms2Feature);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature", SupportedFileType.Ms1Feature);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\mzrt_TopFDv1.6.2.mzrt.csv", SupportedFileType.TopFDMzrt);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\Ms1Tsv_FlashDeconvOpenMs3.0.0_ms1.tsv", SupportedFileType.Ms1Tsv_FlashDeconv);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\Tsv_FlashDeconvOpenMs3.0.0.tsv", SupportedFileType.Tsv_FlashDeconv);
            yield return new TestCaseData(@"FileReadingTests\SearchResults\ExcelEditedPeptide.psmtsv", SupportedFileType.psmtsv);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv", SupportedFileType.ToppicPrsm);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv", SupportedFileType.ToppicProteoform);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv", SupportedFileType.ToppicProteoformSingle);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv", SupportedFileType.ToppicPrsmSingle);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv", SupportedFileType.MsFraggerPsm);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1individual_peptide.tsv", SupportedFileType.MsFraggerPeptide);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1individual_protein.tsv", SupportedFileType.MsFraggerProtein);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1combined_peptide.tsv", SupportedFileType.MsFraggerPeptide);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1combined_protein.tsv", SupportedFileType.MsFraggerProtein);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv", SupportedFileType.FlashLFQQuantifiedPeak);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_TargetResults_IcTarget.tsv", SupportedFileType.MsPathFinderTTargets);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_DecoyResults_IcDecoy.tsv", SupportedFileType.MsPathFinderTDecoys);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_AllResults_IcTda.tsv", SupportedFileType.MsPathFinderTAllResults);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\crux.txt", SupportedFileType.CruxResult);
            yield return new TestCaseData(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults\experiment_annotation.tsv", SupportedFileType.ExperimentAnnotation);
            yield return new TestCaseData(@"FileReadingTests\SearchResults\XL_Intralinks.tsv", SupportedFileType.psmtsv);
            yield return new TestCaseData(@"DataFiles\sliced_ethcd.raw", SupportedFileType.ThermoRaw);
            yield return new TestCaseData(@"DataFiles\SmallCalibratibleYeast.mzml", SupportedFileType.MzML);
            yield return new TestCaseData(@"DataFiles\tester.mgf", SupportedFileType.Mgf);
            yield return new TestCaseData(@"DataFiles\centroid_1x_MS1_4x_autoMS2.d", SupportedFileType.BrukerD);
            yield return new TestCaseData(@"DataFiles\timsTOF_snippet.d", SupportedFileType.BrukerTimsTof);
        }

        private static IEnumerable<SupportedFileType> EnumTestCases() => Enum.GetValues<SupportedFileType>();

        [Test, TestCaseSource(nameof(SupportedFileTypeTestCases))]
        [TestCase(@"FileReadingTests\ExternalFileTypes\XL_Intralinks_MIons.tsv", SupportedFileType.psmtsv)] // Read in intralink XL file even if it has a different extension
        public static void TestSupportedFileTypeExtensions(string filePath, SupportedFileType expectedType)
        {
            var supportedType = filePath.ParseFileType();
            Assert.That(supportedType, Is.EqualTo(expectedType));

            // These files are different ways of naming a psm file and will thus have different extensions when written by metamorpheus. 
            // The extension is not guaranteed to be the same as the file type.
            string[] skipTheseChecks = ["XL_Intralinks.tsv", "XL_Intralinks_MIons.tsv"];
            if (skipTheseChecks.Contains(Path.GetFileName(filePath)))
                return;

            var extension = expectedType.GetFileExtension();
            Assert.That(filePath.EndsWith(extension, StringComparison.InvariantCultureIgnoreCase));
        }

        [Test, TestCaseSource(nameof(SupportedFileTypeTestCases))]
        public static void TestIResultFileReaderWorks(string filePath, SupportedFileType expectedType)
        {
            IResultFile resultFile = FileReader.ReadResultFile(filePath);
            Assert.That(resultFile.FileType, Is.EqualTo(expectedType));
            Type resultFileClass = expectedType.GetResultFileType();
            var convertedFile = Convert.ChangeType(resultFile, resultFileClass);
            Assert.That(convertedFile, Is.Not.Null);
        }

        [Test]
        [TestCase("DatabaseTests/bad4.fasta", typeof(MzLibException))]
        [TestCase("DataFiles/xyzl.psmtsv", typeof(FileNotFoundException))]
        public static void TestIResultFileReaderThrowsException(string filePath, Type expectedExceptionType)
        {
            Assert.Throws(expectedExceptionType, () =>
            {
                FileReader.ReadResultFile(filePath);
            });
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv", SupportedFileType.MsFraggerPsm)]
        public static void TestIQuantifiableResultFileReaderWorks(string filePath, SupportedFileType expectedType)
        {
            IQuantifiableResultFile resultFile = FileReader.ReadQuantifiableResultFile(filePath);
            Assert.That(resultFile.FileType, Is.EqualTo(expectedType));
            Type resultFileClass = expectedType.GetResultFileType();
            var convertedFile = Convert.ChangeType(resultFile, resultFileClass);
            Assert.That(convertedFile, Is.Not.Null);
        }

        [Test]
        [TestCase("DataFiles/sliced_ethcd.raw", typeof(MzLibException))]
        [TestCase("DataFiles/xyzl.psmtsv", typeof(FileNotFoundException))]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature", typeof(MzLibException))]
        public static void TestIQuantifiableResultFileReaderThrowsException(string filePath, Type expectedExceptionType)
        {
            Assert.Throws(expectedExceptionType, () =>
            {
                IResultFile resultFile = FileReader.ReadQuantifiableResultFile(filePath);
            });
        }

        [Test]
        public static void EnsureAllExtensionsAreUnique()
        {
            var allExtensions = Enum.GetValues<SupportedFileType>();
            Assert.That(allExtensions.Length, Is.EqualTo(allExtensions.Distinct().Count()));
        }

        [TestCaseSource(nameof(EnumTestCases))]
        public static void EnsureAllEnumExtensionsAreImplemented(SupportedFileType enumType)
        {
            // Ensure extensions are implemented for all enum values
            string extension = enumType.GetFileExtension();
            Assert.That(extension, Is.Not.Null.Or.Empty, $"Extension for {enumType} is not implemented.");

            // Ensure the enumType determination is implemented for all values 
            Type type = enumType.GetResultFileType();
            Assert.That(type, Is.Not.Null, $"Type for {enumType} is not implemented.");
        }

        [Test]
        public static void TestSupportedFileTypeExtension_Errors()
        {

            string badTest = "badFile.taco";
            Exception e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message, Is.EqualTo($"File type not supported"));

            badTest = "badTaco.feature";
            e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message, Is.EqualTo($"Feature file type not supported"));

            badTest = "badTaco.psm.csv";
            e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message, Is.EqualTo($"Csv file type not supported"));

            badTest = Path.Combine(TestContext.CurrentContext.TestDirectory, "DoubleProtease.tsv");
            e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message, Is.EqualTo($"Tsv file type not supported"));

            var emptyFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "emptyFile.tsv");
            File.Create(emptyFile).Close();
            e = Assert.Throws<MzLibException>(() => emptyFile.ParseFileType());
            Assert.That(e?.Message, Is.EqualTo($"Tsv file is empty"));
            File.Delete(emptyFile);

            // assure all values of enum have a file extension in the swithc method
            foreach (var value in Enum.GetValues<SupportedFileType>())
            {
                _ = value.GetFileExtension();
            }
        }

        [Test]
        public static void TestGetFileExtension_Errors()
        {
            Exception e = Assert.Throws<MzLibException>(() => ((SupportedFileType)100).GetFileExtension());
            Assert.That(e?.Message,
                Is.EqualTo($"File type not supported"));
        }
    }
}
