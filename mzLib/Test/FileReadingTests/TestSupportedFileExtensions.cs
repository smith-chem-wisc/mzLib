using System;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSupportedFileExtensions
    {
        [Test]
        [TestCase("DataFiles/sliced_ethcd.raw", SupportedFileType.ThermoRaw)]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", SupportedFileType.MzML)]
        [TestCase("DataFiles/tester.mgf", SupportedFileType.Mgf)]
        [TestCase("DataFiles/tester.d", SupportedFileType.BrukerD)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvjurkat_td_rep1_fract2_ms2.feature", SupportedFileType.Ms2Feature)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDMs1Feature_jurkat_td_rep1_fract2_ms1.feature", SupportedFileType.Ms1Feature)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDmzrt_jurkat_td_rep1_fract2_frac.mzrt.csv", SupportedFileType.Mzrt_TopFd)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Tsv_FlashDeconvjurkat_td_rep1_fract2_ms1.tsv", SupportedFileType.Ms1Tsv_FlashDeconv)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Tsv_FlashDeconvjurkat_td_rep1_fract2.tsv", SupportedFileType.Tsv_FlashDeconv)]
        [TestCase(@"FileReadingTests\SearchResults\ExcelEditedPeptide.psmtsv", SupportedFileType.psmtsv)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv", SupportedFileType.ToppicPrsm)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv", SupportedFileType.ToppicProteoform)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv", SupportedFileType.ToppicProteoformSingle)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv", SupportedFileType.ToppicPrsmSingle)]
        public static void TestSupportedFileTypeExtensions(string filePath, SupportedFileType expectedType)
        {
            var supportedType = filePath.ParseFileType();
            Assert.That(supportedType, Is.EqualTo(expectedType));

            var extension = expectedType.GetFileExtension();
            Assert.That(filePath.EndsWith(extension, StringComparison.InvariantCultureIgnoreCase));
        }

        [Test]
        public static void TestSupportedFileTypeExtension_Errors()
        {
            string badTest = "badFile.taco";
            Exception e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message,
                Is.EqualTo($"File type not supported"));

            badTest = "badTaco.feature";
            e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message,
                Is.EqualTo($"Feature file type not supported"));

            badTest = "badTaco.psm.csv";
            e = Assert.Throws<MzLibException>(() => badTest.ParseFileType());
            Assert.That(e?.Message,
                Is.EqualTo($"Csv file type not supported"));

            // assure all values of enum have a file extension in the swithc method
            foreach (var value in Enum.GetValues<SupportedFileType>())
            {
                _ = value.GetFileExtension();
            }
        }

    }
}
