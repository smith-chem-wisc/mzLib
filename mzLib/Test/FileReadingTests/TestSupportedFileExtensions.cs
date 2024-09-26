using System;
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
        [Test]
        [TestCase("DataFiles/sliced_ethcd.raw", SupportedFileType.ThermoRaw)]
        [TestCase("DataFiles/SmallCalibratibleYeast.mzml", SupportedFileType.MzML)]
        [TestCase("DataFiles/tester.mgf", SupportedFileType.Mgf)]
        [TestCase("DataFiles/tester.d", SupportedFileType.BrukerD)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms2Feature_FlashDeconvjurkat_td_rep1_fract2_ms2.feature", SupportedFileType.Ms2Feature)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDMs1Feature_jurkat_td_rep1_fract2_ms1.feature", SupportedFileType.Ms1Feature)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\TopFDmzrt_jurkat_td_rep1_fract2_frac.mzrt.csv", SupportedFileType.Mzrt_TopFd)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Tsv_FlashDeconvjurkat_td_rep1_fract2_ms1.tsv", SupportedFileType.Ms1Tsv_FlashDeconv)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Tsv_FlashDeconvOpenMs3.0.0.tsv", SupportedFileType.Tsv_FlashDeconv)]
        [TestCase(@"FileReadingTests\SearchResults\ExcelEditedPeptide.psmtsv", SupportedFileType.psmtsv)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv", SupportedFileType.ToppicPrsm)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv", SupportedFileType.ToppicProteoform)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv", SupportedFileType.ToppicProteoformSingle)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv", SupportedFileType.ToppicPrsmSingle)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv", SupportedFileType.MsFraggerPsm)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1individual_peptide.tsv", SupportedFileType.MsFraggerPeptide)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1individual_protein.tsv", SupportedFileType.MsFraggerProtein)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1combined_peptide.tsv", SupportedFileType.MsFraggerPeptide)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1combined_protein.tsv", SupportedFileType.MsFraggerProtein)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv", SupportedFileType.FlashLFQQuantifiedPeak)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_TargetResults_IcTarget.tsv", SupportedFileType.MsPathFinderTTargets)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_DecoyResults_IcDecoy.tsv", SupportedFileType.MsPathFinderTDecoys)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_AllResults_IcTda.tsv", SupportedFileType.MsPathFinderTAllResults)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\crux.txt", SupportedFileType.CruxResult)]
        public static void TestSupportedFileTypeExtensions(string filePath, SupportedFileType expectedType)
        {
            var supportedType = filePath.ParseFileType();
            Assert.That(supportedType, Is.EqualTo(expectedType));

            var extension = expectedType.GetFileExtension();
            Assert.That(filePath.EndsWith(extension, StringComparison.InvariantCultureIgnoreCase));
        }

        [Test]
        public static void EnsureAllExtensionsAreUnique()
        {
            var allExtensions = Enum.GetValues<SupportedFileType>();
            Assert.That(allExtensions.Length, Is.EqualTo(allExtensions.Distinct().Count()));
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
