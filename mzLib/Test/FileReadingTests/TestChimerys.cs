using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    public class TestChimerys
    {
        private static string directoryPath;

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ChimerysReadingWritingTests");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_psms.tsv", 6)]
        public void ChimerysPsm_LoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPsmFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_peptides.tsv", 5)]
        public void ChimerysPeptides_LoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPeptideFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_modified_peptides.tsv", 6)]
        public void ChimerysModifiedPeptides_LoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysModifiedPeptideFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_precursors.tsv", 10)]
        public void ChimerysPrecursors_LoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPrecursorFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_protein_groups.tsv", 4)]
        public void ChimerysProteinGroups_LoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysProteinGroupFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        public void ChimerysPsm_FirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_psms.tsv");
            ChimerysPsmFile file = new(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.PsmId, Is.EqualTo(0));
            Assert.That(first.BaseSequence, Is.EqualTo("TAGINVR"));
            Assert.That(first.ModifiedSequence, Is.EqualTo("TAGINVR"));
            Assert.That(first.MissedCleavages, Is.EqualTo(0));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(729.4133372));
            Assert.That(first.PrecursorCharge, Is.EqualTo(2));
            Assert.That(first.Length, Is.EqualTo(7));
            Assert.That(first.Mz, Is.EqualTo(365.7139451));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(2649));
            Assert.That(first.RetentionTime, Is.EqualTo(19.29169167));
            Assert.That(first.RetentionTimePrediction, Is.EqualTo(19.23022461));
            Assert.That(first.Coeff, Is.EqualTo(799382.1204));
            Assert.That(first.SpectralAngle, Is.EqualTo(0.799476716));
            Assert.That(first.Ctp, Is.EqualTo(0.773829709));
            Assert.That(first.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.QValue, Is.EqualTo(0.000102037));
            Assert.That(first.SearchEngineScore, Is.EqualTo(4.390117569));
            Assert.That(first.Pep, Is.EqualTo(0.000132497));
            Assert.That(first.IsDecoy, Is.EqualTo(false));
            Assert.That(first.PrecursorId, Is.EqualTo(4292402));
            Assert.That(first.ModifiedPeptideId, Is.EqualTo(42924));
            Assert.That(first.PeptideId, Is.EqualTo(263426));
            Assert.That(first.PositionInProteinIds, Is.EqualTo(new[] { 490 }));
            Assert.That(first.ProteinIds, Is.EqualTo(new long[] { 2483 }));
            Assert.That(first.LocalizationSequence, Is.EqualTo(""));
            Assert.That(first.LocalizationScore, Is.EqualTo(0));
            Assert.That(first.ProteinSites, Is.EqualTo(""));

            Assert.That(last.PsmId, Is.EqualTo(5));
            Assert.That(last.BaseSequence, Is.EqualTo("FAVLHGEAPR"));
            Assert.That(last.ModifiedSequence, Is.EqualTo("FAVLHGEAPR"));
            Assert.That(last.MissedCleavages, Is.EqualTo(0));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(1095.582528));
            Assert.That(last.PrecursorCharge, Is.EqualTo(3));
            Assert.That(last.Length, Is.EqualTo(10));
            Assert.That(last.Mz, Is.EqualTo(366.2014524));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(6100));
            Assert.That(last.RetentionTime, Is.EqualTo(32.520635));
            Assert.That(last.RetentionTimePrediction, Is.EqualTo(33.37529755));
            Assert.That(last.Coeff, Is.EqualTo(4129995.817));
            Assert.That(last.SpectralAngle, Is.EqualTo(0.878042855));
            Assert.That(last.Ctp, Is.EqualTo(0.824411588));
            Assert.That(last.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.QValue, Is.EqualTo(0.000102037));
            Assert.That(last.SearchEngineScore, Is.EqualTo(4.562185625));
            Assert.That(last.Pep, Is.EqualTo(0.000132497));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
            Assert.That(last.PrecursorId, Is.EqualTo(538443403));
            Assert.That(last.ModifiedPeptideId, Is.EqualTo(5384434));
            Assert.That(last.PeptideId, Is.EqualTo(1281602));
            Assert.That(last.PositionInProteinIds, Is.EqualTo(new[] { 577 }));
            Assert.That(last.ProteinIds, Is.EqualTo(new long[] { 11801 }));
            Assert.That(last.LocalizationSequence, Is.EqualTo(""));
            Assert.That(last.LocalizationScore, Is.EqualTo(0));
            Assert.That(last.ProteinSites, Is.EqualTo(""));
        }

        [Test]
        public void ChimerysPeptides_FirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_peptides.tsv");
            ChimerysPeptideFile file = new(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.PeptideId, Is.EqualTo(24));
            Assert.That(first.BaseSequence, Is.EqualTo("AKWLTPK"));
            Assert.That(first.Length, Is.EqualTo(7));
            Assert.That(first.IsAmbiguous, Is.EqualTo(false));
            Assert.That(first.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.QValue, Is.EqualTo(0.647762358));
            Assert.That(first.SearchEngineScore, Is.EqualTo(-7.763616337));
            Assert.That(first.Pep, Is.EqualTo(1));
            Assert.That(first.GlobalQValue, Is.EqualTo(0.647762337));
            Assert.That(first.GlobalSearchEngineScore, Is.EqualTo(-7.763616337));
            Assert.That(first.GlobalPep, Is.EqualTo(1));
            Assert.That(first.IsDecoy, Is.EqualTo(false));
            Assert.That(first.IsIdentifiedByMbr, Is.EqualTo(false));
            Assert.That(first.PsmIds, Is.EqualTo(new long[] { 99003 }));
            Assert.That(first.PrecursorIds, Is.EqualTo(new long[] { 172127102 }));
            Assert.That(first.ModifiedPeptideIds, Is.EqualTo(new long[] { 1721271 }));
            Assert.That(first.PositionInProteinIds, Is.EqualTo(new int[] { 122 }));
            Assert.That(first.ProteinIds, Is.EqualTo(new long[] { 0 }));
            Assert.That(first.CountPsms, Is.EqualTo(0));
            Assert.That(first.CountPrecursors, Is.EqualTo(0));
            Assert.That(first.CountModifiedPeptides, Is.EqualTo(0));
            Assert.That(first.Quantification, Is.EqualTo(0));

            Assert.That(last.PeptideId, Is.EqualTo(380));
            Assert.That(last.BaseSequence, Is.EqualTo("HQIKHAGENLTTK"));
            Assert.That(last.Length, Is.EqualTo(13));
            Assert.That(last.IsAmbiguous, Is.EqualTo(false));
            Assert.That(last.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.QValue, Is.EqualTo(0.486905187));
            Assert.That(last.SearchEngineScore, Is.EqualTo(-4.524063016));
            Assert.That(last.Pep, Is.EqualTo(1));
            Assert.That(last.GlobalQValue, Is.EqualTo(0.486905195));
            Assert.That(last.GlobalSearchEngineScore, Is.EqualTo(-4.524063016));
            Assert.That(last.GlobalPep, Is.EqualTo(1));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
            Assert.That(last.IsIdentifiedByMbr, Is.EqualTo(false));
            Assert.That(last.PsmIds, Is.EqualTo(new long[] { 77470 }));
            Assert.That(last.PrecursorIds, Is.EqualTo(new long[] { 2092141203 }));
            Assert.That(last.ModifiedPeptideIds, Is.EqualTo(new long[] { 20921412 }));
            Assert.That(last.PositionInProteinIds, Is.EqualTo(new int[] { 779, 483 }));
            Assert.That(last.ProteinIds, Is.EqualTo(new long[] { 1, 38917 }));
            Assert.That(last.CountPsms, Is.EqualTo(0));
            Assert.That(last.CountPrecursors, Is.EqualTo(0));
            Assert.That(last.CountModifiedPeptides, Is.EqualTo(0));
            Assert.That(last.Quantification, Is.EqualTo(0));
        }

        [Test]
        public void ChimerysModifiedPeptides_FirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_modified_peptides.tsv");
            ChimerysModifiedPeptideFile file = new(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.ModifiedPeptideId, Is.EqualTo(71));
            Assert.That(first.BaseSequence, Is.EqualTo("GIGAAGI"));
            Assert.That(first.ModifiedSequence, Is.EqualTo("GIGAAGI"));
            Assert.That(first.MissedCleavages, Is.EqualTo(0));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(557.3173115449999));
            Assert.That(first.Length, Is.EqualTo(7));
            Assert.That(first.MaxSpectralAngle, Is.EqualTo(0));
            Assert.That(first.MaxCtp, Is.EqualTo(0));
            Assert.That(first.IsAmbiguous, Is.EqualTo(false));
            Assert.That(first.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.QValue, Is.EqualTo(0.45052990317344666));
            Assert.That(first.SearchEngineScore, Is.EqualTo(-4.15937637));
            Assert.That(first.Pep, Is.EqualTo(1));
            Assert.That(first.GlobalQValue, Is.EqualTo(0.450529916));
            Assert.That(first.GlobalSearchEngineScore, Is.EqualTo(-4.15937637));
            Assert.That(first.GlobalPep, Is.EqualTo(0.997465117));
            Assert.That(first.IsDecoy, Is.EqualTo(true));
            Assert.That(first.IsIdentifiedByMbr, Is.EqualTo(false));
            Assert.That(first.PeptideId, Is.EqualTo(4847746));
            Assert.That(first.PsmIds, Is.EqualTo(new long[] { 170009, 170610 }));
            Assert.That(first.PrecursorIds, Is.EqualTo(new long[] { 7101 }));
            Assert.That(first.PositionInProteinIds, Is.EqualTo(new[] { 23 }));
            Assert.That(first.ProteinIds, Is.EqualTo(new[] {148728}));
            Assert.That(first.LocalizationSequence, Is.EqualTo(""));
            Assert.That(first.LocalizationScore, Is.EqualTo(0));
            Assert.That(first.ProteinSites, Is.EqualTo(""));
            Assert.That(first.CountPsms, Is.EqualTo(0));
            Assert.That(first.CountPrecursors, Is.EqualTo(0));
            Assert.That(first.Quantification, Is.EqualTo(0));

            Assert.That(last.ModifiedPeptideId, Is.EqualTo(408));
            Assert.That(last.BaseSequence, Is.EqualTo("AGPGLGK"));
            Assert.That(last.ModifiedSequence, Is.EqualTo("AGPGLGK"));
            Assert.That(last.MissedCleavages, Is.EqualTo(0));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(598.3438607));
            Assert.That(last.Length, Is.EqualTo(7));
            Assert.That(last.MaxSpectralAngle, Is.EqualTo(0));
            Assert.That(last.MaxCtp, Is.EqualTo(0));
            Assert.That(last.IsAmbiguous, Is.EqualTo(false));
            Assert.That(last.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.QValue, Is.EqualTo(0.22785022854804993));
            Assert.That(last.SearchEngineScore, Is.EqualTo(-2.489934666));
            Assert.That(last.Pep, Is.EqualTo(0.9958520477749483));
            Assert.That(last.GlobalQValue, Is.EqualTo(0.22785022936330668));
            Assert.That(last.GlobalSearchEngineScore, Is.EqualTo(-2.489934666));
            Assert.That(last.GlobalPep, Is.EqualTo(0.981258036));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
            Assert.That(last.IsIdentifiedByMbr, Is.EqualTo(false));
            Assert.That(last.PeptideId, Is.EqualTo(2402820));
            Assert.That(last.PsmIds, Is.EqualTo(new long[] { 124045 }));
            Assert.That(last.PrecursorIds, Is.EqualTo(new long[] { 40802 }));
            Assert.That(last.PositionInProteinIds, Is.EqualTo(new[] { 6 }));
            Assert.That(last.ProteinIds, Is.EqualTo(new long[] { 68644 }));
            Assert.That(last.LocalizationSequence, Is.EqualTo(""));
            Assert.That(last.LocalizationScore, Is.EqualTo(0));
            Assert.That(last.ProteinSites, Is.EqualTo(""));
            Assert.That(last.CountPsms, Is.EqualTo(0));
            Assert.That(last.CountPrecursors, Is.EqualTo(0));
            Assert.That(last.Quantification, Is.EqualTo(0));
        }

        [Test]
        public void ChimerysProteinGroups_FirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_protein_groups.tsv");
            ChimerysProteinGroupFile file = new(filePath);
            var first = file.First();
            var last = file.Last();

          

            Assert.That(first.ProteinIds,Is.EqualTo(new long[]
                {
                    22025, 36547, 39902, 39987, 64766, 68570, 72743, 76011, 1792, 34904, 50610, 53511, 57131, 64236
                }));
            Assert.That(first.ProteinGroupId, Is.EqualTo(0));
            Assert.That(first.FastaHeaders, Is.EqualTo(new string[]
            {
                ">tr|A0A286YF78|A0A286YF78_HUMAN D-3-phosphoglycerate dehydrogenase (Fragment) OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFB2|A0A286YFB2_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YF34|A0A286YF34_HUMAN D-3-phosphoglycerate dehydrogenase (Fragment) OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFC8|A0A286YFC8_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YER3|A0A286YER3_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFF3|A0A286YFF3_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFK5|A0A286YFK5_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFE1|A0A286YFE1_HUMAN D-3-phosphoglycerate dehydrogenase (Fragment) OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">sp|O43175|SERA_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=4",
                ">tr|A0A2C9F2M7|A0A2C9F2M7_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YF22|A0A286YF22_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFL2|A0A286YFL2_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFA2|A0A286YFA2_HUMAN D-3-phosphoglycerate dehydrogenase OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1",
                ">tr|A0A286YFM8|A0A286YFM8_HUMAN D-3-phosphoglycerate dehydrogenase (Fragment) OS=Homo sapiens OX=9606 GN=PHGDH PE=1 SV=1"
            }));
            Assert.That(first.GeneNames, Is.EqualTo(new string[]
            {
                "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH", "PHGDH",
                "PHGDH", "PHGDH", "PHGDH"
            }));
            Assert.That(first.ProteinIdentifiers, Is.EqualTo(new string[]
            {
                "A0A286YF78", "A0A286YFB2", "A0A286YF34", "A0A286YFC8", "A0A286YER3", "A0A286YFF3", "A0A286YFK5",
                "A0A286YFE1", "O43175", "A0A2C9F2M7", "A0A286YF22", "A0A286YFL2", "A0A286YFA2", "A0A286YFM8"
            }));
            Assert.That(first.TaxonomyIds, Is.EqualTo(new string[]
            {
                "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606", "9606",
                "9606"
            }));
            Assert.That(first.Organisms, Is.EqualTo(new string[]
            {
                "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens",
                "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens",
                "Homo sapiens", "Homo sapiens"
            }));
            Assert.That(first.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.GlobalQValue, Is.EqualTo(0.000248139));
            Assert.That(first.GlobalSearchEngineScore, Is.EqualTo(5.846220359).Within(0.0000001));
            Assert.That(first.IsDecoy, Is.EqualTo(false));
            Assert.That(first.PsmIds, Is.EqualTo(new long[]
            {
                50616, 50824, 175471, 198611, 145911, 146235, 73324, 134618, 195346, 159043, 159191, 96157, 96258, 165195,
                165439, 165445, 21829, 34908, 18031, 210843, 210920, 211287, 75472, 60151, 192743, 202821, 29210, 36081
            }));
            Assert.That(first.PrecursorIds, Is.EqualTo(new long[]
            {
                1462667102, 1071882902, 622825102, 1626340602, 486401402, 309281302, 125013502, 40845602, 546063102,
                3505532303, 3505532304, 1205382202, 2058473302, 2058473303, 289754001, 289754002, 663080502, 5528193603,
                3541538504
            }));
            Assert.That(first.ModifiedPeptideIds, Is.EqualTo(new long[]
            {
                14626671, 10718829, 6228251, 16263406, 4864014, 3092813, 1250135, 408456, 5460631, 35055323, 12053822,
                20584733, 2897540, 6630805, 55281936, 35415385
            }));
            Assert.That(first.PeptideIds, Is.EqualTo(new long[]
            {
                189057, 189062, 189070, 189074, 189078, 189085, 189088, 189094, 189099, 189104, 189116, 189118, 189128,
                2141192, 2222334, 2381407
            }));
            Assert.That(first.CountPsms, Is.EqualTo(23));
            Assert.That(first.CountPrecursors, Is.EqualTo(15));
            Assert.That(first.CountModifiedPeptides, Is.EqualTo(12));
            Assert.That(first.CountPeptides, Is.EqualTo(12));
            Assert.That(first.Quantification, Is.EqualTo(93556932.23));

            Assert.That(last.ProteinIds,Is.EqualTo(new long[]
            {
                14056, 14059
            }));
            Assert.That(last.ProteinGroupId, Is.EqualTo(3));
            Assert.That(last.FastaHeaders, Is.EqualTo(new string[]
            {
                ">sp|Q00266|METK1_HUMAN S-adenosylmethionine synthase isoform type-1 OS=Homo sapiens OX=9606 GN=MAT1A PE=1 SV=2", 
                ">sp|P31153|METK2_HUMAN S-adenosylmethionine synthase isoform type-2 OS=Homo sapiens OX=9606 GN=MAT2A PE=1 SV=1"
            }));
            Assert.That(last.GeneNames, Is.EqualTo(new string[]
            {
                "MAT1A", "MAT2A"
            }));
            Assert.That(last.ProteinIdentifiers, Is.EqualTo(new string[]
            {
                "Q00266", "P31153"
            }));
            Assert.That(last.TaxonomyIds, Is.EqualTo(new string[]
            {
                "9606", "9606"
            }));
            Assert.That(last.Organisms, Is.EqualTo(new string[]
            {
                "Homo sapiens", "Homo sapiens"
            }));
            Assert.That(last.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.GlobalQValue, Is.EqualTo(0.000248139));
            Assert.That(last.GlobalSearchEngineScore, Is.EqualTo(5.718183408).Within(0.0000001));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
            Assert.That(last.PsmIds, Is.EqualTo(new long[]
            {
                127769, 172454, 156196, 62433, 21095, 21143, 21153, 21198, 21617, 21714, 109933, 27317, 171647, 111912,
                164198, 22121, 22131, 22752, 34239, 134912, 52156, 52524, 153695, 190864, 191846, 126611, 166675
            }));
            Assert.That(last.PrecursorIds, Is.EqualTo(new long[]
            {
                329200403, 1035831402, 1035831403, 93056502, 1250717602, 2608437103, 5410179603, 5410179604, 1313813402,
                1313813403, 1269284202, 3626002, 4996369504, 1493677702, 54551502, 273856102, 273856103, 558137302
            }));
            Assert.That(last.ModifiedPeptideIds, Is.EqualTo(new long[]
            {
                3292004, 10358314, 930565, 12507176, 26084371, 54101796, 13138134, 12692842, 36260, 49963695, 14936777,
                545515, 2738561, 5581373
            }));
            Assert.That(last.PeptideIds, Is.EqualTo(new long[]
            {
                1513364, 1513415, 1513553, 1513556, 1513567, 1513569, 1513576, 1513580, 1513582, 1513591, 1513594, 1513604,
                1513612, 1513615
            }));
            Assert.That(last.CountPsms, Is.EqualTo(24));
            Assert.That(last.CountPrecursors, Is.EqualTo(15));
            Assert.That(last.CountModifiedPeptides, Is.EqualTo(11));
            Assert.That(last.CountPeptides, Is.EqualTo(11));
            Assert.That(last.Quantification, Is.EqualTo(80669868.84));
        }

        [Test]
        public void ChimerysPrecursors_FirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_precursors.tsv");
            ChimerysPrecursorFile file = new(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.PrecursorId, Is.EqualTo(7101));
            Assert.That(first.BaseSequence, Is.EqualTo("GIGAAGI"));
            Assert.That(first.ModifiedSequence, Is.EqualTo("GIGAAGI"));
            Assert.That(first.MissedCleavages, Is.EqualTo(0));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(557.3173115));
            Assert.That(first.PrecursorCharge, Is.EqualTo(1));
            Assert.That(first.Length, Is.EqualTo(7));
            Assert.That(first.Mz, Is.EqualTo(558.324588));
            Assert.That(first.MinRetentionTime, Is.EqualTo(0));
            Assert.That(first.MaxRetentionTime, Is.EqualTo(0));
            Assert.That(first.MaxSpectralAngle, Is.EqualTo(0));
            Assert.That(first.MaxCtp, Is.EqualTo(0));
            Assert.That(first.IsAmbiguous, Is.False);
            Assert.That(first.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(first.QValue, Is.EqualTo(0.426193237));
            Assert.That(first.SearchEngineScore, Is.EqualTo(-4.15937637));
            Assert.That(first.Pep, Is.EqualTo(1));
            Assert.That(first.GlobalQValue, Is.EqualTo(0.426193242));
            Assert.That(first.GlobalSearchEngineScore, Is.EqualTo(-4.15937637));
            Assert.That(first.GlobalPep, Is.EqualTo(0.991515389));
            Assert.That(first.IsDecoy, Is.True);
            Assert.That(first.IsIdentifiedByMbr, Is.False);
            Assert.That(first.ModifiedPeptideId, Is.EqualTo(71));
            Assert.That(first.PeptideId, Is.EqualTo(4847746));
            Assert.That(first.PsmIds, Is.EqualTo(new long[] { 170009, 170610 }));
            Assert.That(first.PositionInProteinIds, Is.EqualTo(new[] {23}));
            Assert.That(first.ProteinIds, Is.EqualTo(new long[] { 148728 }));
            Assert.That(first.LocalizationSequence, Is.EqualTo(""));
            Assert.That(first.LocalizationScore, Is.EqualTo(0));
            Assert.That(first.ProteinSites, Is.EqualTo(""));
            Assert.That(first.CountPsms, Is.EqualTo(0));
            Assert.That(first.Quantification, Is.EqualTo(0));

            Assert.That(last.PrecursorId, Is.EqualTo(56502));
            Assert.That(last.BaseSequence, Is.EqualTo("GLAGSAK"));
            Assert.That(last.ModifiedSequence, Is.EqualTo("GLAGSAK"));
            Assert.That(last.MissedCleavages, Is.EqualTo(0));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(602.3387753));
            Assert.That(last.PrecursorCharge, Is.EqualTo(2));
            Assert.That(last.Length, Is.EqualTo(7));
            Assert.That(last.Mz, Is.EqualTo(302.1766641));
            Assert.That(last.MinRetentionTime, Is.EqualTo(0));
            Assert.That(last.MaxRetentionTime, Is.EqualTo(0));
            Assert.That(last.MaxSpectralAngle, Is.EqualTo(0));
            Assert.That(last.MaxCtp, Is.EqualTo(0));
            Assert.That(last.IsAmbiguous, Is.False);
            Assert.That(last.RawFileName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.SampleName, Is.EqualTo("20100604_Velos1_TaGe_SA_A549_1.raw"));
            Assert.That(last.QValue, Is.EqualTo(0.383429915));
            Assert.That(last.SearchEngineScore, Is.EqualTo(-3.797142909));
            Assert.That(last.Pep, Is.EqualTo(1));
            Assert.That(last.GlobalQValue, Is.EqualTo(0.383429904));
            Assert.That(last.GlobalSearchEngineScore, Is.EqualTo(-3.797142909));
            Assert.That(last.GlobalPep, Is.EqualTo(0.991515389));
            Assert.That(last.IsDecoy, Is.True);
            Assert.That(last.IsIdentifiedByMbr, Is.False);
            Assert.That(last.ModifiedPeptideId, Is.EqualTo(565));
            Assert.That(last.PeptideId, Is.EqualTo(3969410));
            Assert.That(last.PsmIds, Is.EqualTo(new long[] { 124224 }));
            Assert.That(last.PositionInProteinIds, Is.EqualTo(new[] {26}));
            Assert.That(last.ProteinIds, Is.EqualTo(new long[] { 93434 }));
            Assert.That(last.LocalizationSequence, Is.EqualTo(""));
            Assert.That(last.LocalizationScore, Is.EqualTo(0));
            Assert.That(last.ProteinSites, Is.EqualTo(""));
            Assert.That(last.CountPsms, Is.EqualTo(0));
            Assert.That(last.Quantification, Is.EqualTo(0));
        }
    }
}
