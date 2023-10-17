using System;
using System.Collections.Generic;
using System.Data.Entity.Core.Objects;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Newtonsoft.Json;
using NUnit.Framework;
using pepXML.Generated;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    public class TestToppicResultFiles
    {
        private static string directoryPath;

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ReadingWritingTests");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv", 4)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv", 37)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv", 7)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv", 10)]
        public void TestToppicProteoformAndPrsmLoadAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));

            file = FileReader.ReadFile<ToppicSearchResultFile>(path);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv")]
        public void TestToppicHeaderParsing(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            file.LoadResults();

            Assert.That(file.ProteinDatabasePath, Is.EqualTo("D:/Databases/Human_uniprotkb_proteome_UP000005640_AND_revi_2023_09_29.fasta"));
            Assert.That(file.SpectrumFilePath, Is.EqualTo("B:/Users/Nic/SharedWithMe/targetBias/Centroided\\HumanFasta_TopPicDecoys_Small_ms2.msalign"));
            Assert.That(file.NumberOfCombinedSpectra, Is.EqualTo(1));
            Assert.That(file.FragmentationMethod, Is.EqualTo(DissociationType.HCD));
            Assert.That(file.SearchType, Is.EqualTo("TARGET+DECOY"));

            CollectionAssert.AreEquivalent(file.FixedModifications,
                new List<string> { "Carbamidomethylation 57.0215 C" });
            CollectionAssert.AreEquivalent(file.AllowedNTerminalForms,
                new List<string> { "NONE", "NME", "NME_ACETYLATION", "M_ACETYLATION" });

            Assert.That(file.NumberOfMaxUnexpectedModifications, Is.EqualTo(1));
            Assert.That(file.MaximumMassShift, Is.EqualTo(500.0));
            Assert.That(file.MinimumMassShift, Is.EqualTo(-500.0));
            Assert.That(file.SpectrumLevelCutOffType, Is.EqualTo("FDR"));
            Assert.That(file.SpectrumLevelCutOffValue, Is.EqualTo(100.0));
            Assert.That(file.ProteoformLevelCutOffType, Is.EqualTo("FDR"));
            Assert.That(file.ProteoformLevelCutOffValue, Is.EqualTo(100.0));
            Assert.That(file.PrecursorErrorTolerance, Is.EqualTo(10));
            Assert.That(file.PrsmClusterErrorTolerance, Is.EqualTo(1.2));
            Assert.That(file.UseToppicFeatureFile);
            Assert.That(file.EValueComputation, Is.EqualTo("Generating function"));
            Assert.That(!file.LocalizationWithMIScore);
            Assert.That(file.ThreadNumber, Is.EqualTo(4));
            Assert.That(file.ExecutableFileDirectory, Is.EqualTo(@"C:\Users\Nic\Downloads\toppic-windows-1.6.2\toppic-windows-1.6.2"));
            Assert.That(file.StartTime, Is.EqualTo(new DateTime(2023, 9, 29, 10, 14, 12)));
            Assert.That(file.EndTime, Is.EqualTo(new DateTime(2023, 9, 30, 3, 0, 8)));
            Assert.That(file.Version, Is.EqualTo("1.6.2"));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv")]
        public void TestToppicProteoformFileResultReading(string path)
        {
            string filePath = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            Assert.That(file.FileType,
                path.Contains("single")
                    ? Is.EqualTo(SupportedFileType.ToppicProteoformSingle)
                    : Is.EqualTo(SupportedFileType.ToppicProteoform));
            file.LoadResults();

            var first = file.First();
            Assert.That(first.FilePath, Is.EqualTo("B:/Users/Nic/SharedWithMe/targetBias/Centroided/05-26-17_B7A_yeast_td_fract5_rep1_ms2.msalign"));
            Assert.That(first.PrsmID, Is.EqualTo(2));
            Assert.That(first.ScanNum, Is.EqualTo(10));
            Assert.That(first.DissociationType, Is.EqualTo(DissociationType.HCD));
            Assert.That(first.Scans, Is.EqualTo(354));
            Assert.That(first.RetentionTime, Is.EqualTo(465.09));
            Assert.That(first.PeakCount, Is.EqualTo(16));
            Assert.That(first.PrecursorCharge, Is.EqualTo(1));
            Assert.That(first.PrecursorMass, Is.EqualTo(535.1598));
            Assert.That(first.AdjustedPrecursorMass, Is.EqualTo(534.1675).Within(0.001));
            Assert.That(first.ProteoformId, Is.EqualTo(1346).Within(0.001));
            Assert.That(first.FeatureIntensity, Is.EqualTo(1.038E+07).Within(0.001));
            Assert.That(first.FeatureScore, Is.EqualTo(-1000));
            Assert.That(first.FeatureApexTime, Is.EqualTo(548.52).Within(0.001));
            Assert.That(first.ProteinHitsCount, Is.EqualTo(1));
            Assert.That(first.ProteinAccession, Is.EqualTo("DECOY_sp|Q9NRR8|C42S1_HUMAN"));
            Assert.That(first.ProteinDescription, Is.EqualTo("CDC42 small effector protein 1 OS=Homo sapiens OX=9606 GN=CDC42SE1 PE=1 SV=1"));
            Assert.That(first.FirstResidue, Is.EqualTo(2));
            Assert.That(first.LastResidue, Is.EqualTo(9));
            Assert.That(first.SpecialAminoAcids, Is.EqualTo(""));
            Assert.That(first.BaseSequence, Is.EqualTo("SAECTKRF"));
            Assert.That(first.FullSequence, Is.EqualTo("M.SAE(C)[Carbamidomethylation](TKRF)[-463.2977].D"));
            Assert.That(first.FullSequenceMass, Is.EqualTo(534.1675).Within(0.001));
            Assert.That(first.ProteinNTerminalForm, Is.EqualTo("NME"));
            Assert.That(first.FixedPTMs, Is.EqualTo("Carbamidomethylation:[4]"));
            Assert.That(first.UnexpectedModificationsCount, Is.EqualTo(1));
            Assert.That(first.UnexpectedModifications, Is.EqualTo("-463.2977:[5-8]"));
            Assert.That(first.VariableModificationsCount, Is.EqualTo(0));
            Assert.That(first.VariableModifications, Is.EqualTo(string.Empty));
            Assert.That(first.MIScore, Is.EqualTo(null));
            Assert.That(first.MatchedPeaksCount, Is.EqualTo(1));
            Assert.That(first.MatchedFragmentIonsCount, Is.EqualTo(2));
            Assert.That(first.EValue, Is.EqualTo(1E+300).Within(0.001));
            Assert.That(first.QValueSpectrumLevel, Is.EqualTo(1));
            Assert.That(first.QValueProteoformLevel, Is.EqualTo(1));
        }

        [Test]
        public void TestToppicProteoformsFileAlternativeResults()
        {
            string path = @"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv";
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.ToppicProteoform));
            file.LoadResults();

            var containsAlternatives = file.First(p => p.PrsmID == 37);
            Assert.That(containsAlternatives.AlternativeIdentifications.Count, Is.EqualTo(3));

            var alternative = containsAlternatives.AlternativeIdentifications[0];
            Assert.That(alternative.PrsmId, Is.EqualTo(37));
            Assert.That(alternative.Accession, Is.EqualTo("sp|Q14162|SREC_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("Scavenger receptor class F member 1 OS=Homo sapiens OX=9606 GN=SCARF1 PE=1 SV=3"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(281));
            Assert.That(alternative.LastResidue, Is.EqualTo(285));

            alternative = containsAlternatives.AlternativeIdentifications[1];
            Assert.That(alternative.PrsmId, Is.EqualTo(37));
            Assert.That(alternative.Accession, Is.EqualTo("DECOY_sp|Q2MKA7|RSPO1_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("R-spondin-1 OS=Homo sapiens OX=9606 GN=RSPO1 PE=1 SV=1"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(175));
            Assert.That(alternative.LastResidue, Is.EqualTo(179));

            alternative = containsAlternatives.AlternativeIdentifications[2];
            Assert.That(alternative.PrsmId, Is.EqualTo(37));
            Assert.That(alternative.Accession, Is.EqualTo("DECOY_sp|Q58EX7|PKHG4_HUMAN"));
            Assert.That(alternative.ProteinDescription, 
                Is.EqualTo("Puratrophin-1 OS=Homo sapiens OX=9606 GN=PLEKHG4 PE=1 SV=1"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(1174));
            Assert.That(alternative.LastResidue, Is.EqualTo(1178));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv")]
        public void TestToppicPrsmFileResultReading(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            Assert.That(file.FileType,
                path.Contains("single")
                    ? Is.EqualTo(SupportedFileType.ToppicPrsmSingle)
                    : Is.EqualTo(SupportedFileType.ToppicPrsm));

            file.LoadResults();

            var first = file.First();
            Assert.That(first.FilePath, Is.EqualTo("B:/Users/Nic/SharedWithMe/targetBias/Centroided/05-26-17_B7A_yeast_td_fract5_rep1_ms2.msalign"));
            Assert.That(first.PrsmID, Is.EqualTo(0));
            Assert.That(first.ScanNum, Is.EqualTo(7));
            Assert.That(first.DissociationType, Is.EqualTo(DissociationType.HCD));
            Assert.That(first.Scans, Is.EqualTo(259));
            Assert.That(first.RetentionTime, Is.EqualTo(340.12));
            Assert.That(first.PeakCount, Is.EqualTo(14));
            Assert.That(first.PrecursorCharge, Is.EqualTo(1));
            Assert.That(first.PrecursorMass, Is.EqualTo(576.1219).Within(0.001));
            Assert.That(first.AdjustedPrecursorMass, Is.EqualTo(575.1295).Within(0.001));
            Assert.That(first.ProteoformId, Is.EqualTo(1029).Within(0.001));
            Assert.That(first.FeatureIntensity, Is.EqualTo(3.604E+06).Within(0.001));
            Assert.That(first.FeatureScore, Is.EqualTo(-1000));
            Assert.That(first.FeatureApexTime, Is.EqualTo(536.51).Within(0.001));
            Assert.That(first.ProteinHitsCount, Is.EqualTo(1));
            Assert.That(first.ProteinAccession, Is.EqualTo("DECOY_sp|Q8N5H7|SH2D3_HUMAN"));
            Assert.That(first.ProteinDescription, Is.EqualTo("SH2 domain-containing protein 3C OS=Homo sapiens OX=9606 GN=SH2D3C PE=1 SV=1"));
            Assert.That(first.FirstResidue, Is.EqualTo(719));
            Assert.That(first.LastResidue, Is.EqualTo(724));
            Assert.That(first.SpecialAminoAcids, Is.EqualTo(""));
            Assert.That(first.BaseSequence, Is.EqualTo("GYASDS"));
            Assert.That(first.FullSequence, Is.EqualTo("R.GYA(S)[-23.0940]DS.P"));
            Assert.That(first.FullSequenceMass, Is.EqualTo(575.1295).Within(0.001));
            Assert.That(first.ProteinNTerminalForm, Is.EqualTo("NONE"));
            Assert.That(first.FixedPTMs, Is.EqualTo(""));
            Assert.That(first.UnexpectedModificationsCount, Is.EqualTo(1));
            Assert.That(first.UnexpectedModifications, Is.EqualTo("-23.0940:[4]"));
            Assert.That(first.VariableModificationsCount, Is.EqualTo(0));
            Assert.That(first.VariableModifications, Is.EqualTo(string.Empty));
            Assert.That(first.MIScore, Is.EqualTo(null));
            Assert.That(first.MatchedPeaksCount, Is.EqualTo(2));
            Assert.That(first.MatchedFragmentIonsCount, Is.EqualTo(3));
            Assert.That(first.EValue, Is.EqualTo(1E+300).Within(0.001));
            Assert.That(first.QValueSpectrumLevel, Is.EqualTo(1));
            Assert.That(first.QValueProteoformLevel, Is.EqualTo(1));
        }

        [Test]
        public void TestToppicPrsmsFileAlternativeResults()
        {
            string path = @"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv";
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ToppicSearchResultFile file = new ToppicSearchResultFile(filePath);
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.ToppicPrsm));
            file.LoadResults();

            var containsAlternatives = file.First(p => p.PrsmID == 1);
            Assert.That(containsAlternatives.AlternativeIdentifications.Count, Is.EqualTo(4));

            var alternative = containsAlternatives.AlternativeIdentifications[0];
            Assert.That(alternative.PrsmId, Is.EqualTo(1));
            Assert.That(alternative.Accession, Is.EqualTo("DECOY_sp|P21781|FGF7_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("Fibroblast growth factor 7 OS=Homo sapiens OX=9606 GN=FGF7 PE=1 SV=1"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(101));
            Assert.That(alternative.LastResidue, Is.EqualTo(105));

            alternative = containsAlternatives.AlternativeIdentifications[1];
            Assert.That(alternative.PrsmId, Is.EqualTo(1));
            Assert.That(alternative.Accession, Is.EqualTo("sp|Q13423|NNTM_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("NAD(P) transhydrogenase, mitochondrial OS=Homo sapiens OX=9606 GN=NNT PE=1 SV=3"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(369));
            Assert.That(alternative.LastResidue, Is.EqualTo(373));

            alternative = containsAlternatives.AlternativeIdentifications[2];
            Assert.That(alternative.PrsmId, Is.EqualTo(1));
            Assert.That(alternative.Accession, Is.EqualTo("DECOY_sp|Q6ZT12|UBR3_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("E3 ubiquitin-protein ligase UBR3 OS=Homo sapiens OX=9606 GN=UBR3 PE=2 SV=2"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(1805));
            Assert.That(alternative.LastResidue, Is.EqualTo(1809));

            alternative = containsAlternatives.AlternativeIdentifications[3];
            Assert.That(alternative.PrsmId, Is.EqualTo(1));
            Assert.That(alternative.Accession, Is.EqualTo("sp|Q96RW7|HMCN1_HUMAN"));
            Assert.That(alternative.ProteinDescription,
                Is.EqualTo("Hemicentin-1 OS=Homo sapiens OX=9606 GN=HMCN1 PE=1 SV=2"));
            Assert.That(alternative.FirstResidue, Is.EqualTo(179));
            Assert.That(alternative.LastResidue, Is.EqualTo(183));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofrom_TopPICv1.6.2_proteoform.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicProteofromSingle_TopPICv1.6.2_proteoform_single.tsv")]
        public static void TestTopicProteoformsReadWrite(string path)
        {
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var testOutputPath = path.Contains("single")
                ? Path.Combine(directoryPath, "toppic_proteoform_single.tsv")
                : Path.Combine(directoryPath, "toppic_proteoform.tsv");

            ToppicSearchResultFile file = new ToppicSearchResultFile(filepath);
            file.LoadResults();
            file.WriteResults(testOutputPath);
            var writtenFile = FileReader.ReadFile<ToppicSearchResultFile>(testOutputPath);
            writtenFile.LoadResults();
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < file.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(file.Results[i]);
                var written = JsonConvert.SerializeObject(writtenFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            var originalLines = File.ReadAllLines(filepath);
            var writtenLines = File.ReadAllLines(testOutputPath);
            Assert.That(writtenLines.Length, Is.EqualTo(originalLines.Length));
            int paramCount = 0;
            for (int i = 0; i < originalLines.Length; i++)
            {
                if (originalLines[i].Contains("********"))
                    paramCount++;
                if (paramCount >= 2)
                    break;
                Assert.That(writtenLines[i], Is.EqualTo(originalLines[i]));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "toppic");
            file.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsm_TopPICv1.6.2_prsm.tsv")]
        [TestCase(@"FileReadingTests\ExternalFileTypes\ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv")]
        public static void TestToppicPrsmReadWrite(string path)
        {
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            var testOutputPath = path.Contains("single")
                ? Path.Combine(directoryPath, "toppic_prsm_single.tsv")
                : Path.Combine(directoryPath, "toppic_prsm.tsv");

            ToppicSearchResultFile file = new ToppicSearchResultFile(filepath);
            file.LoadResults();
            file.WriteResults(testOutputPath);
            var writtenFile = FileReader.ReadFile<ToppicSearchResultFile>(testOutputPath);
            writtenFile.LoadResults();
            Assert.That(File.Exists(testOutputPath));

            // check are equivalent
            for (int i = 0; i < file.Results.Count; i++)
            {
                var original = JsonConvert.SerializeObject(file.Results[i]);
                var written = JsonConvert.SerializeObject(writtenFile.Results[i]);
                Assert.That(original, Is.EqualTo(written));
            }

            var originalLines = File.ReadAllLines(filepath);
            var writtenLines = File.ReadAllLines(testOutputPath);
            Assert.That(writtenLines.Length, Is.EqualTo(originalLines.Length));
            int paramCount = 0;
            for (int i = 0; i < originalLines.Length; i++)
            {
                if (originalLines[i].Contains("********"))
                    paramCount++;
                if (paramCount >= 2)
                    break;
                Assert.That(writtenLines[i], Is.EqualTo(originalLines[i]));
            }

            // test writer still works without specifying extensions
            File.Delete(testOutputPath);
            var testOutputPathWithoutExtension = Path.Combine(directoryPath, "toppic");
            file.WriteResults(testOutputPathWithoutExtension);
            Assert.That(File.Exists(testOutputPath));
        }
    }
}
