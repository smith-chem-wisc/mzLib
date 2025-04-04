using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using Readers;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    public class TestMsPathFinderT
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_TargetResults_IcTarget.tsv", 6)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_DecoyResults_IcDecoy.tsv", 6)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_AllResults_IcTda.tsv", 6)]
        public void TestMsPathFinderTLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsPathFinderTResultFile file = new MsPathFinderTResultFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        public static void TestMsPathFinderTAllResultsFirstAndLastCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                               @"FileReadingTests\ExternalFileTypes\MsPathFinderT_AllResults_IcTda.tsv");
            MsPathFinderTResultFile file = new MsPathFinderTResultFile(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.OneBasedScanNumber, Is.EqualTo(1180));
            Assert.That(first.PreviousResidue, Is.EqualTo('M'));
            Assert.That(first.BaseSequence, Is.EqualTo("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD"));
            Assert.That(first.NextResidue, Is.EqualTo('-'));
            Assert.That(first.Modifications, Is.Empty);
            Assert.That(first.ChemicalFormula.Formula, Is.EqualTo("C442H756N140O156"));
            Assert.That(first.ProteinName, Is.EqualTo("sp|P05114|HMGN1_HUMAN"));
            Assert.That(first.ProteinDescription, Is.EqualTo("Non-histone chromosomal protein HMG-14 OS=Homo sapiens OX=9606 GN=HMGN1 PE=1 SV=3"));
            Assert.That(first.Length, Is.EqualTo(101));
            Assert.That(first.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(first.OneBasedEndResidue, Is.EqualTo(100));
            Assert.That(first.Charge, Is.EqualTo(15));
            Assert.That(first.MostAbundantIsotopeMz, Is.EqualTo(702.8454697));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(10521.55277));
            Assert.That(first.Ms1Features, Is.EqualTo(0));
            Assert.That(first.NumberOfMatchedFragments, Is.EqualTo(61));
            Assert.That(first.Probability, Is.EqualTo(1.0));
            Assert.That(first.SpecEValue, Is.EqualTo(9.99E-308));
            Assert.That(first.EValue, Is.EqualTo(9.99E-308));
            Assert.That(first.QValue, Is.EqualTo(0));
            Assert.That(first.PepQValue, Is.EqualTo(0));
            Assert.That(first.Accession, Is.EqualTo("P05114"));
            Assert.That(first.IsDecoy, Is.EqualTo(false));
            Assert.That(first.FileNameWithoutExtension, Is.EqualTo("MsPathFinderT_AllResults"));

            Assert.That(last.OneBasedScanNumber, Is.EqualTo(1181));
            Assert.That(last.PreviousResidue, Is.EqualTo('M'));
            Assert.That(last.BaseSequence, Is.EqualTo("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD"));
            Assert.That(last.NextResidue, Is.EqualTo('-'));
            Assert.That(last.Modifications, Is.Empty);
            Assert.That(last.ChemicalFormula.Formula, Is.EqualTo("C442H756N140O156"));
            Assert.That(last.ProteinName, Is.EqualTo("sp|P05114|HMGN1_HUMAN"));
            Assert.That(last.ProteinDescription, Is.EqualTo("Non-histone chromosomal protein HMG-14 OS=Homo sapiens OX=9606 GN=HMGN1 PE=1 SV=3"));
            Assert.That(last.Length, Is.EqualTo(101));
            Assert.That(last.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(last.OneBasedEndResidue, Is.EqualTo(100));
            Assert.That(last.Charge, Is.EqualTo(14));
            Assert.That(last.MostAbundantIsotopeMz, Is.EqualTo(752.9767692));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(10521.55277));
            Assert.That(last.Ms1Features, Is.EqualTo(0));
            Assert.That(last.NumberOfMatchedFragments, Is.EqualTo(60));
            Assert.That(last.Probability, Is.EqualTo(1.0));
            Assert.That(last.SpecEValue, Is.EqualTo(9.99E-308));
            Assert.That(last.EValue, Is.EqualTo(9.99E-308));
            Assert.That(last.QValue, Is.EqualTo(0));
            Assert.That(last.PepQValue, Is.EqualTo(0));
            Assert.That(last.Accession, Is.EqualTo("P05114"));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
            Assert.That(last.FileNameWithoutExtension, Is.EqualTo("MsPathFinderT_AllResults"));
        }

        [Test]
        public static void TestMsPathFinderTTargetResultsFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\MsPathFinderT_TargetResults_IcTarget.tsv");
            MsPathFinderTResultFile file = new MsPathFinderTResultFile(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.OneBasedScanNumber, Is.EqualTo(592));
            Assert.That(first.PreviousResidue, Is.EqualTo('M'));
            Assert.That(first.BaseSequence, Is.EqualTo("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK"));
            Assert.That(first.NextResidue, Is.EqualTo('-'));
            Assert.That(first.Modifications, Is.Empty);
            Assert.That(first.ChemicalFormula.Formula, Is.EqualTo("C395H673N129O127"));
            Assert.That(first.ProteinName, Is.EqualTo("sp|P05204|HMGN2_HUMAN"));
            Assert.That(first.ProteinDescription, Is.EqualTo("Non-histone chromosomal protein HMG-17 OS=Homo sapiens OX=9606 GN=HMGN2 PE=1 SV=3"));
            Assert.That(first.Length, Is.EqualTo(91));
            Assert.That(first.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(first.OneBasedEndResidue, Is.EqualTo(90));
            Assert.That(first.Charge, Is.EqualTo(14));
            Assert.That(first.MostAbundantIsotopeMz, Is.EqualTo(662.5096855));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(9256.016953));
            Assert.That(first.Ms1Features, Is.EqualTo(531));
            Assert.That(first.NumberOfMatchedFragments, Is.EqualTo(17));
            Assert.That(first.Probability, Is.EqualTo(0.9999));
            Assert.That(first.SpecEValue, Is.EqualTo(1.642469E-36));
            Assert.That(first.EValue, Is.EqualTo(6.258791E-32));
            Assert.That(first.QValue, Is.EqualTo(0));
            Assert.That(first.PepQValue, Is.EqualTo(0));
            Assert.That(first.Accession, Is.EqualTo("P05204"));
            Assert.That(first.IsDecoy, Is.EqualTo(false));

            Assert.That(last.OneBasedScanNumber, Is.EqualTo(1169));
            Assert.That(last.PreviousResidue, Is.EqualTo('M'));
            Assert.That(last.BaseSequence, Is.EqualTo("PKRKSPENTEGKDGSKVTKQEPTRRSARLSAKPAPPKPEPKPRKTSAKKEPGAKISRGAKGKKEEKQEAGKEGTAPSENGETKAEEAQKTESVDNEGE"));
            Assert.That(last.NextResidue, Is.EqualTo('-'));
            Assert.That(last.Modifications, Is.Empty);
            Assert.That(last.ChemicalFormula.Formula, Is.EqualTo("C442H749N141O156"));
            Assert.That(last.ProteinName, Is.EqualTo("sp|Q15651|HMGN3_HUMAN"));
            Assert.That(last.ProteinDescription, Is.EqualTo("High mobility group nucleosome-binding domain-containing protein 3 OS=Homo sapiens OX=9606 GN=HMGN3 PE=1 SV=2"));
            Assert.That(last.Length, Is.EqualTo(100));
            Assert.That(last.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(last.OneBasedEndResidue, Is.EqualTo(99));
            Assert.That(last.Charge, Is.EqualTo(15));
            Assert.That(last.MostAbundantIsotopeMz, Is.EqualTo(703.3086896));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(10528.50107));
            Assert.That(last.Ms1Features, Is.EqualTo(666));
            Assert.That(last.NumberOfMatchedFragments, Is.EqualTo(7));
            Assert.That(last.Probability, Is.EqualTo(0.161));
            Assert.That(last.SpecEValue, Is.EqualTo(3.706717E-14));
            Assert.That(last.EValue, Is.EqualTo(1.412482E-09));
            Assert.That(last.QValue, Is.EqualTo(0));
            Assert.That(last.PepQValue, Is.EqualTo(0));
            Assert.That(last.Accession, Is.EqualTo("Q15651"));
            Assert.That(last.IsDecoy, Is.EqualTo(false));
        }

        [Test]
        public static void TestMsPathFinderTDecoysFirstAndLastIsCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                               @"FileReadingTests\ExternalFileTypes\MsPathFinderT_DecoyResults_IcDecoy.tsv");
            MsPathFinderTResultFile file = new MsPathFinderTResultFile(filePath);
            var first = file.First();
            var last = file.Last();

            Assert.That(first.OneBasedScanNumber, Is.EqualTo(1180));
            Assert.That(first.PreviousResidue, Is.EqualTo('-'));
            Assert.That(first.BaseSequence, Is.EqualTo("EKTEDKAQDTVAITQIFINGQYGEDAHKDYERGSDRIKKYKVRSGKKLEGTMDEYNIISAKKNVHKATKNRQKYNKINELIGQSNLYSRGF"));
            Assert.That(first.NextResidue, Is.EqualTo('T'));
            Assert.That(first.Modifications, Is.Empty);
            Assert.That(first.ChemicalFormula.Formula, Is.EqualTo("C460H740N136O146S"));
            Assert.That(first.ProteinName, Is.EqualTo("XXX_sp|Q2G1S6|SSL5_STAA8"));
            Assert.That(first.ProteinDescription, Is.EqualTo("Staphylococcal superantigen-like 5 OS=Staphylococcus aureus (strain NCTC 8325"));
            Assert.That(first.Length, Is.EqualTo(235));
            Assert.That(first.OneBasedStartResidue, Is.EqualTo(1));
            Assert.That(first.OneBasedEndResidue, Is.EqualTo(91));
            Assert.That(first.Charge, Is.EqualTo(15));
            Assert.That(first.MostAbundantIsotopeMz, Is.EqualTo(703.9044982));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(10537.4382));
            Assert.That(first.Ms1Features, Is.EqualTo(670));
            Assert.That(first.NumberOfMatchedFragments, Is.EqualTo(5));
            Assert.That(first.Probability, Is.EqualTo(0.0017));
            Assert.That(first.SpecEValue, Is.EqualTo(2.399925E-07));
            Assert.That(first.EValue, Is.EqualTo(0.009145));
            Assert.That(first.QValue, Is.EqualTo(0));
            Assert.That(first.PepQValue, Is.EqualTo(0));
            Assert.That(first.Accession, Is.EqualTo("Q2G1S6"));
            Assert.That(first.IsDecoy, Is.EqualTo(true));

            Assert.That(last.OneBasedScanNumber, Is.EqualTo(1181));
            Assert.That(last.PreviousResidue, Is.EqualTo('N'));
            Assert.That(last.BaseSequence, Is.EqualTo("PAELTKRVQLDMVCRIKEPLRSVKPLPSLAAPSVSSPQLASKSAATKVSMSFSAYMLKMISPVNRAQRCSPQSYRTHVTSPMKRSAVVSEYLSHDP"));
            Assert.That(last.NextResidue, Is.EqualTo('T'));
            Assert.That(last.Modifications, Is.Empty);
            Assert.That(last.ChemicalFormula.Formula, Is.EqualTo("C459H763N133O136S7"));
            Assert.That(last.ProteinName, Is.EqualTo("XXX_sp|Q9PTQ7|DMRT1_CHICK"));
            Assert.That(last.ProteinDescription, Is.EqualTo("Doublesex- and mab-3-related transcription factor 1 OS=Gallus gallus OX=9031 GN=DMRT1 PE=2 SV=2"));
            Assert.That(last.Length, Is.EqualTo(366));
            Assert.That(last.OneBasedStartResidue, Is.EqualTo(2));
            Assert.That(last.OneBasedEndResidue, Is.EqualTo(97));
            Assert.That(last.Charge, Is.EqualTo(14));
            Assert.That(last.MostAbundantIsotopeMz, Is.EqualTo(754.1867306));
            Assert.That(last.MonoisotopicMass, Is.EqualTo(10538.49223));
            Assert.That(last.Ms1Features, Is.EqualTo(671));
            Assert.That(last.NumberOfMatchedFragments, Is.EqualTo(8));
            Assert.That(last.Probability, Is.EqualTo(1.8379E-04));
            Assert.That(last.SpecEValue, Is.EqualTo(2.215391E-06));
            Assert.That(last.EValue, Is.EqualTo(0.08442));
            Assert.That(last.QValue, Is.EqualTo(0));
            Assert.That(last.PepQValue, Is.EqualTo(0));
            Assert.That(last.Accession, Is.EqualTo("Q9PTQ7"));
            Assert.That(last.IsDecoy, Is.EqualTo(true));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_DecoyResults_IcDecoy.tsv", 2)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_TargetResults_IcTarget.tsv", 1)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\MsPathFinderT_AllResults_IcTda.tsv", 3)]
        public static void TestMsPathFinderTAllResultsReadWrite(string path, int fileNum)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            string outpath = Path.Combine(directoryPath, $"MsPathFinderT_Out{fileNum}");

            MsPathFinderTResultFile file = new MsPathFinderTResultFile(filePath);
            file.WriteResults(outpath);

            MsPathFinderTResultFile outFile = new MsPathFinderTResultFile(outpath + path.ParseFileType().GetFileExtension());
            Assert.That(outFile.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file.Count(); i++)
            {
                var original = System.Text.Json.JsonSerializer.Serialize(file.ElementAt(i));
                var written = System.Text.Json.JsonSerializer.Serialize(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }

            MsPathFinderTResultFile file2 = FileReader.ReadFile<MsPathFinderTResultFile>(filePath);
            Assert.That(file2.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file2.Count(); i++)
            {
                var original = System.Text.Json.JsonSerializer.Serialize(file.ElementAt(i));
                var read = System.Text.Json.JsonSerializer.Serialize(file2.ElementAt(i));
                Assert.That(original, Is.EqualTo(read));
            }
        }

        [Test]
        public static void ModificationReading()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\MsPathFinderT_WithMods_IcTda.tsv");
            MsPathFinderTResultFile file = new MsPathFinderTResultFile(path);

            Assert.That(file.Count(), Is.EqualTo(5));

            var result = file.First();
            Assert.That(result.AllModsOneIsNterminus.Count, Is.EqualTo(1));
            Assert.That(result.AllModsOneIsNterminus.First().Value.Target.ToString()[0], Is.EqualTo('M'));
            Assert.That(result.AllModsOneIsNterminus.First().Value.IdWithMotif, Does.Contain("Oxidation"));
            Assert.That(result.AllModsOneIsNterminus.First().Key, Is.EqualTo(37));
            Assert.That(result.AllModsOneIsNterminus.First().Value.MonoisotopicMass, Is.EqualTo(15.99491463).Within(0.0001));
            Assert.That(result.FullSequence, Is.EqualTo("KVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRM[Common Variable:Oxidation on M]QYNRRFVNVVPTFGKKKGPNANS"));

            result = file[1];
            Assert.That(result.AllModsOneIsNterminus.Count, Is.EqualTo(0));
            Assert.That(result.FullSequence, Is.EqualTo("KVHGSLARAGKVRGQTPKVAKQEKKKKKTGRAKRRMQYNRRFVNVVPTFGKKKGPNANS"));

            result = file[2];
            Assert.That(result.AllModsOneIsNterminus.Count, Is.EqualTo(2));
            Assert.That(result.AllModsOneIsNterminus.First().Value.Target.ToString()[0], Is.EqualTo('X'));
            Assert.That(result.AllModsOneIsNterminus.First().Value.IdWithMotif, Does.Contain("Acetyl"));
            Assert.That(result.AllModsOneIsNterminus.First().Key, Is.EqualTo(1));
            Assert.That(result.AllModsOneIsNterminus.First().Value.MonoisotopicMass, Is.EqualTo(42.010565).Within(0.0001));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.Target.ToString()[0], Is.EqualTo('K'));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.IdWithMotif, Does.Contain("Acetyl"));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Key, Is.EqualTo(27));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.MonoisotopicMass, Is.EqualTo(42.010565).Within(0.0001));
            Assert.That(result.FullSequence, Is.EqualTo("[Common Biological:Acetylation on X]TEIKEKSVAELNALLKEKKVLLFTLK[Common Biological:Acetylation on K]QKLKTMQLSNPNEIRALKKEIARINTAISASK"));

            result = file[3];
            Assert.That(result.AllModsOneIsNterminus.Count, Is.EqualTo(2));
            Assert.That(result.AllModsOneIsNterminus.First().Value.Target.ToString()[0], Is.EqualTo('K'));
            Assert.That(result.AllModsOneIsNterminus.First().Value.IdWithMotif, Does.Contain("Methyl"));
            Assert.That(result.AllModsOneIsNterminus.First().Key, Is.EqualTo(25));
            Assert.That(result.AllModsOneIsNterminus.First().Value.MonoisotopicMass, Is.EqualTo(14.01565).Within(0.0001));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.Target.ToString()[0], Is.EqualTo('K'));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.IdWithMotif, Does.Contain("Acetyl"));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Key, Is.EqualTo(51));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.MonoisotopicMass, Is.EqualTo(42.010565).Within(0.0001));
            Assert.That(result.FullSequence, Is.EqualTo("CNRQAVAGQLLPSTWSLHAHGLAK[Common Biological:Methylation on K]EAPILPVKKISRSCSVNNKVSKKTTK[Common Biological:Acetylation on K]PPTLRSFLSPI"));

            result = file[4];
            Assert.That(result.AllModsOneIsNterminus.Count, Is.EqualTo(2));
            Assert.That(result.AllModsOneIsNterminus.First().Value.Target.ToString()[0], Is.EqualTo('K'));
            Assert.That(result.AllModsOneIsNterminus.First().Value.IdWithMotif, Does.Contain("Acetyl"));
            Assert.That(result.AllModsOneIsNterminus.First().Key, Is.EqualTo(9));
            Assert.That(result.AllModsOneIsNterminus.First().Value.MonoisotopicMass, Is.EqualTo(42.010565).Within(0.0001));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.Target.ToString()[0], Is.EqualTo('K'));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.IdWithMotif, Does.Contain("Acetyl"));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Key, Is.EqualTo(53));
            Assert.That(result.AllModsOneIsNterminus.Skip(1).First().Value.MonoisotopicMass, Is.EqualTo(42.010565).Within(0.0001));
            Assert.That(result.FullSequence, Is.EqualTo("SFFDHLQK[Common Biological:Acetylation on K]KGVGAIQAQKVQIRKVERKPSKVVSSSSSSSIAKPQRRLDTVSK[Common Biological:Acetylation on K]PVAARRSA"));
        }
    }
}
