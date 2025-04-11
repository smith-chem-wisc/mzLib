using NUnit.Framework;
using Proteomics;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using System.Diagnostics.Metrics;
using System.Threading;
using MzLibUtil.NoiseEstimation;
using System.Windows.Data;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    internal class TestMsFraggerResultFiles
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv", 5)]
        public void TestMsFraggerPsmLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv")]
        public void TestAddProteinGroupInfoCountCorrect (string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            var allResults = file.ToList();

            // one protein associated with given results, list should only contain this one element 
            Assert.That(allResults[0].ProteinGroupInfos.Count, Is.EqualTo(1));
            // two proteins associated with given results, list should contain two elements
            Assert.That(allResults[2].ProteinGroupInfos.Count, Is.EqualTo(2));

        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1individual_peptide.tsv", 7)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1combined_peptide.tsv", 6)]
        public static void TestMsFraggerPeptideLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPeptideFile file = new MsFraggerPeptideFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(filePath));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1individual_protein.tsv", 8)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1combined_protein.tsv", 4)]
        public static void TestMsFraggerProteinLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerProteinFile file = new MsFraggerProteinFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(filePath));
        }

        [Test]
        public void TestMsFraggerPsmFirstAndLastCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv");
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            MsFraggerPsm first = file.First();
            MsFraggerPsm last = file.Last();
            
            Assert.That(first.Spectrum, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1.00003.00003.2"));
            Assert.That(first.SpectrumFilePath, Is.EqualTo(@"D:\Projects\Chimeras\Mann_11cell_analysis\Hela\MsFragger\Hela_1_1\interact-20100611_Velos1_TaGe_SA_Hela_1.pep.xml"));
            Assert.That(first.BaseSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(first.ModifiedPeptide, Is.EqualTo(""));
            Assert.That(first.ExtendedSequence, Is.EqualTo("KAGGTKPK.KPVGAAK.KPKKAAGG"));
            Assert.That(first.PreviousAminoAcid, Is.EqualTo('K'));
            Assert.That(first.NextAminoAcid, Is.EqualTo('K'));
            Assert.That(first.PeptideLength, Is.EqualTo(7));
            Assert.That(first.Charge, Is.EqualTo(2));
            Assert.That(first.RetentionTime, Is.EqualTo(1.9398));
            Assert.That(first.ObservedMass, Is.EqualTo(669.4164));
            Assert.That(first.CalibratedObservedMass, Is.EqualTo(669.4123));
            Assert.That(first.ObservedMz, Is.EqualTo(335.7155));
            Assert.That(first.CalibratedObservedMz, Is.EqualTo(335.7135));
            Assert.That(first.CalculatedPeptideMass, Is.EqualTo(669.4173));
            Assert.That(first.CalculatedMz, Is.EqualTo(335.7159));
            Assert.That(first.DeltaMass, Is.EqualTo(-0.0049));
            Assert.That(first.Expectation, Is.EqualTo(0.01675));
            Assert.That(first.HyperScore, Is.EqualTo(20.604));
            Assert.That(first.NextScore, Is.EqualTo(12.521));
            Assert.That(first.PeptideProphetProbability, Is.EqualTo(0.9176));
            Assert.That(first.NumberOfEnzymaticTermini, Is.EqualTo(2));
            Assert.That(first.NumberOfMissedCleavages, Is.EqualTo(1));
            Assert.That(first.ProteinStart, Is.EqualTo(130));
            Assert.That(first.ProteinEnd, Is.EqualTo(136));
            Assert.That(first.Intensity, Is.EqualTo(0));
            Assert.That(first.AssignedModifications, Is.EqualTo(""));
            Assert.That(first.ObservedModifications, Is.EqualTo(""));
            Assert.That(first.Purity, Is.EqualTo(0));
            Assert.That(first.IsUnique, Is.EqualTo(true));
            Assert.That(first.Protein, Is.EqualTo("sp|P16403|H12_HUMAN"));
            Assert.That(first.ProteinAccession, Is.EqualTo("P16403"));
            Assert.That(first.EntryName, Is.EqualTo("H12_HUMAN"));
            Assert.That(first.Gene, Is.EqualTo("H1-2"));
            Assert.That(first.ProteinDescription, Is.EqualTo("Histone H1.2"));
            Assert.That(first.MappedGenes, Is.EqualTo(""));
            Assert.That(first.MappedProteins, Is.EqualTo(""));
            Assert.That(first.FileNameWithoutExtension, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1"));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(3));
            Assert.That(first.OneBasedScanNumber, Is.EqualTo(3));

            Assert.That(last.Spectrum, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1.00018.00018.2"));
            Assert.That(last.SpectrumFilePath, Is.EqualTo(@"D:\Projects\Chimeras\Mann_11cell_analysis\Hela\MsFragger\Hela_1_1\interact-20100611_Velos1_TaGe_SA_Hela_1.pep.xml"));
            Assert.That(last.BaseSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(last.ModifiedPeptide, Is.EqualTo(""));
            Assert.That(last.ExtendedSequence, Is.EqualTo("GTAIKNGK.VVTHGGR.VIAVTAIR"));
            Assert.That(last.PreviousAminoAcid, Is.EqualTo('K'));
            Assert.That(last.NextAminoAcid, Is.EqualTo('V'));
            Assert.That(last.PeptideLength, Is.EqualTo(7));
            Assert.That(last.Charge, Is.EqualTo(2));
            Assert.That(last.RetentionTime, Is.EqualTo(19.114));
            Assert.That(last.ObservedMass, Is.EqualTo(724.3984));
            Assert.That(last.CalibratedObservedMass, Is.EqualTo
            (724.3949));
            Assert.That(last.ObservedMz, Is.EqualTo(363.2065));
            Assert.That(last.CalibratedObservedMz, Is.EqualTo(363.2047));
            Assert.That(last.CalculatedPeptideMass, Is.EqualTo
            (724.398));
            Assert.That(last.CalculatedMz, Is.EqualTo(363.2063));
            Assert.That(last.DeltaMass, Is.EqualTo(-0.0031));
            Assert.That(last.Expectation, Is.EqualTo(0.7595));
            Assert.That(last.HyperScore, Is.EqualTo(11.978));
            Assert.That(last.NextScore, Is.EqualTo(0));
            Assert.That(last.PeptideProphetProbability, Is.EqualTo(0.8144));
            Assert.That(last.NumberOfEnzymaticTermini, Is.EqualTo(2));
            Assert.That(last.NumberOfMissedCleavages, Is.EqualTo(0));
            Assert.That(last.ProteinStart, Is.EqualTo(379));
            Assert.That(last.ProteinEnd, Is.EqualTo(385));
            Assert.That(last.Intensity, Is.EqualTo(0));
            Assert.That(last.AssignedModifications, Is.EqualTo(""));
            Assert.That(last.ObservedModifications, Is.EqualTo(""));
            Assert.That(last.Purity, Is.EqualTo(0));
            Assert.That(last.IsUnique, Is.EqualTo(true));
            Assert.That(last.Protein, Is.EqualTo("sp|P22102|PUR2_HUMAN"));
            Assert.That(last.ProteinAccession, Is.EqualTo("P22102"));
            Assert.That(last.EntryName, Is.EqualTo("PUR2_HUMAN"));
            Assert.That(last.Gene, Is.EqualTo("GART"));
            Assert.That(last.ProteinDescription, Is.EqualTo("Trifunctional purine biosynthetic protein adenosine-3"));
            Assert.That(last.MappedGenes, Is.EqualTo(""));
            Assert.That(last.MappedProteins, Is.EqualTo(""));
            Assert.That(last.FileNameWithoutExtension, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1"));
            Assert.That(last.OneBasedScanNumber, Is.EqualTo(18));
        }

        [Test]
        public void TestMsFraggerIndividualPeptideFirstAndLastAreCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                               @"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1individual_peptide.tsv");
            MsFraggerPeptideFile file = new MsFraggerPeptideFile(filePath);
            MsFraggerPeptide first = file.First();
            MsFraggerPeptide last = file.Last();

            Assert.That(first.BaseSequence, Is.EqualTo("AAAAVVEFQR"));
            Assert.That(first.PreviousAminoAcid, Is.EqualTo('M'));
            Assert.That(first.NextAminoAcid, Is.EqualTo('A'));
            Assert.That(first.PeptideLength, Is.EqualTo(10));
            Assert.That(first.OneBasedStartResidueInProtein, Is.EqualTo(2));
            Assert.That(first.OneBasedEndResidueInProtein, Is.EqualTo(11));
            Assert.That(first.Charge, Is.EqualTo(new int[] { 2 }));
            Assert.That(first.Probability, Is.EqualTo(0.9473));
            Assert.That(first.SpectralCount, Is.EqualTo(1));
            Assert.That(first.Intensity, Is.EqualTo(0));
            Assert.That(first.AssignedModifications, Is.EqualTo(new string[] { "N-term(42.0106)" }));
            Assert.That(first.ObservedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(first.Protein, Is.EqualTo("sp|O00231|PSD11_HUMAN"));
            Assert.That(first.ProteinAccession, Is.EqualTo("O00231"));
            Assert.That(first.ProteinName, Is.EqualTo("PSD11_HUMAN"));
            Assert.That(first.Gene, Is.EqualTo("PSMD11"));
            Assert.That(first.ProteinDescription, Is.EqualTo("26S proteasome non-ATPase regulatory subunit 11"));
            Assert.That(first.MappedGenes, Is.EqualTo(""));
            Assert.That(first.MappedProteins, Is.EqualTo(""));

            Assert.That(last.BaseSequence, Is.EqualTo("AACLCFR"));
            Assert.That(last.PreviousAminoAcid, Is.EqualTo('R'));
            Assert.That(last.NextAminoAcid, Is.EqualTo('S'));
            Assert.That(last.PeptideLength, Is.EqualTo(7));
            Assert.That(last.OneBasedStartResidueInProtein, Is.EqualTo(21));
            Assert.That(last.OneBasedEndResidueInProtein, Is.EqualTo(27));
            Assert.That(last.Charge, Is.EqualTo(new int[] { 2 }));
            Assert.That(last.Probability, Is.EqualTo(0.9553));
            Assert.That(last.SpectralCount, Is.EqualTo(1));
            Assert.That(last.Intensity, Is.EqualTo(0));
            Assert.That(last.AssignedModifications, Is.EqualTo(new [] { "3C(57.0214)", "5C(57.0214)" }));
            Assert.That(last.ObservedModifications, Is.EqualTo(new string[] { }));
            Assert.That(last.Protein, Is.EqualTo("sp|Q9NZJ9|NUDT4_HUMAN"));
            Assert.That(last.ProteinAccession, Is.EqualTo("Q9NZJ9"));
            Assert.That(last.ProteinName, Is.EqualTo("NUDT4_HUMAN"));
            Assert.That(last.Gene, Is.EqualTo("NUDT4"));
            Assert.That(last.ProteinDescription, Is.EqualTo("Diphosphoinositol polyphosphate phosphohydrolase 2"));
            Assert.That(last.MappedGenes, Is.EqualTo(new [] {"NUDT10", "NUDT11", "NUDT3", "NUDT4B"}));
            Assert.That(last.MappedProteins, Is.EqualTo(new []{"sp|A0A024RBG1|NUD4B_HUMAN", "sp|O95989|NUDT3_HUMAN", "sp|Q8NFP7|NUD10_HUMAN", "sp|Q96G61|NUD11_HUMAN"}));
        }

        [Test]
        public static void TestMsFraggerIndividualProteinFirstAndLastAreCorrect()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                               @"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1individual_protein.tsv");
            MsFraggerProteinFile file = new MsFraggerProteinFile(path);
            MsFraggerProtein first = file.First();
            MsFraggerProtein last = file.Last();

            Assert.That(first.Protein, Is.EqualTo("contam_sp|P00761|TRYP_PIG"));
            Assert.That(first.Accession, Is.EqualTo("P00761"));
            Assert.That(first.AccessionOrganism, Is.EqualTo("TRYP_PIG"));
            Assert.That(first.Gene, Is.EqualTo(""));
            Assert.That(first.Length, Is.EqualTo(231));
            Assert.That(first.Organism, Is.EqualTo("Sus scrofa"));
            Assert.That(first.Description, Is.EqualTo("Trypsin"));
            Assert.That(first.ProteinExistence, Is.EqualTo("1:Experimental evidence at protein level"));
            Assert.That(first.Coverage, Is.EqualTo(31.60));
            Assert.That(first.ProteinProbability, Is.EqualTo(1.0000));
            Assert.That(first.TopPeptideProbability, Is.EqualTo(0.9990));
            Assert.That(first.TotalPeptides, Is.EqualTo(6));
            Assert.That(first.UniquePeptides, Is.EqualTo(6));
            Assert.That(first.RazorPeptides, Is.EqualTo(6));
            Assert.That(first.TotalSpectralCount, Is.EqualTo(14));
            Assert.That(first.UniqueSpectralCount, Is.EqualTo(14));
            Assert.That(first.RazorSpectralCount, Is.EqualTo(14));
            Assert.That(first.TotalIntensity, Is.EqualTo(0));
            Assert.That(first.UniqueIntensity, Is.EqualTo(0));
            Assert.That(first.RazorIntensity, Is.EqualTo(0));
            Assert.That(first.RazorAssignedModifications, Is.EqualTo(new string[] { "17M(15.9949)", "2C(57.0214)" }));
            Assert.That(first.RazorObservedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(first.IndistinguishableProteins, Is.EqualTo(new string[] {}));

            Assert.That(last.Protein, Is.EqualTo("sp|A0MZ66|SHOT1_HUMAN"));
            Assert.That(last.Accession, Is.EqualTo("A0MZ66"));
            Assert.That(last.AccessionOrganism, Is.EqualTo("SHOT1_HUMAN"));
            Assert.That(last.Gene, Is.EqualTo("SHTN1"));
            Assert.That(last.Length, Is.EqualTo(631));
            Assert.That(last.Organism, Is.EqualTo("Homo sapiens"));
            Assert.That(last.Description, Is.EqualTo("Shootin-1"));
            Assert.That(last.ProteinExistence, Is.EqualTo("1:Experimental evidence at protein level"));
            Assert.That(last.Coverage, Is.EqualTo(2.38));
            Assert.That(last.ProteinProbability, Is.EqualTo(1.0000));
            Assert.That(last.TopPeptideProbability, Is.EqualTo(0.9990));
            Assert.That(last.TotalPeptides, Is.EqualTo(2));
            Assert.That(last.UniquePeptides, Is.EqualTo(2));
            Assert.That(last.RazorPeptides, Is.EqualTo(2));
            Assert.That(last.TotalSpectralCount, Is.EqualTo(2));
            Assert.That(last.UniqueSpectralCount, Is.EqualTo(2));
            Assert.That(last.RazorSpectralCount, Is.EqualTo(2));
            Assert.That(last.TotalIntensity, Is.EqualTo(0));
            Assert.That(last.UniqueIntensity, Is.EqualTo(0));
            Assert.That(last.RazorIntensity, Is.EqualTo(0));
            Assert.That(last.RazorAssignedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(last.RazorObservedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(last.IndistinguishableProteins, Is.EqualTo(new string[] { }));
        }

        [Test]
        public void TestMsFraggerCombinedPeptideFirstAndLastAreCorrect()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                                              @"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1combined_peptide.tsv");
            MsFraggerPeptideFile file = new MsFraggerPeptideFile(path);
            MsFraggerPeptide first = file.First();
            MsFraggerPeptide last = file.Last();

            Assert.That(first.BaseSequence, Is.EqualTo("AAAAAAAAAAAAAAAGAGAGAK"));
            Assert.That(first.PreviousAminoAcid, Is.EqualTo(default(char)));
            Assert.That(first.NextAminoAcid, Is.EqualTo(default(char)));
            Assert.That(first.PeptideLength, Is.EqualTo(22));
            Assert.That(first.OneBasedStartResidueInProtein, Is.EqualTo(0));
            Assert.That(first.OneBasedEndResidueInProtein, Is.EqualTo(0));
            Assert.That(first.Charge, Is.EqualTo(new int[] { 2, 3 }));
            Assert.That(first.Probability, Is.EqualTo(1.000000));
            Assert.That(first.SpectralCount, Is.EqualTo(0));
            Assert.That(first.Intensity, Is.EqualTo(0.0000));
            Assert.That(first.AssignedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(first.ObservedModifications, Is.Null);
            Assert.That(first.Protein, Is.EqualTo("sp|P55011|S12A2_HUMAN"));
            Assert.That(first.ProteinAccession, Is.EqualTo("P55011"));
            Assert.That(first.ProteinName, Is.EqualTo("S12A2_HUMAN"));
            Assert.That(first.Gene, Is.EqualTo("SLC12A2"));
            Assert.That(first.ProteinDescription, Is.EqualTo("Solute carrier family 12 member 2"));
            Assert.That(first.MappedGenes, Is.Null);
            Assert.That(first.MappedProteins, Is.Null);

            Assert.That(last.BaseSequence, Is.EqualTo("AAAAAAGAASGLPGPVAQGLK"));
            Assert.That(last.PreviousAminoAcid, Is.EqualTo(default(char)));
            Assert.That(last.NextAminoAcid, Is.EqualTo(default(char)));
            Assert.That(last.PeptideLength, Is.EqualTo(21));
            Assert.That(last.OneBasedStartResidueInProtein, Is.EqualTo(0));
            Assert.That(last.OneBasedEndResidueInProtein, Is.EqualTo(0));
            Assert.That(last.Charge, Is.EqualTo(new int[] { 3,2 }));
            Assert.That(last.Probability, Is.EqualTo(1.000000));
            Assert.That(last.SpectralCount, Is.EqualTo(0));
            Assert.That(last.Intensity, Is.EqualTo(0.0000));
            Assert.That(last.AssignedModifications, Is.EqualTo(new string[] { "42.010600" }));
            Assert.That(last.ObservedModifications, Is.Null);
            Assert.That(last.Protein, Is.EqualTo("sp|Q96P70|IPO9_HUMAN"));
            Assert.That(last.ProteinAccession, Is.EqualTo("Q96P70"));
            Assert.That(last.ProteinName, Is.EqualTo("IPO9_HUMAN"));
            Assert.That(last.Gene, Is.EqualTo("IPO9"));
            Assert.That(last.ProteinDescription, Is.EqualTo("Importin-9"));
            Assert.That(last.MappedGenes, Is.Null);
            Assert.That(last.MappedProteins, Is.Null);
        }

        [Test]
        public static void TestMsFraggerCombinedProteinFirstAndLastAreCorrect()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                                                             @"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1combined_protein.tsv");
            MsFraggerProteinFile file = new MsFraggerProteinFile(path);
            MsFraggerProtein first = file.First();
            MsFraggerProtein last = file.Last();

            Assert.That(first.Protein, Is.EqualTo("contam_sp|P00761|TRYP_PIG"));
            Assert.That(first.Accession, Is.EqualTo("P00761")); 
            Assert.That(first.AccessionOrganism, Is.EqualTo("TRYP_PIG"));
            Assert.That(first.Gene, Is.EqualTo(""));
            Assert.That(first.Length, Is.EqualTo(231));
            Assert.That(first.Organism, Is.EqualTo("Sus scrofa"));
            Assert.That(first.Description, Is.EqualTo("Trypsin"));
            Assert.That(first.ProteinExistence, Is.EqualTo("1:Experimental evidence at protein level"));
            Assert.That(first.Coverage, Is.EqualTo(0));
            Assert.That(first.ProteinProbability, Is.EqualTo(1.0000));
            Assert.That(first.TopPeptideProbability, Is.EqualTo(0.9990));
            Assert.That(first.TotalPeptides, Is.EqualTo(12));
            Assert.That(first.UniquePeptides, Is.EqualTo(0));
            Assert.That(first.RazorPeptides, Is.EqualTo(0));
            Assert.That(first.TotalSpectralCount, Is.EqualTo(550));
            Assert.That(first.UniqueSpectralCount, Is.EqualTo(549));
            Assert.That(first.RazorSpectralCount, Is.EqualTo(0));
            Assert.That(first.TotalIntensity, Is.EqualTo(0));
            Assert.That(first.UniqueIntensity, Is.EqualTo(0));
            Assert.That(first.RazorIntensity, Is.EqualTo(0));
            Assert.That(first.RazorAssignedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(first.RazorObservedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(first.IndistinguishableProteins, Is.EqualTo(new string[] { }));

            Assert.That(last.Protein, Is.EqualTo("sp|A0AV96|RBM47_HUMAN"));
            Assert.That(last.Accession, Is.EqualTo("A0AV96"));
            Assert.That(last.AccessionOrganism, Is.EqualTo("RBM47_HUMAN"));
            Assert.That(last.Gene, Is.EqualTo("RBM47"));
            Assert.That(last.Length, Is.EqualTo(593));
            Assert.That(last.Organism, Is.EqualTo("Homo sapiens"));
            Assert.That(last.Description, Is.EqualTo("RNA-binding protein 47"));
            Assert.That(last.ProteinExistence, Is.EqualTo("1:Experimental evidence at protein level"));
            Assert.That(last.Coverage, Is.EqualTo(0));
            Assert.That(last.ProteinProbability, Is.EqualTo(1.0000));
            Assert.That(last.TopPeptideProbability, Is.EqualTo(0.9990));
            Assert.That(last.TotalPeptides, Is.EqualTo(15));
            Assert.That(last.UniquePeptides, Is.EqualTo(0));
            Assert.That(last.RazorPeptides, Is.EqualTo(0));
            Assert.That(last.TotalSpectralCount, Is.EqualTo(45));
            Assert.That(last.UniqueSpectralCount, Is.EqualTo(38));
            Assert.That(last.RazorSpectralCount, Is.EqualTo(0));
            Assert.That(last.TotalIntensity, Is.EqualTo(0));
            Assert.That(last.UniqueIntensity, Is.EqualTo(0));
            Assert.That(last.RazorIntensity, Is.EqualTo(0));
            Assert.That(last.RazorAssignedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(last.RazorObservedModifications, Is.EqualTo(new string[] {  }));
            Assert.That(last.IndistinguishableProteins, Is.EqualTo(new string[] { }));
        }




        [Test]
        public static void TestMsFraggerPsmReadWrite()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv");
            string outPath = Path.Combine(directoryPath, "testFraggerPsm");

            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            file.WriteResults(outPath);
            
            MsFraggerPsmFile outFile = new MsFraggerPsmFile(outPath + "psm.tsv");
            Assert.That(outFile.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }

            MsFraggerPsmFile file2 = FileReader.ReadFile<MsFraggerPsmFile>(filePath);
            Assert.That(file2.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file2.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file2.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }
        }


        [Test]
        public static void TestMsFraggerPeptideReadWrite()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\FraggerPeptide_FragPipev21.1individual_peptide.tsv");
            string outPath = Path.Combine(directoryPath, "testFraggerPeptide");

            var file = new MsFraggerPeptideFile(filePath);
            file.WriteResults(outPath);

            var outFile = new MsFraggerPeptideFile(outPath + "peptide.tsv");
            Assert.That(outFile.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }

            MsFraggerPeptideFile file2 = FileReader.ReadFile<MsFraggerPeptideFile>(filePath);
            Assert.That(file2.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file2.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file2.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }
        }

        [Test]
        public static void TestMsFraggerProteinReadWrite()
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\FraggerProtein_FragPipev21.1individual_protein.tsv");
            var outpath = Path.Combine(directoryPath, "testFraggerProtein");

            var file = new MsFraggerProteinFile(filePath);
            file.WriteResults(outpath);

            var outFile = new MsFraggerProteinFile(outpath + "protein.tsv");
            Assert.That(outFile.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }

            MsFraggerProteinFile file2 = FileReader.ReadFile<MsFraggerProteinFile>(filePath);
            Assert.That(file2.Count(), Is.EqualTo(file.Count()));
            for (int i = 0; i < file2.Count(); i++)
            {
                var original = JsonConvert.SerializeObject(file2.ElementAt(i));
                var written = JsonConvert.SerializeObject(outFile.ElementAt(i));
                Assert.That(original, Is.EqualTo(written));
            }
        }
    }

   
}
