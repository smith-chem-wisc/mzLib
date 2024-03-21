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

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    internal class TestMsFraggerPsm
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
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1.tsv", 5)]
        public void TestMsFraggerPsmeLoadsAndCountCorrect(string path, int count)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        public void TestMsFraggerPsmFirstAndLastCorrect()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1.tsv");
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            MsFraggerPsm first = file.First();
            MsFraggerPsm last = file.Last();
            
            Assert.That(first.Spectrum, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1.00003.00003.2"));
            Assert.That(first.SpectrumFilePath, Is.EqualTo(@"D:\Projects\Chimeras\Mann_11cell_analysis\Hela\MsFragger\Hela_1_1\interact-20100611_Velos1_TaGe_SA_Hela_1.pep.xml"));
            Assert.That(first.BaseSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(first.FullSequence, Is.EqualTo(""));
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

            Assert.That(last.Spectrum, Is.EqualTo("20100611_Velos1_TaGe_SA_Hela_1.00018.00018.2"));
            Assert.That(last.SpectrumFilePath, Is.EqualTo(@"D:\Projects\Chimeras\Mann_11cell_analysis\Hela\MsFragger\Hela_1_1\interact-20100611_Velos1_TaGe_SA_Hela_1.pep.xml"));
            Assert.That(last.BaseSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(last.FullSequence, Is.EqualTo(""));
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
        public static void TestMsFraggerReadWrite()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1.tsv");
            string outPath = Path.Combine(directoryPath, "testFragger");

            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            file.WriteResults(outPath);
            
            MsFraggerPsmFile outFile = new MsFraggerPsmFile(outPath + ".tsv");
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
    }

   
}
