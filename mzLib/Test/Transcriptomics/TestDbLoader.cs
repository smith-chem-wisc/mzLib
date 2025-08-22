using NUnit.Framework;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases.Transcriptomics;
using UsefulProteomicsDatabases;
using Transcriptomics;
using Omics;

namespace Test.Transcriptomics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestDbLoader
    {
        public static string ModomicsUnmodifedFastaPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta");

        /// <summary>
        /// Detect the headertype of the test cases
        /// </summary>
        private static IEnumerable<(string, RnaFastaHeaderType)> DetectHeaderTestCases =>
            new List<(string, RnaFastaHeaderType)>
            {
                (Path.Combine(TestContext.CurrentContext.TestDirectory, "DoubleProtease.tsv"), RnaFastaHeaderType.Unknown),
                (ModomicsUnmodifedFastaPath, RnaFastaHeaderType.Modomics),
                (Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta"), RnaFastaHeaderType.Modomics),
                
            };

        /// <summary>
        /// Test the correctness of checking headertype
        /// </summary>
        /// <param name="testData"></param>
        [Test]
        [TestCaseSource(nameof(DetectHeaderTestCases))]
        public static void TestDetectHeaderType((string dbPath, RnaFastaHeaderType headerType) testData)
        {
            string line = File.ReadLines(testData.dbPath).First();
            if (char.IsDigit(line.First()))
            {
                line = File.ReadLines(testData.dbPath).Skip(1).First();
            }
            var type = RnaDbLoader.DetectRnaFastaHeaderType(line);
            Assert.That(testData.headerType, Is.EqualTo(type));
        }


        [Test]
        [TestCase("ModomicsUnmodifiedTrimmed.fasta")]
        [TestCase("ModomicsUnmodifiedTrimmed.fasta.gz")]
        public static void TestModomicsUnmodifiedFasta(string databaseFileName)
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData",
                databaseFileName);
            var oligos = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false,
                out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(5));
            Assert.That(oligos.First().BaseSequence,
                Is.EqualTo("GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA"));
            Assert.That(oligos.First().Name, Is.EqualTo("tdbR00000010"));
            Assert.That(oligos.First().Accession, Is.EqualTo("SO:0000254"));
            Assert.That(oligos.First().Organism, Is.EqualTo("Escherichia coli"));
            Assert.That(oligos.First().DatabaseFilePath, Is.EqualTo(dbPath));
            Assert.That(oligos.First().IsContaminant, Is.False);
            Assert.That(oligos.First().IsDecoy, Is.False);
            Assert.That(oligos.First().AdditionalDatabaseFields!.Count, Is.EqualTo(5));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Id"], Is.EqualTo("1"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Type"], Is.EqualTo("tRNA"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Subtype"], Is.EqualTo("Ala"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Feature"], Is.EqualTo("VGC"));
            Assert.That(oligos.First().AdditionalDatabaseFields!["Cellular Localization"], Is.EqualTo("prokaryotic cytosol"));
        }

        [Test]
        public static void TestContaminantFollowsThrough()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, true,
                               out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(5));
            Assert.That(oligos.First().BaseSequence,
                               Is.EqualTo("GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA"));
            Assert.That(oligos.All(p => p.IsContaminant));
            Assert.That(oligos.All(p => !p.IsDecoy));
        }

        [Test]
        public static void TestNotGeneratingTargetsOrDecoys()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, false, DecoyType.None, true,
                out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestFastaWithCustomIdentifier()
        {
            var oligos = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.Reverse, true,
                out var errors, decoyIdentifier: "rev");

            foreach (var rna in oligos)
            {
                if (!rna.IsDecoy) continue;

                Assert.That(rna.Accession, Does.StartWith("rev"));
                Assert.That(rna.Accession, Does.Not.StartWith("DECOY"));
            }
        }

        [Test]
        public static void TestXmlWriterReader()
        {
            var rna = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, false, out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));

            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   G\r\nCF   C1H2\r\n" + @"//";
            var methylG = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> modsOut).First();

            Dictionary<string, HashSet<Tuple<int, Modification>>> mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            mods.Add("SO:0000254", new HashSet<Tuple<int, Modification>>()
            {
                new Tuple<int, Modification>(1, methylG),
                new Tuple<int, Modification>(3, methylG)
            });

            IDictionary<int, List<Modification>> simpleModDictionary = new Dictionary<int, List<Modification>>();
            simpleModDictionary.Add(3, new List<Modification>() { methylG });
            simpleModDictionary.Add(4, new List<Modification>() { methylG });

            RNA newRna = (RNA)rna.First().CloneWithNewSequenceAndMods(
                "GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA",
                simpleModDictionary);
            rna.RemoveAt(0);
            rna.Add(newRna);
            string outpath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.xml");

            var xml = ProteinDbWriter.WriteXmlDatabase(rna, outpath);

            var temp = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.None, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods);

            Assert.That(unknownMods.Count, Is.EqualTo(0));
            Assert.That(temp.Count, Is.EqualTo(5));
            var first = temp.Last();
            var loadedMods = first.OneBasedPossibleLocalizedModifications;
            Assert.That(loadedMods.Count, Is.EqualTo(2));
            Assert.That(loadedMods[3].Count, Is.EqualTo(1));
            Assert.That(loadedMods[4].Count, Is.EqualTo(1));
            Assert.That(loadedMods[3].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
            Assert.That(loadedMods[4].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
        }

        [Test]
        public static void TestXmlWriterReaderAsBioPolymer()
        {
            var rna = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, false, out var errors)
                .Cast<IBioPolymer>().ToList();
            Assert.That(errors.Count, Is.EqualTo(0));

            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   G\r\nCF   C1H2\r\n" + @"//";
            var methylG = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> modsOut).First();

            Dictionary<string, HashSet<Tuple<int, Modification>>> mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            mods.Add("SO:0000254", new HashSet<Tuple<int, Modification>>()
            {
                new Tuple<int, Modification>(1, methylG),
                new Tuple<int, Modification>(3, methylG)
            });
            IDictionary<int, List<Modification>> simpleModDictionary = new Dictionary<int, List<Modification>>();
            simpleModDictionary.Add(3, new List<Modification>() { methylG });
            simpleModDictionary.Add(4, new List<Modification>() { methylG });

            RNA newRna = (RNA)rna.First().CloneWithNewSequenceAndMods(
                "GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA",
                simpleModDictionary);
            rna.RemoveAt(0);
            rna.Add(newRna);
            string outpath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed2.xml");

            var xml = ProteinDbWriter.WriteXmlDatabase(rna, outpath);
            var temp = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.None, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods);

            Assert.That(unknownMods.Count, Is.EqualTo(0));
            Assert.That(temp.Count, Is.EqualTo(5));
            var first = temp.Last();
            var loadedMods = first.OneBasedPossibleLocalizedModifications;
            Assert.That(loadedMods.Count, Is.EqualTo(2));
            Assert.That(loadedMods[3].Count, Is.EqualTo(1));
            Assert.That(loadedMods[4].Count, Is.EqualTo(1));
            Assert.That(loadedMods[3].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
            Assert.That(loadedMods[4].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
        }

        [Test]
        public static void TestXmlWithCustomIdentifier()
        {
            var rna = RnaDbLoader.LoadRnaFasta(ModomicsUnmodifedFastaPath, true, DecoyType.None, false, out var errors);
            Assert.That(errors.Count, Is.EqualTo(0));

            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   G\r\nCF   C1H2\r\n" + @"//";
            var methylG = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> modsOut).First();

            Dictionary<string, HashSet<Tuple<int, Modification>>> mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            mods.Add("SO:0000254", new HashSet<Tuple<int, Modification>>()
            {
                new Tuple<int, Modification>(1, methylG),
                new Tuple<int, Modification>(3, methylG)
            });

            string outpath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed2.xml");
            var xml = ProteinDbWriter.WriteXmlDatabase(rna, outpath);
            rna = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.Reverse, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods, decoyIdentifier: "rev");

            foreach (var oligo in rna)
            {
                if (!oligo.IsDecoy) continue;

                Assert.That(oligo.Accession, Does.StartWith("rev"));
                Assert.That(oligo.Accession, Does.Not.StartWith("DECOY"));
            }
        }

        [Test]
        [TestCase("ATCG", "AUCG", true)]
        [TestCase("ATCG", "UAGC", false)]
        [TestCase("ATCGZ", "AUCGZ", true)]
        [TestCase("ATCGZ", "UAGCZ", false)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACT", "AUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACU", true)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACT", "UAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGA", false)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACT", "AUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACU", true)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACT", "UAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGA", false)]
        public static void TestTranscribe(string input, string expected, bool isCodingStrand)
        {
            Assert.That(input.Transcribe(isCodingStrand), Is.EqualTo(expected));
        }

        [Test]
        [TestCase("20mer1.fasta")]
        [TestCase("20mer1.fasta.gz")]
        [TestCase("20mer1.xml")]
        [TestCase("20mer1.xml.gz")]
        public static void TestDbReadingDifferentExtensions(string databaseFileName)
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData",
                databaseFileName);

            List<RNA> rna;
            if (dbPath.Contains("fasta"))
                rna = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false,
                    out var errors);
            else
                rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false,
                    new List<Modification>(), new List<string>(), out _);
            
            Assert.That(rna.Count, Is.EqualTo(1));
            Assert.That(rna.First().BaseSequence, Is.EqualTo("GUACUGCCUCUAGUGAAGCA"));
        }

        [Test]
        public static void TestEnsemblFastaParsing()
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TestDatabase_Ensembl.GRCh38.ncrna.fa");
            var oligos = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false, out var errors);

            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.EqualTo(6));

            var first = oligos.First();
            Assert.That(first.Accession, Is.EqualTo("ENST00000616830.1"));
            Assert.That(first.Name, Is.EqualTo("GRCh38:KI270744.1:51009:51114:-1"));
            Assert.That(first.BaseSequence, Is.EqualTo("GUACUUAUUUCAACAGCACAUAUUUUAAAUUGGAUCAAUACAGAGCAGAUAAGCAUGGUU" +
                "ACUGCCUAGGGAUGGCACACAAAUUCAGAAAGCAUUCCAUAUUUUG"));
            Assert.That(first.Organism, Is.EqualTo("GRCh38"));

            // GeneNames: should contain the gene ENSG00000278625.1
            Assert.That(first.GeneNames.Count, Is.EqualTo(1));
            Assert.That(first.GeneNames.First().Item1, Is.EqualTo("ENSG00000278625.1"));

            // Additional fields
            var fields = first.AdditionalDatabaseFields!;
            Assert.That(fields["GeneBiotype"], Is.EqualTo("snRNA"));
            Assert.That(fields["TranscriptBiotype"], Is.EqualTo("snRNA"));
            Assert.That(fields["GeneSymbol"], Is.EqualTo("U6"));
            Assert.That(fields["Description"], Does.Contain("U6 spliceosomal RNA"));

        }

        [Test]
        public static void TestNcbiAssemblyFastaParsing()
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TestDatabase_NcbiAssembly.fna");
            var oligos = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false, out var errors);

            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.GreaterThanOrEqualTo(1));

            var first = oligos.First();
            Assert.That(first.Accession, Is.EqualTo("NM_000014.6"));
            Assert.That(first.Name, Is.EqualTo("Homo sapiens alpha-2-macroglobulin "));
            Assert.That(first.BaseSequence, Does.StartWith("GGGACCAGAUGGAUUGUAGGGAGUAGGGUACAAUACAGUCUGUUCUCCUCCAGCUCCUUCUUUCUGCAACAUGGGGAAGA"));
            Assert.That(first.Organism, Is.EqualTo("Homo sapiens"));
            Assert.That(first.GeneNames.Count, Is.EqualTo(1));
            Assert.That(first.GeneNames.First().Item1, Is.EqualTo("alpha-2-macroglobulin (A2M)"));
        }

        [Test]
        public static void TestNcbiRefSeqGeneFastaParsing()
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TestDatabase_NcbiRefSeq_gene.fna");
            var oligos = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false, out var errors);

            Assert.That(errors.Count, Is.EqualTo(0));
            Assert.That(oligos.Count, Is.GreaterThanOrEqualTo(1));

            var first = oligos.First();
            Assert.That(first.Accession, Is.EqualTo("NC_051336.1"));
            Assert.That(first.Name, Is.EqualTo("Muc2 "));
            Assert.That(first.BaseSequence, Does.StartWith("GCUCUUCUGUGCCACCCUCGUGAGCCACCAUGGGGCUGCCACUAGCUCGCCUGGUGGCUGUGUGCCUAGU"));
            Assert.That(first.Organism, Is.EqualTo("Rattus norvegicus"));
            Assert.That(first.GeneNames.Count, Is.EqualTo(1));
            Assert.That(first.GeneNames.First().Item1, Is.EqualTo("24572"));
            Assert.That(first.AdditionalDatabaseFields!["Chromosome"], Is.EqualTo("1"));
        }

    }
}
