using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

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

            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData");
            Directory.CreateDirectory(outDir);
            var outpath = Path.Combine(outDir, $"ModomicsUnmodifiedTrimmed_{Guid.NewGuid():N}.xml");

            try
            {
                var xml = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), rna, outpath);

                var temp = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.None, false,
                    new List<Modification>() { methylG }, new List<string>(), out var unknownMods);

                Assert.That(unknownMods.Count, Is.EqualTo(0));
                Assert.That(temp.Count, Is.EqualTo(5));

                // Select the modified entry explicitly (accession SO:0000254), not by list order
                var modified = temp.FirstOrDefault(t => string.Equals(t.Accession, "SO:0000254", StringComparison.Ordinal))
                               ?? temp.FirstOrDefault(t => t.OneBasedPossibleLocalizedModifications?.Count == 2);
                Assert.That(modified, Is.Not.Null, "Modified RNA entry not found after round-trip.");

                var loadedMods = modified!.OneBasedPossibleLocalizedModifications;
                Assert.That(loadedMods.Count, Is.EqualTo(2));
                Assert.That(loadedMods[3].Count, Is.EqualTo(1));
                Assert.That(loadedMods[4].Count, Is.EqualTo(1));
                Assert.That(loadedMods[3].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
                Assert.That(loadedMods[4].First().IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
            }
            finally
            {
                try { if (File.Exists(outpath)) File.Delete(outpath); } catch { /* ignore cleanup errors */ }
            }
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

            var xml = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), rna, outpath);
            var temp = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.None, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods);

            Assert.That(unknownMods.Count, Is.EqualTo(0));
            Assert.That(temp.Count, Is.EqualTo(5));

            // Select modified entry explicitly
            var modified = temp.FirstOrDefault(t => string.Equals(t.Accession, "SO:0000254", StringComparison.Ordinal))
                           ?? temp.FirstOrDefault(t => t.OneBasedPossibleLocalizedModifications?.Count == 2);
            Assert.That(modified, Is.Not.Null, "Modified RNA entry not found after round-trip.");

            var loadedMods = modified!.OneBasedPossibleLocalizedModifications;
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
            var xml = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), rna, outpath);
            rna = RnaDbLoader.LoadRnaXML(outpath, true, DecoyType.Reverse, false,
                new List<Modification>() { methylG }, new List<string>(), out var unknownMods, decoyIdentifier: "rev");

            foreach (var oligo in rna)
            {
                if (!oligo.IsDecoy) continue;

                Assert.That(oligo.Accession, Does.StartWith("rev"));
                Assert.That(oligo.Accession, Does.Not.StartWith("DECOY"));
            }
        }

        // Helper to compute expected transcription for long inputs
        private static string ExpectedTranscription(string dna, bool isCodingStrand)
        {
            if (isCodingStrand)
            {
                // Coding strand: replace T with U
                return dna.Replace('T', 'U');
            }

            // Template strand: nucleotide complement with RNA bases (A->U, T->A, C->G, G->C)
            var sb = new StringBuilder(dna.Length);
            foreach (char c in dna)
            {
                sb.Append(c switch
                {
                    'A' => 'U',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    _   => c
                });
            }
            return sb.ToString();
        }

        [Test]
        public static void TestTranscribe_Long_Coding()
        {
            var input =
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT";
            var expected = ExpectedTranscription(input, true);
            Assert.That(input.Transcribe(true), Is.EqualTo(expected));
        }

        [Test]
        public static void TestTranscribe_Long_Template()
        {
            var input =
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT" +
                "ATCGACGAATCACGATCAGTCATGCATTGCTAACT";
            var expected = ExpectedTranscription(input, false);
            Assert.That(input.Transcribe(false), Is.EqualTo(expected));
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
        [Test]
        public static void TestLoadRnaXmlWithSequenceVariation_ExpandsAppliedVariants()
        {
            // Create a simple RNA with one sequence variant: position 3 G->A
            // Canonical: ACGUACGU  -> Variant: ACAUACGU
            var seq = "ACGUACGU";
            var variants = new List<SequenceVariation>
            {
                new SequenceVariation(
                    oneBasedPosition: 3,
                    originalSequence: "G",
                    variantSequence: "A",
                    description: "SNP:G3A")
            };

            var rnaWithVar = new RNA(
                sequence: seq,
                accession: "TEST-RNA-1",
                oneBasedPossibleModifications: null,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "Test RNA with 1 variant",
                organism: "UnitTestus",
                databaseFilePath: null,
                isContaminant: false,
                isDecoy: false,
                geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("primary", "GENE1") },
                databaseAdditionalFields: null,
                truncationProducts: null,
                sequenceVariations: variants,
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                fullName: "Test RNA with 1 variant (full)");

            // Write to a temporary XML under test data folder
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData");
            Directory.CreateDirectory(outDir);
            var outPath = Path.Combine(outDir, "RnaWithSeqVar.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<RNA> { rnaWithVar }, outPath);

            // Load with variant expansion enabled:
            var loaded = RnaDbLoader.LoadRnaXML(
                rnaDbLocation: outPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                isContaminant: false,
                allKnownModifications: Array.Empty<Modification>(),
                modTypesToExclude: Array.Empty<string>(),
                unknownModifications: out var unknownMods,
                maxThreads: 1,
                maxSequenceVariantsPerIsoform: 1,
                minAlleleDepth: 0,
                maxSequenceVariantIsoforms: 2);

            Assert.That(unknownMods.Count, Is.EqualTo(0), "No unknown modifications expected.");
            Assert.That(loaded.Count, Is.GreaterThanOrEqualTo(2), "Expected canonical and at least one applied-variant RNA.");

            // Find canonical (same accession, no applied variants)
            var canonical = loaded.FirstOrDefault(r =>
                r.Accession == "TEST-RNA-1" &&
                (r.AppliedSequenceVariations == null || r.AppliedSequenceVariations.Count == 0));

            // Find applied (has applied variants; accession starts with canonical accession + variant tag)
            var applied = loaded.FirstOrDefault(r =>
                r.AppliedSequenceVariations != null &&
                r.AppliedSequenceVariations.Count > 0 &&
                r.Accession.StartsWith("TEST-RNA-1", StringComparison.Ordinal));

            Assert.That(canonical, Is.Not.Null, "Canonical RNA should be present.");
            Assert.That(applied, Is.Not.Null, "Applied-variant RNA should be present.");

            // Canonical assertions
            Assert.That(canonical!.Accession, Is.EqualTo("TEST-RNA-1"));
            Assert.That(canonical.BaseSequence, Is.EqualTo(seq), "Canonical base sequence should match input.");
            Assert.That(canonical.SequenceVariations, Is.Not.Null);
            Assert.That(canonical.SequenceVariations.Count, Is.EqualTo(1), "Canonical should carry the candidate variant annotation.");

            var cv = canonical.SequenceVariations[0];
            Assert.That(cv.OneBasedBeginPosition, Is.EqualTo(3));
            Assert.That(cv.OneBasedEndPosition, Is.EqualTo(3));
            Assert.That(cv.OriginalSequence, Is.EqualTo("G"));
            Assert.That(cv.VariantSequence, Is.EqualTo("A"));

            // Applied variant assertions
            // The variant-applied base sequence must reflect G(3)->A substitution
            Assert.That(applied!.BaseSequence, Is.EqualTo("ACAUACGU"), "Applied variant base sequence should be mutated at position 3.");
            Assert.That(applied.Accession, Does.StartWith("TEST-RNA-1"), "Applied accession should be based on the canonical accession.");
            Assert.That(applied.Accession, Does.Contain("_"), "Applied accession should include a variant tag suffix.");

            // This test did not add any variant-specific modifications; ensure none exist
            Assert.That(applied.OneBasedPossibleLocalizedModifications == null
                        || applied.OneBasedPossibleLocalizedModifications.Count == 0,
                        Is.True, "No base-level modifications expected in this test.");
        }
        [Test]
        public static void TestLoadRnaXmlWithSequenceVariation_CanonicalOnlyByDefault()
        {
            // Reuse the same XML as previous test to avoid duplication
            var outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "RnaWithSeqVar.xml");
            Assert.That(File.Exists(outPath), "Expected RnaWithSeqVar.xml to exist from prior test.");

            // Load with default variant parameters:
            // Defaults are maxSequenceVariantsPerIsoform = 0 and maxSequenceVariantIsoforms = 1,
            // which should produce only the canonical entry (no variant-applied isoforms).
            var loaded = RnaDbLoader.LoadRnaXML(
                rnaDbLocation: outPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                isContaminant: false,
                allKnownModifications: Array.Empty<Modification>(),
                modTypesToExclude: Array.Empty<string>(),
                unknownModifications: out var unknownMods);

            Assert.That(unknownMods.Count, Is.EqualTo(0), "No unknown modifications expected.");

            // Expect exactly one entry (canonical only)
            Assert.That(loaded.Count, Is.EqualTo(1), "Default parameters should not emit applied-variant isoforms.");

            var canonical = loaded[0];
            Assert.That(canonical.Accession, Is.EqualTo("TEST-RNA-1"));
            Assert.That(canonical.BaseSequence, Is.EqualTo("ACGUACGU"));

            // The candidate variant should be present on the canonical entry as an annotation
            Assert.That(canonical.SequenceVariations, Is.Not.Null);
            Assert.That(canonical.SequenceVariations.Count, Is.EqualTo(1));
            Assert.That(canonical.AppliedSequenceVariations == null || canonical.AppliedSequenceVariations.Count == 0, Is.True,
                "No applied variants expected under default parameters.");
        }
        [Test]
        public static void TestVariantSpecificModification_PromotedAndPersistsThroughXml()
        {
            // Create a variant-specific modification (targets G)
            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   G\r\nCF   C1H2\r\n//";
            var methylG = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> _).First();

            // Canonical RNA has no base (consensus) modifications, but it has 1 candidate sequence variation:
            // Position 2: A -> G, with a variant-specific methylG at absolute position 2 (post-variation coordinate system)
            var canonicalSeq = "AACU";
            var variantPosition = 2;
            var svMods = new Dictionary<int, List<Modification>> { [variantPosition] = new List<Modification> { methylG } };
            var seqVar = new SequenceVariation(
                oneBasedPosition: variantPosition,
                originalSequence: "A",
                variantSequence: "G",
                description: "A2G with methylG",
                variantCallFormatDataString: null,
                oneBasedModifications: svMods);

            var rnaCanonical = new RNA(
                sequence: canonicalSeq,
                accession: "TEST-RNA-2",
                oneBasedPossibleModifications: null,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "ConsRNA_NoBaseMods_OneVariantWithMod",
                organism: "UnitTestus",
                databaseFilePath: null,
                isContaminant: false,
                isDecoy: false,
                geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("primary", "GENE2") },
                databaseAdditionalFields: null,
                truncationProducts: null,
                sequenceVariations: new List<SequenceVariation> { seqVar },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                fullName: "Consensus RNA with variant-specific mod");

            // Write canonical to XML
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData");
            Directory.CreateDirectory(outDir);
            var xmlPath = Path.Combine(outDir, "RnaVarWithVariantMod.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<RNA> { rnaCanonical }, xmlPath);

            // Load with variant expansion enabled to generate an applied-variant RNA
            var loaded = RnaDbLoader.LoadRnaXML(
                rnaDbLocation: xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                isContaminant: false,
                allKnownModifications: new List<Modification> { methylG },
                modTypesToExclude: Array.Empty<string>(),
                unknownModifications: out var unknownMods,
                maxThreads: 1,
                maxSequenceVariantsPerIsoform: 1, // allow applying the variant
                minAlleleDepth: 0,
                maxSequenceVariantIsoforms: 2);  // emit canonical + applied-variant

            Assert.That(unknownMods.Count, Is.EqualTo(0), "No unknown modifications expected.");
            Assert.That(loaded.Count, Is.GreaterThanOrEqualTo(2), "Expected canonical and applied-variant RNAs.");

            // Find canonical (same accession, no applied variants)
            var canonical = loaded.FirstOrDefault(r =>
                r.Accession == "TEST-RNA-2" &&
                (r.AppliedSequenceVariations == null || r.AppliedSequenceVariations.Count == 0));

            // Find applied (has applied variants; accession is prefixed by the canonical accession + variant tag)
            var applied = loaded.FirstOrDefault(r =>
                r.AppliedSequenceVariations != null &&
                r.AppliedSequenceVariations.Count > 0 &&
                r.Accession.StartsWith("TEST-RNA-2", StringComparison.Ordinal));

            Assert.That(canonical, Is.Not.Null, "Canonical RNA should be present.");
            Assert.That(applied, Is.Not.Null, "Applied-variant RNA should be present.");

            // Canonical assertions
            Assert.That(canonical!.BaseSequence, Is.EqualTo(canonicalSeq));
            Assert.That(canonical.OneBasedPossibleLocalizedModifications == null || canonical.OneBasedPossibleLocalizedModifications.Count == 0, Is.True);

            // Applied assertions...
            var expectedAppliedSeq = "AGCU";
            Assert.That(applied!.BaseSequence, Is.EqualTo(expectedAppliedSeq), "Applied variant base sequence should reflect A2G at position 2.");
            // Accessions for applied variants should include a variant suffix (e.g., "_A2G")
            Assert.That(applied.Accession, Does.StartWith("TEST-RNA-2"), "Applied accession should be based on the canonical accession.");
            Assert.That(applied.Accession, Does.Contain("_"), "Applied accession should include a variant tag suffix.");
            Assert.That(applied.OneBasedPossibleLocalizedModifications, Is.Not.Null);
            Assert.That(applied.OneBasedPossibleLocalizedModifications.ContainsKey(variantPosition), Is.True, "Variant mod should be promoted to RNA at pos 2.");
            Assert.That(applied.OneBasedPossibleLocalizedModifications[variantPosition].Count, Is.EqualTo(1));
            Assert.That(applied.OneBasedPossibleLocalizedModifications[variantPosition][0].IdWithMotif, Is.EqualTo(methylG.IdWithMotif));

            // Now write ONLY the applied-variant RNA back to XML and re-load to ensure the mod persists through IO
            var appliedOnlyPath = Path.Combine(outDir, "RnaVarWithVariantMod_AppliedOnly.xml");
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                new List<RNA> { applied },
                appliedOnlyPath,
                includeAppliedVariantEntries: true); // write applied variant entries, too

            var roundtrip = RnaDbLoader.LoadRnaXML(
                rnaDbLocation: appliedOnlyPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                isContaminant: false,
                allKnownModifications: new List<Modification> { methylG },
                modTypesToExclude: Array.Empty<string>(),
                unknownModifications: out var unknown2);

            Assert.That(unknown2.Count, Is.EqualTo(0), "Roundtrip: no unknown modifications expected.");
            Assert.That(roundtrip.Count, Is.GreaterThanOrEqualTo(1), "Roundtrip should load at least one entry.");

            // Find the applied isoform we wrote (accession prefix + mutated sequence)
            var rt = roundtrip.FirstOrDefault(r =>
                r.Accession.StartsWith("TEST-RNA-2", StringComparison.Ordinal) &&
                r.BaseSequence == expectedAppliedSeq);

            Assert.That(rt, Is.Not.Null, "Roundtrip applied-variant RNA not found.");

            // The roundtrip RNA should keep the applied sequence and the promoted modification
            Assert.That(rt!.BaseSequence, Is.EqualTo(expectedAppliedSeq), "Roundtrip base sequence should match applied variant.");
            Assert.That(rt.OneBasedPossibleLocalizedModifications, Is.Not.Null);
            Assert.That(rt.OneBasedPossibleLocalizedModifications.ContainsKey(variantPosition), Is.True);
            Assert.That(rt.OneBasedPossibleLocalizedModifications[variantPosition].Count, Is.EqualTo(1));
            Assert.That(rt.OneBasedPossibleLocalizedModifications[variantPosition][0].IdWithMotif, Is.EqualTo(methylG.IdWithMotif));
        }
        [Test]
        public static void TestTruncationVariant_RemovesDownstreamModification_PersistsThroughXml()
        {
            // Base sequence (length 13). We will delete positions 10..13 (truncate tail).
            var baseSeq = "GUACUGUAGCCUA";
            // Place a consensus modification at position 12 (this site will be removed by the truncation)
            var modString = "ID   Methylation\r\nMT   Biological\r\nPP   Anywhere.\r\nTG   U\r\nCF   C1H2\r\n//";
            var methylU = PtmListLoader.ReadModsFromString(modString, out List<(Modification, string)> _).First();

            var consensusMods = new Dictionary<int, List<Modification>>
            {
                [12] = new List<Modification> { methylU }
            };

            // Define a deletion variant: remove positions 10..13 (inclusive).
            // For correctness, set OriginalSequence to the actual substring being removed.
            int delBegin = 10, delEnd = 13;
            string originalSpan = baseSeq.Substring(delBegin - 1, delEnd - delBegin + 1);
            var truncation = new SequenceVariation(
                oneBasedPosition: delBegin,
                originalSequence: originalSpan,
                variantSequence: "",
                description: "deletion(10..13)");

            var canonical = new RNA(
                sequence: baseSeq,
                accession: "TRUNC-RNA-1",
                oneBasedPossibleModifications: consensusMods,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "TruncationTest",
                organism: "UnitTestus",
                databaseFilePath: null,
                isContaminant: false,
                isDecoy: false,
                geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("primary", "GENE-T") },
                databaseAdditionalFields: null,
                truncationProducts: null,
                sequenceVariations: new List<SequenceVariation> { truncation },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                fullName: "RNA with tail-deletion variant");

            // Expand to get applied variant isoform
            var isoforms = canonical.GetVariantBioPolymers(
                maxSequenceVariantsPerIsoform: 1,
                minAlleleDepth: 0,
                maxSequenceVariantIsoforms: 2);

            Assert.That(isoforms.Count, Is.GreaterThanOrEqualTo(2), "Expected canonical + applied variant.");

            var applied = isoforms.FirstOrDefault(r => r.AppliedSequenceVariations.Count > 0);
            var refLike = isoforms.FirstOrDefault(r => r.AppliedSequenceVariations.Count == 0);

            Assert.That(applied, Is.Not.Null, "Applied truncation isoform not found.");
            Assert.That(refLike, Is.Not.Null, "Canonical isoform not found.");

            // Expected applied sequence (remove 10..13)
            var expectedAppliedSeq = baseSeq.Substring(0, delBegin - 1);
            Assert.That(applied!.BaseSequence, Is.EqualTo(expectedAppliedSeq), "Applied sequence should be truncated.");

            // Precondition: consensus has the mod at position 12
            Assert.That(refLike!.OneBasedPossibleLocalizedModifications.ContainsKey(12), Is.True,
                "Consensus should have a modification at position 12.");

            // After truncation, mod at 12 must be gone (position out of range)
            Assert.That(applied.OneBasedPossibleLocalizedModifications.ContainsKey(12), Is.False,
                "Applied truncation isoform should not retain a modification at removed position 12.");

            // Also ensure no modification key exceeds applied length
            int appliedLen = applied.Length;
            Assert.That(applied.OneBasedPossibleLocalizedModifications.Keys.All(k => k >= 1 && k <= appliedLen), Is.True,
                "Applied isoform contains a modification indexed outside its new length.");

            // Roundtrip: write consensus + applied, including applied entries, then reload
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData");
            Directory.CreateDirectory(outDir);
            var outPath = Path.Combine(outDir, $"TruncVar_{Guid.NewGuid():N}.xml");

            try
            {
                ProteinDbWriter.WriteXmlDatabase(
                    new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    new List<RNA> { canonical, applied },
                    outPath,
                    includeAppliedVariantEntries: true);

                var reloaded = RnaDbLoader.LoadRnaXML(
                    rnaDbLocation: outPath,
                    generateTargets: true,
                    decoyType: DecoyType.None,
                    isContaminant: false,
                    allKnownModifications: new List<Modification> { methylU },
                    modTypesToExclude: Array.Empty<string>(),
                    unknownModifications: out var unknownMods);

                Assert.That(unknownMods.Count, Is.EqualTo(0), "No unknown mods expected on reload.");
                Assert.That(reloaded.Count, Is.GreaterThanOrEqualTo(2), "Reloaded set should contain canonical and applied.");

                var reApplied = reloaded.FirstOrDefault(r =>
                    r.Accession.StartsWith("TRUNC-RNA-1", StringComparison.Ordinal) &&
                    string.Equals(r.BaseSequence, expectedAppliedSeq, StringComparison.Ordinal));

                var reCanon = reloaded.FirstOrDefault(r =>
                    r.Accession == "TRUNC-RNA-1" &&
                    (r.AppliedSequenceVariations == null || r.AppliedSequenceVariations.Count == 0));

                Assert.That(reApplied, Is.Not.Null, "Reloaded applied truncation isoform not found.");
                Assert.That(reCanon, Is.Not.Null, "Reloaded canonical isoform not found.");

                // Verify applied is still truncated and lacks the removed-site modification
                Assert.That(reApplied!.BaseSequence, Is.EqualTo(expectedAppliedSeq));
                Assert.That(reApplied.OneBasedPossibleLocalizedModifications.ContainsKey(12), Is.False,
                    "Reloaded applied truncation isoform should not have mod at removed position 12.");

                // Verify canonical retains the original site modification
                Assert.That(reCanon!.OneBasedPossibleLocalizedModifications.ContainsKey(12), Is.True,
                    "Reloaded canonical should retain the mod at position 12.");
                Assert.That(reCanon.OneBasedPossibleLocalizedModifications[12][0].IdWithMotif, Is.EqualTo(methylU.IdWithMotif));
            }
            finally
            {
                try { if (File.Exists(outPath)) File.Delete(outPath); } catch { /* ignore */ }
            }
        }
    }
}
