using MassSpectrometry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.BioPolymer;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestProteomicsReadWrite
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public void ReadXmlNulls()
        {
            var ok = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"),
                true, DecoyType.None, null, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
        }
        [Test]
        public void ReadSomeOldXmlWithLongSubstitutionThatHasAConflict()
        {
            //In this case, we have two different sequence variants. One is a long substitution, the other is a point mutation.
            //If their positions didn't overlap, we should end up with four total protein sequences: the base protein, the protein with the long substitution,
            //the protein with the point mutation, and the protein with both the long substitution and the point mutation.
            //but, because the point mutation falls within the range of the long substitution, we should only end up with three total protein sequences:
            string oldXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"longSubstitution.xml");
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(oldXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 2,
                maxSequenceVariantIsoforms: 100);
            Assert.IsTrue(ok.Count == 3);
        }
        [Test]
        public void SequenceVariantRefersToAlternateIsoform()
        {
            //In this case, we have a sequence variant that refers to an alternate isoform.
            //We should still be able to load the protein, even if we don't have the alternate isoform sequence.
            //for now we are ignoring the sequence variant if we don't have the alternate isoform sequence.
            string oldXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"sequenceVariantOnAlternateIsoform.xml");
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(
                oldXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            Assert.IsTrue(ok.Count == 1);
        }
        [Test]
        public void ReadXmlSkipVariants()
        {
            //In this case, we have a couple different sequence variants. But, we don't want to apply any of them.
            //instead, we just want the base protein sequence with mods.
            string oldXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"longSubstitution.xml");
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(oldXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un, maxSequenceVariantIsoforms: 1);
            Assert.IsTrue(ok.Count == 1);
        }
        [Test]
        public void Test_readUniProtXML_writeProteinXml()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null)
            };

            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"),
                true, DecoyType.None, uniprotPtms.Concat(nice), false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            // Write and read back
            string outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, outPath);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(
                outPath, true, DecoyType.None, nice, false, new List<string>(),
                out un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            // Count equality
            Assert.AreEqual(ok.Count, ok2.Count);

            // Compare order-independently by accession
            var byAcc1 = ok.ToDictionary(p => p.Accession, p => p);
            var byAcc2 = ok2.ToDictionary(p => p.Accession, p => p);

            CollectionAssert.AreEquivalent(byAcc1.Keys, byAcc2.Keys);

            foreach (var acc in byAcc1.Keys)
            {
                // Base sequence round-trip
                Assert.AreEqual(byAcc1[acc].BaseSequence, byAcc2[acc].BaseSequence, $"BaseSequence mismatch for {acc}");

                // Gene name (first)
                var g1 = byAcc1[acc].GeneNames.First().Item2;
                var g2 = byAcc2[acc].GeneNames.First().Item2;
                Assert.AreEqual(g1, g2, $"Gene name mismatch for {acc}");

                // Full name
                Assert.AreEqual(byAcc1[acc].FullName, byAcc2[acc].FullName, $"FullName mismatch for {acc}");
            }

            // Keep detailed checks but anchor them to the same protein as ok[0]
            var anchorAcc = ok[0].Accession;

            Assert.AreEqual(9, byAcc1[anchorAcc].DatabaseReferences.Count(dbRef => dbRef.Type == "GO"));
            Assert.AreEqual(1, byAcc1[anchorAcc].DatabaseReferences.Count(dbRef => dbRef.Type == "GeneID"));
            Assert.AreEqual(3, byAcc1[anchorAcc].DatabaseReferences.First(dbRef => dbRef.Type == "GO").Properties.Count());
            Assert.AreEqual(3, byAcc1[anchorAcc].GeneNames.Count());
            Assert.AreEqual("primary", byAcc1[anchorAcc].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", byAcc1[anchorAcc].GeneNames.First().Item2);
            Assert.AreEqual("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", byAcc1[anchorAcc].Organism);
            Assert.AreEqual(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), byAcc1[anchorAcc].DatabaseFilePath);

            Assert.AreEqual(9, byAcc2[anchorAcc].DatabaseReferences.Count(dbRef => dbRef.Type == "GO"));
            Assert.AreEqual(3, byAcc2[anchorAcc].DatabaseReferences.First(dbRef => dbRef.Type == "GO").Properties.Count());
            Assert.AreEqual(3, byAcc2[anchorAcc].GeneNames.Count());
            Assert.AreEqual("primary", byAcc2[anchorAcc].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", byAcc2[anchorAcc].GeneNames.First().Item2);
            Assert.AreEqual("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", byAcc2[anchorAcc].Organism);
            Assert.AreEqual(outPath, byAcc2[anchorAcc].DatabaseFilePath);

            // Truncation product bounds remain valid
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }
        [Test]
        public void Test_readUniProtXML_writeProteinXmlCheckEntryUpdated()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            string inputXmlPath = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml");
            bool lineModified = false;
            foreach (var line in File.ReadLines(inputXmlPath))
            {
                if (line.Contains("<entry"))
                {
                    string modified = "modified=\"" + "2017-02-15";
                    // Checks if the line contains the exact entry text
                    lineModified = line.Contains(modified);
                    
                    break;
                }
            }
            Assert.IsTrue(lineModified);
            lineModified = false; // Reset for the next check
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(
                inputXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml");

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, outputPath, true);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"),
                true, DecoyType.None, uniprotPtms, false, new List<string>(),
                out un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            foreach (var line in File.ReadLines(outputPath))
            {
                if (line.Contains("<entry"))
                {
                    //make sure that when the xml is written, the modified date is updated to today's date
                    string todaysDate = DateTime.Now.ToString("yyyy-MM-dd");
                    string modified = "modified=\"" + todaysDate;
                    // Checks if the line contains the exact entry text
                    lineModified = line.Contains(modified);
                    
                }
            }
            Assert.IsTrue(lineModified);
            // Clean up the output file after the test
            if (File.Exists(outputPath))
            {
                File.Delete(outputPath);
            }
        }

        [Test]
        public void Test_readUniProtXML_featureBeginEndPosition()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            string inputXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"unknownStatus.xml");
            string targetLine = "<begin status=\"unknown\"/>";
            bool foundTargetLine = false;
            foreach (var line in File.ReadLines(inputXmlPath))
            {
                if (line.Equals(targetLine))
                {
                    foundTargetLine = true;
                    Assert.IsTrue(foundTargetLine); //"Found the target line with unknown status for begin position."
                    foundTargetLine = false; // Reset for the next check
                    break;
                }
            }

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(
                inputXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_unknownStatus.xml");

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, outputPath, true);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(
                outputPath, true, DecoyType.None, uniprotPtms, false, new List<string>(),
                out un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            foreach (var line in File.ReadLines(outputPath))
            {
                if (line.Equals(targetLine))
                {
                    foundTargetLine = true;
                    Assert.IsTrue(foundTargetLine); //"Found the target line with unknown status for begin position."
                    break;
                }
            }
            // Clean up the output file after the test
            if (File.Exists(outputPath))
            {
                File.Delete(outputPath);
            }
        }
        [Test]
        public void Test_read_Ensembl_pepAllFasta()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            string fastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta");
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml");

            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(
                fastaPath, true, DecoyType.None, false, out var a,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null);

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, xmlPath);

            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(
                xmlPath, true, DecoyType.None, nice, false, null, out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            // Counts equal
            Assert.AreEqual(ok.Count, ok2.Count);

            // Compare by accession (order-independent)
            var okByAcc = ok.ToDictionary(p => p.Accession, p => p);
            var ok2ByAcc = ok2.ToDictionary(p => p.Accession, p => p);
            CollectionAssert.AreEquivalent(okByAcc.Keys, ok2ByAcc.Keys);

            // Validate per-accession equality for sequence, gene name (first), and full name
            foreach (var acc in okByAcc.Keys)
            {
                Assert.AreEqual(okByAcc[acc].BaseSequence, ok2ByAcc[acc].BaseSequence, $"BaseSequence mismatch for {acc}");

                var okGene = okByAcc[acc].GeneNames.First().Item2;
                var ok2Gene = ok2ByAcc[acc].GeneNames.First().Item2;
                Assert.AreEqual(okGene, ok2Gene, $"Gene name mismatch for {acc}");

                Assert.AreEqual(okByAcc[acc].FullName, ok2ByAcc[acc].FullName, $"FullName mismatch for {acc}");
            }

            // Explicit content checks (still order-independent)
            var expectedAccs = new[] { "ENSP00000381386", "ENSP00000215773" };
            CollectionAssert.IsSubsetOf(expectedAccs, okByAcc.Keys);
            CollectionAssert.IsSubsetOf(expectedAccs, ok2ByAcc.Keys);

            Assert.AreEqual("ENSG00000099977", okByAcc["ENSP00000381386"].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", okByAcc["ENSP00000215773"].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok2ByAcc["ENSP00000381386"].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok2ByAcc["ENSP00000215773"].GeneNames.First().Item2);

            Assert.AreEqual(
                "pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding",
                okByAcc["ENSP00000381386"].FullName);
            Assert.AreEqual(
                "pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding",
                okByAcc["ENSP00000215773"].FullName);
            Assert.AreEqual(
                "pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding",
                ok2ByAcc["ENSP00000381386"].FullName);
            Assert.AreEqual(
                "pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding",
                ok2ByAcc["ENSP00000215773"].FullName);

            // File paths (apply to all entries rather than a single index)
            Assert.True(ok.All(p => p.DatabaseFilePath == fastaPath));
            Assert.True(ok2.All(p => p.DatabaseFilePath == xmlPath));

            // Truncation product bounds remain valid
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }
        [Test]
        public static void FastaTest()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fasta.fasta"), true, DecoyType.Reverse, false, out var a,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);
            ProteinDbWriter.WriteFastaDatabase(prots, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_fasta.fasta"), "|");
            List<Protein> prots2 = ProteinDbLoader.LoadProteinFasta(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_fasta.fasta"),
                true,
                DecoyType.None,
                false,
                out var un,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);

            Assert.AreEqual("P62805", prots.First().Accession);
            Assert.AreEqual("H4_HUMAN", prots.First().Name);
            Assert.AreEqual("Histone H4", prots.First().FullName);
            Assert.AreEqual("HIST1H4A", prots.First().GeneNames.First().Item2);
            Assert.AreEqual("Homo sapiens", prots.First().Organism);

            Assert.AreEqual("P62805", prots2.First().Accession);
            Assert.AreEqual("H4_HUMAN", prots2.First().Name);
            Assert.AreEqual("Histone H4", prots2.First().FullName);
            Assert.AreEqual("HIST1H4A", prots2.First().GeneNames.First().Item2);
            Assert.AreEqual("Homo sapiens", prots2.First().Organism);
        }

        [Test]
        public void Test_read_write_read_fasta()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false, out var a,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), true, DecoyType.None, false, out var b,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));

            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void AddModsDirectlyToProteinDbWriter()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m = new Modification("mod", null, "mt", null, motif, "Anywhere.", null, 1, null, null, null, new Dictionary<DissociationType, List<double>>() { { DissociationType.AnyActivationType, new List<double> { -1 } } }, null, null);
            Dictionary<string, HashSet<Tuple<int, Modification>>> new_mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>
            {
                {  "P62805", new HashSet<Tuple<int, Modification>> {new Tuple<int, Modification>(6, m ) } }
            };
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fasta.fasta"), true, DecoyType.None, false, out var a,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(new_mods, ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_fasta.xml"));
            Assert.AreEqual(1, newModResEntries.Count);
            var key = newModResEntries.Keys.First();
            var value = newModResEntries[key];
            Assert.AreEqual("P62805", new_mods.Keys.First());
            Assert.AreEqual("mod on K", key);
            Assert.AreEqual(1, value);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_fasta.xml"), true, DecoyType.None,
                new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(0, ok[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, ok2[0].OneBasedPossibleLocalizedModifications.Count);
        }

        [Test]
        public void Test_read_xml_write_read_fasta()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"),
                true, DecoyType.None, nice, false, null,
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml_test.fasta"), "|");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml_test.fasta"), true, DecoyType.None, false, out var b,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].Name == ok2[i].Name));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].Organism == ok2[i].Organism));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].GeneNames.First().Item2 == ok2[i].GeneNames.First().Item2));

            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.TruncationProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void Test_accession_regex_weird()
        {
            FastaHeaderFieldRegex bad = new FastaHeaderFieldRegex("", @"/()/", 0, 1);
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false, out var a,
                bad, bad, bad, bad, bad);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), true, DecoyType.None, false, out var b,
                bad, bad, bad, bad, bad);

            Assert.AreEqual("ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[0].Accession);
            Assert.AreEqual("ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[0].Accession);
            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
        }
        [Test]
        public void Test_write_with_custom_mods()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif m1);
            ModificationMotif.TryGetMotif("T", out ModificationMotif m2);
            ModificationMotif.TryGetMotif("X", out ModificationMotif motiff);

            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motiff, "Anywhere.", null, 10, null, null, null, null, null, null),
                new Modification("Phosphoserine", null, "mt", null, m1, "Anywhere.", null, 80, null, null, null, null, null, null),
                new Modification("Phosphothreonine", null, "mt", null,  m2, "Anywhere.", null, 80, null, null, null, null, null, null)
            };

            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            Modification m = new Modification("mod", null, "mt", null, motif, "Anywhere.", null, 1, null, null, null, new Dictionary<DissociationType, List<double>>() { { DissociationType.AnyActivationType, new List<double> { -1 } } }, null, null);

            Dictionary<string, HashSet<Tuple<int, Modification>>> new_mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>
            {
                {  "P53863", new HashSet<Tuple<int, Modification>> {new Tuple<int, Modification>(2, m ) } }
            };

            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            // Load, write, reload
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, uniprotPtms.Concat(nice), false, new List<string>(),
                out Dictionary<string, Modification> un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"));
            Assert.AreEqual(0, newModResEntries.Count);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out un,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 1);

            // Count equality
            Assert.AreEqual(ok.Count, ok2.Count);

            // Compare order-independently by accession
            var byAcc1 = ok.ToDictionary(p => p.Accession, p => p);
            var byAcc2 = ok2.ToDictionary(p => p.Accession, p => p);

            CollectionAssert.AreEquivalent(byAcc1.Keys, byAcc2.Keys);

            // Base sequences must match per accession
            foreach (var acc in byAcc1.Keys)
            {
                Assert.AreEqual(byAcc1[acc].BaseSequence, byAcc2[acc].BaseSequence, $"BaseSequence mismatch for {acc}");
            }

            // The original test expected 2 possible localized mods on ok[0]; anchor by that accession
            var anchorAcc = ok[0].Accession;
            Assert.AreEqual(2, byAcc1[anchorAcc].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(2, byAcc2[anchorAcc].OneBasedPossibleLocalizedModifications.Count);
        }
        [Test]
        public void SmallXml_VariantTokens_And_Lengths()
        {
            // Arrange
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "small.xml");

            // Load with single-variant expansion (base + each single variant)
            var proteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: Enumerable.Empty<Modification>(),
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var _,
                maxSequenceVariantsPerIsoform: 1,
                maxSequenceVariantIsoforms: 50);

            // Expect: 1 base + 6 single-variant proteoforms
            Assert.AreEqual(7, proteins.Count, "Unexpected proteoform count (expected base + 6 variants).");

            // Collect base (no underscore) and variant proteoforms (underscore suffix)
            var baseProteins = proteins.Where(p => !p.Accession.Contains('_')).ToList();
            Assert.AreEqual(1, baseProteins.Count, "Should have exactly one base (non-suffixed) accession.");
            var baseProt = baseProteins.Single();
            int baseLength = baseProt.Length;

            // Expected variant tokens (SimpleString forms)
            var expectedTokens = new HashSet<string>
            {
                "S70N",
                "S311L",
                "C337CS",
                "AHMPC369-373VHMPY",
                "H383R",
                "K428E"
            };

            // Pull variant proteoforms
            var variantProteins = proteins.Where(p => p.Accession.Contains('_')).ToList();
            Assert.AreEqual(expectedTokens.Count, variantProteins.Count, "Mismatch in variant isoform count.");

            // Map accession suffix to proteoform
            var tokenToProtein = new Dictionary<string, Protein>(StringComparer.Ordinal);
            foreach (var vp in variantProteins)
            {
                string suffix = vp.Accession[(vp.Accession.IndexOf('_') + 1)..];
                tokenToProtein[suffix] = vp;
            }

            // Ensure all expected tokens present
            foreach (var token in expectedTokens)
            {
                Assert.IsTrue(tokenToProtein.ContainsKey(token), $"Missing variant accession token {token}");
            }

            // Insertion variant (C337CS) should have length +1
            Assert.AreEqual(baseLength + 1, tokenToProtein["C337CS"].Length, "Insertion variant length incorrect.");

            // All other variants should retain base length
            foreach (var kv in tokenToProtein.Where(kv => kv.Key != "C337CS"))
            {
                Assert.AreEqual(baseLength, kv.Value.Length, $"Length mismatch for {kv.Key}");
            }

            // UniProtSequenceAttributes integrity (present and matching length if available)
            foreach (var p in proteins)
            {
                if (p.UniProtSequenceAttributes != null)
                {
                    Assert.AreEqual(p.Length, p.UniProtSequenceAttributes.Length,
                        $"UniProtSequenceAttributes.Length mismatch for {p.Accession}");
                }
            }

            // AppliedSequenceVariations: base has none; each variant exactly one applied
            Assert.IsTrue(baseProt.AppliedSequenceVariations == null || baseProt.AppliedSequenceVariations.Count == 0,
                "Base protein should have no applied sequence variations.");

            foreach (var kv in tokenToProtein)
            {
                var ap = kv.Value.AppliedSequenceVariations;
                Assert.IsNotNull(ap, $"AppliedSequenceVariations null for {kv.Key}");
                Assert.AreEqual(1, ap.Count, $"Expected exactly 1 applied variant for {kv.Key}");
                Assert.AreEqual(kv.Key, ap[0].SimpleString(), $"Applied variant token mismatch for {kv.Key}");
            }

            // Base protein should enumerate all 6 defined variants (original annotations)
            Assert.IsNotNull(baseProt.SequenceVariations, "Base SequenceVariations null.");
            Assert.AreEqual(6, baseProt.SequenceVariations.Count(), "Base protein should define 6 sequence variants.");
            var baseVariantTokens = new HashSet<string>(baseProt.SequenceVariations.Select(v => v.SimpleString()));
            foreach (var token in expectedTokens)
            {
                Assert.IsTrue(baseVariantTokens.Contains(token), $"Base variant list missing {token}");
            }

            // Variant name tagging (variant:token present in Name for variants)
            foreach (var kv in tokenToProtein)
            {
                string name = kv.Value.Name ?? "";
                Assert.IsTrue(name.Contains(kv.Key) || name.Contains("variant:"), $"Variant name missing token hint for {kv.Key}");
            }

            // Accession uniqueness
            Assert.AreEqual(proteins.Count, proteins.Select(p => p.Accession).Distinct().Count(), "Duplicate accessions detected.");

            // Sequence uniqueness sanity: at least insertion differs in length; substitutions differ in sequence
            var seqSet = new HashSet<string>(proteins.Select(p => p.BaseSequence));
            Assert.IsTrue(seqSet.Count >= 2, "Expected at least two distinct sequences (insertion must differ).");
            Assert.IsTrue(tokenToProtein["C337CS"].BaseSequence.Length == baseProt.BaseSequence.Length + 1,
                "Insertion sequence length delta not observed.");

            // No zero-length sequences
            Assert.IsFalse(proteins.Any(p => string.IsNullOrEmpty(p.BaseSequence)), "Found empty BaseSequence.");

            // Final safety: all applied variants' coordinates are within sequence bounds
            foreach (var vp in variantProteins)
            {
                foreach (var sv in vp.AppliedSequenceVariations)
                {
                    Assert.IsTrue(sv.OneBasedBeginPosition >= 1 && sv.OneBasedBeginPosition <= vp.Length,
                        $"Begin out of range in {vp.Accession}");
                    Assert.IsTrue(sv.OneBasedEndPosition >= sv.OneBasedBeginPosition && sv.OneBasedEndPosition <= vp.Length,
                        $"End out of range in {vp.Accession}");
                }
            }
        }
        [Test]
        public void SmallXml_TwoVariantCombinations()
        {
            // Arrange
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "small.xml");

            var proteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: Enumerable.Empty<Modification>(),
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var _,
                maxSequenceVariantsPerIsoform: 2,
                maxSequenceVariantIsoforms: 200);

            var baseProt = proteins.Single(p => !p.Accession.Contains('_'));
            int baseLength = baseProt.Length;

            // Explicit expected single variant tokens (SimpleString forms)
            var expectedSingles = new List<string>
            {
                "S70N",
                "S311L",
                "C337CS",
                "AHMPC369-373VHMPY",
                "H383R",
                "K428E"
            };
            Assert.AreEqual(6, expectedSingles.Count, "Expected 6 single variant tokens.");

            // Explicit expected pair tokens (canonical: lower coordinate variant first)
            var expectedPairTokensOrdered = new List<string>
            {
                "S70N_S311L",
                "S70N_C337CS",
                "S70N_AHMPC369-373VHMPY",
                "S70N_H383R",
                "S70N_K428E",
                "S311L_C337CS",
                "S311L_AHMPC369-373VHMPY",
                "S311L_H383R",
                "S311L_K428E",
                "C337CS_AHMPC369-373VHMPY",
                "C337CS_H383R",
                "C337CS_K428E",
                "AHMPC369-373VHMPY_H383R",
                "AHMPC369-373VHMPY_K428E",
                "H383R_K428E"
            };
            Assert.AreEqual(15, expectedPairTokensOrdered.Count, "Expected 15 two-variant combinations.");

            var expectedSinglesSet = new HashSet<string>(expectedSingles);
            var expectedPairsCanonical = new HashSet<string>(expectedPairTokensOrdered);

            // Helper: extract first coordinate for ordering
            int ExtractBegin(string token)
            {
                for (int i = 0; i < token.Length; i++)
                {
                    if (char.IsDigit(token[i]))
                    {
                        int j = i;
                        while (j < token.Length && char.IsDigit(token[j])) j++;
                        return int.Parse(token[i..j]);
                    }
                }
                return int.MaxValue;
            }

            string CanonicalPair(string a, string b)
            {
                var ordered = new[] { a, b }
                    .OrderBy(t => ExtractBegin(t))
                    .ThenBy(t => t, StringComparer.Ordinal)
                    .ToArray();
                return $"{ordered[0]}_{ordered[1]}";
            }

            // Expected total: 1 base + 6 singles + 15 pairs = 22
            int expectedTotal = 1 + expectedSinglesSet.Count + expectedPairsCanonical.Count;
            Assert.AreEqual(expectedTotal, proteins.Count, "Unexpected total proteoform count.");

            var singleIsoforms = proteins.Where(p => p.Accession.Contains('_') && p.AppliedSequenceVariations.Count == 1).ToList();
            var pairIsoforms = proteins.Where(p => p.AppliedSequenceVariations.Count == 2).ToList();

            Assert.AreEqual(expectedSinglesSet.Count, singleIsoforms.Count, "Mismatch in single-variant isoform count.");
            Assert.AreEqual(expectedPairsCanonical.Count, pairIsoforms.Count, "Mismatch in pair-variant isoform count.");

            // Validate singles
            foreach (var iso in singleIsoforms)
            {
                string suffix = iso.Accession[(iso.Accession.IndexOf('_') + 1)..];
                Assert.IsTrue(expectedSinglesSet.Contains(suffix), $"Unexpected single variant accession suffix {suffix}");
                Assert.AreEqual(1, iso.AppliedSequenceVariations.Count, "Single isoform must have exactly one applied variant.");
                Assert.AreEqual(suffix, iso.AppliedSequenceVariations[0].SimpleString(), $"Applied variant token mismatch for {suffix}");

                // Length rule: only insertion C337CS adds +1
                if (suffix == "C337CS")
                    Assert.AreEqual(baseLength + 1, iso.Length, "Insertion single variant length incorrect.");
                else
                    Assert.AreEqual(baseLength, iso.Length, $"Length mismatch for single {suffix}");

                if (iso.UniProtSequenceAttributes != null)
                    Assert.AreEqual(iso.Length, iso.UniProtSequenceAttributes.Length, $"Attribute length mismatch (single) {suffix}");
            }

            // Track coverage of pairs
            var seenPairs = new HashSet<string>();

            // Validate pairs (order-insensitive)
            foreach (var iso in pairIsoforms)
            {
                var appliedTokens = iso.AppliedSequenceVariations
                    .Select(v => v.SimpleString())
                    .ToList();
                Assert.AreEqual(2, appliedTokens.Count, $"Applied variant count mismatch for {iso.Accession}");

                string canonical = CanonicalPair(appliedTokens[0], appliedTokens[1]);
                seenPairs.Add(canonical);

                Assert.IsTrue(expectedPairsCanonical.Contains(canonical),
                    $"Unexpected pair combination canonical={canonical} accession={iso.Accession}");

                bool containsInsertion = appliedTokens.Contains("C337CS");
                int expectedLen = containsInsertion ? baseLength + 1 : baseLength;
                Assert.AreEqual(expectedLen, iso.Length, $"Length mismatch for pair {canonical}");

                if (iso.UniProtSequenceAttributes != null)
                    Assert.AreEqual(iso.Length, iso.UniProtSequenceAttributes.Length, $"Attribute length mismatch (pair) {canonical}");

                string name = iso.Name ?? "";
                foreach (var t in appliedTokens)
                {
                    Assert.IsTrue(name.Contains(t) || name.Contains("variant:"),
                        $"Variant name missing token {t} for pair {canonical}");
                }

                // Non-overlap guarantee (data has disjoint variants)
                var spans = iso.AppliedSequenceVariations
                    .Select(v => (v.OneBasedBeginPosition, v.OneBasedEndPosition))
                    .OrderBy(s => s.OneBasedBeginPosition)
                    .ToList();
                Assert.IsTrue(spans[0].OneBasedEndPosition < spans[1].OneBasedBeginPosition,
                    $"Unexpected coordinate overlap in pair {canonical}");
            }

            // Report any missing / extra pairs explicitly
            var missingPairs = expectedPairsCanonical.Except(seenPairs).ToList();
            var unexpectedPairs = seenPairs.Except(expectedPairsCanonical).ToList();

            Assert.IsTrue(missingPairs.Count == 0,
                "Missing expected pair tokens: " + string.Join(", ", missingPairs));
            Assert.IsTrue(unexpectedPairs.Count == 0,
                "Found unexpected pair tokens: " + string.Join(", ", unexpectedPairs));

            // Global accession uniqueness
            Assert.AreEqual(proteins.Count, proteins.Select(p => p.Accession).Distinct().Count(), "Duplicate accessions detected.");

            // Coordinate sanity
            foreach (var iso in proteins.Where(p => p.AppliedSequenceVariations.Any()))
            {
                foreach (var sv in iso.AppliedSequenceVariations)
                {
                    Assert.That(sv.OneBasedBeginPosition, Is.InRange(1, iso.Length),
                        $"Begin out of range ({sv.OneBasedBeginPosition}) in {iso.Accession}");
                    Assert.That(sv.OneBasedEndPosition, Is.InRange(sv.OneBasedBeginPosition, iso.Length),
                        $"End out of range ({sv.OneBasedEndPosition}) in {iso.Accession}");
                }
            }
        }

        //[Test]
        //[Explicit("Long-running diagnostic; generates protein_variant_log.txt with per-protein variant expansion results.")]
        //public void LargeXml_VariantExpansion_Logging_NoCrash()
        //{
        //    // Preferred explicit large XML path (user-specified)
        //    const string preferredLargeXml = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";
        //    const string preferredOutputDir = @"E:\Projects\Mann_11cell_lines\A549\A549_1"; // Force all output here

        //    // Ensure output directory exists
        //    try
        //    {
        //        if (!Directory.Exists(preferredOutputDir))
        //        {
        //            Directory.CreateDirectory(preferredOutputDir);
        //        }
        //    }
        //    catch
        //    {
        //        Assert.Inconclusive($"Cannot create/access output directory: {preferredOutputDir}");
        //        return;
        //    }
             
        //    string dbDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests");
        //    string overridePath = Environment.GetEnvironmentVariable("MZLIB_LARGE_XML") ?? "";
        //    string chosenPath = null;

        //    if (File.Exists(preferredLargeXml))
        //    {
        //        chosenPath = preferredLargeXml;
        //    }
        //    else if (!string.IsNullOrWhiteSpace(overridePath) && File.Exists(overridePath))
        //    {
        //        chosenPath = overridePath;
        //    }
        //    else if (Directory.Exists(dbDir))
        //    {
        //        chosenPath = Directory.GetFiles(dbDir, "*.xml")
        //            .OrderByDescending(f => new FileInfo(f).Length)
        //            .FirstOrDefault();
        //    }

        //    if (chosenPath == null)
        //    {
        //        Assert.Inconclusive("No XML database file found to run large variant logging diagnostic.");
        //        return;
        //    }

        //    string logPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "protein_variant_log.txt");
        //    var sb = new StringBuilder(1 << 16);
        //    sb.AppendLine("=== Protein Variant Expansion Diagnostic ===");
        //    sb.AppendLine($"Timestamp: {DateTime.Now:O}");
        //    sb.AppendLine($"InputFile: {chosenPath}");
        //    var fi = new FileInfo(chosenPath);
        //    sb.AppendLine($"FileSize: {fi.Length:N0} bytes  LastWrite: {fi.LastWriteTime}");
        //    sb.AppendLine("Parameters: maxVariantsPerIsoform=4 maxVariantIsoforms=400");
        //    sb.AppendLine();

        //    List<Protein> proteins = null;
        //    try
        //    {
        //        proteins = ProteinDbLoader.LoadProteinXML(
        //            chosenPath,
        //            generateTargets: true,
        //            decoyType: DecoyType.None,
        //            allKnownModifications: Enumerable.Empty<Modification>(),
        //            isContaminant: false,
        //            modTypesToExclude: null,
        //            unknownModifications: out var _,
        //            maxSequenceVariantsPerIsoform: 0,      // load base entries only first
        //            maxSequenceVariantIsoforms: 1);
        //    }
        //    catch (Exception ex)
        //    {
        //        sb.AppendLine("[FATAL] Exception during initial XML load:");
        //        sb.AppendLine(ex.ToString());
        //        File.WriteAllText(logPath, sb.ToString());
        //        Assert.Fail("Failed to load base XML. See log.");
        //        return;
        //    }

        //    if (proteins == null || proteins.Count == 0)
        //    {
        //        sb.AppendLine("[WARN] No proteins loaded; aborting variant expansion.");
        //        File.WriteAllText(logPath, sb.ToString());
        //        Assert.Inconclusive("No proteins loaded from selected XML.");
        //        return;
        //    }

        //    sb.AppendLine($"[INFO] Base proteins loaded: {proteins.Count}");
        //    sb.AppendLine();

        //    int proteinsAttempted = 0;
        //    int proteinsWithVariants = 0;
        //    int totalVariantIsoforms = 0;
        //    int totalExceptions = 0;

        //    foreach (var prot in proteins)
        //    {
        //        proteinsAttempted++;
        //        try
        //        {
        //            var varList = prot.GetVariantBioPolymers(
        //                maxSequenceVariantsPerIsoform: 4,
        //                minAlleleDepth: 1,
        //                maxSequenceVariantIsoforms: 400);

        //            // GetVariantBioPolymers returns list including base if combinatorics > 0; filter strict variants
        //            var distinct = varList
        //                .GroupBy(v => v.Accession)
        //                .Select(g => g.First())
        //                .ToList();

        //            int variantCount = distinct.Count - 1; // subtract base
        //            if (variantCount > 0)
        //            {
        //                proteinsWithVariants++;
        //                totalVariantIsoforms += variantCount;
        //            }

        //            sb.Append($"[OK] {prot.Accession} Len:{prot.Length} VariantsDefined:{prot.SequenceVariations?.Count ?? 0} Generated:{variantCount}");

        //            // Quick audit of each generated variant (length & attribute agreement, error markers)
        //            if (variantCount > 0)
        //            {
        //                var audits = new List<string>();
        //                foreach (var iso in distinct.Where(v => !ReferenceEquals(v, prot)))
        //                {
        //                    bool lenAttrMismatch = iso.UniProtSequenceAttributes != null &&
        //                                           iso.UniProtSequenceAttributes.Length != iso.Length;
        //                    string token = string.Join("+",
        //                        iso.AppliedSequenceVariations.Select(v => v.SimpleString()));
        //                    if (string.IsNullOrEmpty(token))
        //                        token = "NO_TOKEN";

        //                    audits.Add(token +
        //                               (lenAttrMismatch ? "(LenAttrMismatch)" : "") +
        //                               (iso.BaseSequence.Length == prot.BaseSequence.Length ? "" : "(SeqLenΔ)"));
        //                }
        //                if (audits.Count > 0)
        //                    sb.Append(" [" + string.Join(", ", audits.Take(15)) + (audits.Count > 15 ? ", ..." : "") + "]");
        //            }

        //            sb.AppendLine();
        //        }
        //        catch (Exception ex)
        //        {
        //            totalExceptions++;
        //            sb.AppendLine($"[ERR] {prot.Accession} Exception: {ex.GetType().Name} - {ex.Message}");
        //        }

        //        // Periodically flush to disk for very large sets
        //        if (proteinsAttempted % 250 == 0)
        //        {
        //            File.WriteAllText(logPath, sb.ToString());
        //        }
        //    }

        //    sb.AppendLine();
        //    sb.AppendLine("=== Summary ===");
        //    sb.AppendLine($"ProteinsAttempted: {proteinsAttempted}");
        //    sb.AppendLine($"ProteinsWithVariants: {proteinsWithVariants}");
        //    sb.AppendLine($"TotalVariantIsoforms (excl. bases): {totalVariantIsoforms}");
        //    sb.AppendLine($"Exceptions: {totalExceptions}");
        //    sb.AppendLine("================");

        //    File.WriteAllText(logPath, sb.ToString());

        //    // Soft assertions: test passes as long as no catastrophic failure
        //    Assert.That(File.Exists(logPath), "Log file not created.");
        //    Assert.That(proteinsAttempted, Is.GreaterThan(0), "No proteins processed.");
        //    // Do not fail on variant exceptions; log is the artifact for inspection.
        //}
    }
}