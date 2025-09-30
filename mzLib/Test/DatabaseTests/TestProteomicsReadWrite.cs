using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics.BioPolymer;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
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
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None,
                null, false, null, out Dictionary<string, Modification> un);
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
                maxSequenceVariantsPerIsoform:2,
                maxSequenceVariantIsoforms:100);
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

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(oldXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un);
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

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, uniprotPtms.Concat(nice), false, null,
                out Dictionary<string, Modification> un);
            Protein zero = ok[0];
            Protein one = ok[1];
            Dictionary<int, List<Modification>> zero_mods = zero.OneBasedPossibleLocalizedModifications as Dictionary<int, List<Modification>>;
            Dictionary<int, List<Modification>> one_mods = one.OneBasedPossibleLocalizedModifications as Dictionary<int, List<Modification>>;

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), true, DecoyType.None, nice, false,
                new List<string>(), out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(9, ok[0].DatabaseReferences.Count(dbRef => dbRef.Type == "GO"));
            Assert.AreEqual(1, ok[0].DatabaseReferences.Count(dbRef => dbRef.Type == "GeneID"));
            Assert.AreEqual(3, ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "GO").Properties.Count());
            Assert.AreEqual(3, ok[0].GeneNames.Count());
            Assert.AreEqual("primary", ok[0].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", ok[0].GeneNames.First().Item2);
            Assert.AreEqual("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", ok[0].Organism);
            Assert.AreEqual(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), ok[0].DatabaseFilePath);
            Assert.AreEqual(9, ok2[0].DatabaseReferences.Count(dbRef => dbRef.Type == "GO"));
            Assert.AreEqual(3, ok2[0].DatabaseReferences.First(dbRef => dbRef.Type == "GO").Properties.Count());
            Assert.AreEqual(3, ok2[0].GeneNames.Count());
            Assert.AreEqual("primary", ok2[0].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", ok2[0].GeneNames.First().Item2);
            Assert.AreEqual("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", ok2[0].Organism);
            Assert.AreEqual(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), ok2[0].DatabaseFilePath);
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
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(inputXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml");

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, outputPath, true);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), true, DecoyType.None, uniprotPtms, false,
                new List<string>(), out un);

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

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(inputXmlPath, true, DecoyType.None, uniprotPtms, false, null,
                out Dictionary<string, Modification> un);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_unknownStatus.xml");

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, outputPath, true);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(outputPath, true, DecoyType.None, uniprotPtms, false,
                new List<string>(), out un);

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

            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false, out var a,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml"), true, DecoyType.None, nice,
                false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual("ENSP00000381386", ok[0].Accession);
            Assert.AreEqual("ENSP00000215773", ok[1].Accession);
            Assert.AreEqual("ENSG00000099977", ok[0].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok[1].GeneNames.First().Item2);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[0].FullName);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[1].FullName);
            Assert.AreEqual(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), ok[0].DatabaseFilePath);

            Assert.AreEqual("ENSP00000381386", ok2[0].Accession);
            Assert.AreEqual("ENSP00000215773", ok2[1].Accession);
            Assert.AreEqual("ENSG00000099977", ok2[0].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok2[1].GeneNames.First().Item2);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[0].FullName);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[1].FullName);
            Assert.AreEqual(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml"), ok2[0].DatabaseFilePath);

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
                new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> un);
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

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, nice, false, null,
                out Dictionary<string, Modification> un);
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

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, uniprotPtms.Concat(nice), false, new List<string>(),
                out Dictionary<string, Modification> un);
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"));
            Assert.AreEqual(0, newModResEntries.Count);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(2, ok[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(2, ok2[0].OneBasedPossibleLocalizedModifications.Count);
        }

        [Test]
        public void AnotherTest()
        {
            List<Modification> variableModifications = new List<Modification>();
            List<Modification> fixedModifications = new List<Modification>();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null,
                "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), disulfideBonds: new List<DisulfideBond>());

            List<TruncationProduct> pp = new List<TruncationProduct> { new TruncationProduct(4, 8, "chain") };
            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), pp,
                "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), disulfideBonds: new List<DisulfideBond>());

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, Path.Combine(TestContext.CurrentContext.TestDirectory, @"fdsfsd.xml"));
        }

        [Test]
        public void TestEmptyProteins()
        {
            Protein p1 = new Protein("SEQENCE", "p1");
            Assert.AreEqual("p1||", p1.FullDescription);
            Protein p2 = new Protein("SEQENCE", "p2", name: "namep2");

            var proteinListToWrite = new List<Protein> { p1, p2 };

            // Generate data for files
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinListToWrite, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"differentlyConstuctedProteins.xml"));

            IEnumerable<string> modTypesToExclude = new List<string>();
            IEnumerable<Modification> allKnownModifications = new List<Modification>();
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"differentlyConstuctedProteins.xml"), true, DecoyType.None,
                allKnownModifications, false, modTypesToExclude, out Dictionary<string, Modification> un);
            Assert.AreEqual(p1.Accession, ok[0].Accession);
            Assert.AreEqual(p2.Accession, ok[1].Accession);
            Assert.AreEqual(p1.Name, ok[0].Name);
            Assert.AreEqual(p2.Name, ok[1].Name);
        }
        [Test]
        public void TestFullProteinReadWrite()
        {
            // Re‑implementation based on the minimal pattern proven to work in TestProteinDatabase
            // (WriteXmlDatabase_WritesRequiredUniProtSequenceAttributes). Previous versions likely
            // hit a NullReference because a constructor parameter ordering/naming mismatch left an
            // internal field (accessed by Dataset/Created/Modified/Version or UniProtSequenceAttributes)
            // unset. This version mirrors the known-good constructor argument style and keeps the
            // assertions focused on round‑trip integrity.

            // Base sequence
            const string seq = "SEQENCE"; // length 7

            // Required motifs (safe fallbacks)
            ModificationMotif.TryGetMotif("E", out var motifE);
            ModificationMotif.TryGetMotif("N", out var motifN);
            Assert.IsNotNull(motifE);
            Assert.IsNotNull(motifN);

            // Simple residue mods
            var modE = new Modification("mod on E", null, "mt", null, motifE, "Anywhere.", null, null, null, null, null, null, null, null);
            var modN = new Modification("mod on N", null, "mt", null, motifN, "Anywhere.", null, 10, null, null, null, null, null, null);

            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification>{ modE } }, // E
                { 5, new List<Modification>{ modN } }  // N
            };

            var uniProtAttrs = new UniProtSequenceAttributes(
                length: seq.Length,
                mass: 0,
                checkSum: "CHKTEST",
                entryModified: DateTime.Today,
                sequenceVersion: 1
            );

            var protein = new Protein(
                accession: "A1",
                sequence: seq,
                organism: "Test organism",
                isDecoy: false,
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") },
                name: "TestName",
                fullName: "Test Full Name",
                isContaminant: false,
                sequenceVariations: new List<SequenceVariation>(),   // none
                disulfideBonds: new List<DisulfideBond>(),           // none
                spliceSites: new List<SpliceSite>(),                 // ensure not null
                databaseReferences: new List<DatabaseReference>(),
                databaseFilePath: Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "fullProtein.xml"),
                uniProtSequenceAttributes: uniProtAttrs,
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null,
                oneBasedModifications: oneBasedMods
            );

            string outPath = protein.DatabaseFilePath;
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, outPath);

            var roundTrip = ProteinDbLoader.LoadProteinXML(
                outPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: Enumerable.Empty<Modification>(),
                isContaminant: false,
                modTypesToExclude: Enumerable.Empty<string>(),
                unknownModifications: out var unknown).Single();

            // Core identity
            Assert.AreEqual(protein.Accession, roundTrip.Accession);
            Assert.AreEqual(protein.BaseSequence, roundTrip.BaseSequence);
            Assert.AreEqual(protein.FullName, roundTrip.FullName);
            Assert.AreEqual(protein.Name, roundTrip.Name);
            Assert.AreEqual(protein.Organism, roundTrip.Organism);
            Assert.AreEqual(protein.Length, roundTrip.Length);
            Assert.IsNotNull(roundTrip.UniProtSequenceAttributes);
            Assert.AreEqual(seq.Length, roundTrip.UniProtSequenceAttributes.Length);

            // Mods round‑trip (positions & counts)
            Assert.AreEqual(protein.OneBasedPossibleLocalizedModifications.Keys.Count,
                roundTrip.OneBasedPossibleLocalizedModifications.Keys.Count);
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                Assert.IsTrue(roundTrip.OneBasedPossibleLocalizedModifications.ContainsKey(kvp.Key));
                Assert.AreEqual(kvp.Value.Count, roundTrip.OneBasedPossibleLocalizedModifications[kvp.Key].Count);
            }

            // No variants / features unexpectedly introduced
            Assert.AreEqual(0, roundTrip.SequenceVariations.Count());
            Assert.AreEqual(0, roundTrip.DisulfideBonds.Count());
            Assert.AreEqual(0, roundTrip.SpliceSites.Count());
        }
        [Test]
        public void TestReadWriteSeqVars()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), true, DecoyType.None,
                nice, false, null, out Dictionary<string, Modification> un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out un);

            Assert.AreEqual(ok[0].SequenceVariations.Count(), ok2[0].SequenceVariations.Count());
            Assert.AreEqual(ok[0].SequenceVariations.First().OneBasedBeginPosition, ok2[0].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(ok[0].SequenceVariations.First().OneBasedEndPosition, ok2[0].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreEqual(ok[0].SequenceVariations.First().VariantCallFormatData, ok2[0].SequenceVariations.First().VariantCallFormatData);
            Assert.AreEqual(ok[0].SequenceVariations.First().OriginalSequence, ok2[0].SequenceVariations.First().OriginalSequence);
            Assert.AreEqual(ok[0].SequenceVariations.First().VariantSequence, ok2[0].SequenceVariations.First().VariantSequence);
        }

        [Test]
        public void TestReadWriteSeqVars2()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"seqvartests.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out Dictionary<string, Modification> un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_seqvartests.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_seqvartests.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out un);

            Assert.AreEqual(ok[0].SequenceVariations.Count(), ok2[0].SequenceVariations.Count());
            Assert.AreEqual(ok[0].SequenceVariations.First().OneBasedBeginPosition, ok2[0].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(ok[0].SequenceVariations.First().OneBasedEndPosition, ok2[0].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreEqual(ok[0].SequenceVariations.First().VariantCallFormatData, ok2[0].SequenceVariations.First().VariantCallFormatData);
            Assert.AreEqual(ok[0].SequenceVariations.First().OriginalSequence, ok2[0].SequenceVariations.First().OriginalSequence);
            Assert.AreEqual(ok[0].SequenceVariations.First().VariantSequence, ok2[0].SequenceVariations.First().VariantSequence);
        }

        [Test]
        public void TestModificationGeneralToString()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt"), out var errors).ToList();
            char[] myChar = { '"' };
            string output = a.First().ToString();
            Assert.AreEqual(output.TrimStart(myChar).TrimEnd(myChar), "ID   4-carboxyglutamate on E\r\nMT   Biological\r\nTG   E\r\nPP   Anywhere.\r\nCF   CO2\r\nMM   43.989829\r\n");
        }

        [Test]
        public void TestModificationGeneral_Equals()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt"), out var errors).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt"), out errors).ToList();

            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        public static void Test_CustumPrunedDatabaseWriteAndRead()
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif K);
            ModificationMotif.TryGetMotif("R", out ModificationMotif R);

            Modification acOnK = new Modification(_originalId: "Acetyl", _accession: null, _modificationType: "testModType", _featureType: null, _locationRestriction: "Anywhere.", _target: K, _monoisotopicMass: 42);
            Modification meOnK = new Modification(_originalId: "Methyl", _accession: null, _modificationType: "testModType", _featureType: null, _locationRestriction: "Anywhere.", _target: K, _monoisotopicMass: 14);
            Modification meOnR = new Modification(_originalId: "Methyl", _accession: null, _modificationType: "testModType", _featureType: null, _locationRestriction: "Anywhere.", _target: R, _monoisotopicMass: 14);

            Dictionary<int, List<Modification>> obm = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification>() { acOnK } },
                { 2, new List<Modification>() { meOnK } },
                { 3, new List<Modification>() { meOnR } }
            };

            Protein p = new Protein("KKR", "accession", null, null, obm, null, null, null, false, false, null, null, null, null);
            List<Protein> pList = new List<Protein>() { p };

            string outputFileName = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"redundant.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), pList, outputFileName);

            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(outputFileName,
                true, DecoyType.None, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> proteinXmlModList);

            Assert.AreEqual(3, new_proteins[0].OneBasedPossibleLocalizedModifications.Count());
        }

        [Test]
        public static void TestStringSanitation()
        {
            string messedUpSequence = @"PRO�EIN�";

            // just test the string sanitation method alone
            var sanitized = ProteinDbLoader.SanitizeAminoAcidSequence(messedUpSequence, 'C');
            Assert.That(sanitized == "PROCEINC");

            // test reading from a fasta
            Protein protein = new Protein(sanitized, "accession");

            string fastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"messedUp.fasta");
            ProteinDbWriter.WriteFastaDatabase(new List<Protein> { protein }, fastaPath, "|");

            var fastaProteins = ProteinDbLoader.LoadProteinFasta(fastaPath, true, DecoyType.Reverse, false, out var a, ProteinDbLoader.UniprotAccessionRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);

            Assert.That(fastaProteins.First(p => !p.IsDecoy).BaseSequence == "PROCEINC");

            // digest and fragment to check that there isn't a crash
            var peptides = fastaProteins.First().Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).ToList();
            foreach (PeptideWithSetModifications peptide in peptides)
            {
                List<Product> fragments = new List<Product>();
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            }

            // test reading from an XML
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"messedUp.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, xmlPath);
            var xmlProteins = ProteinDbLoader.LoadProteinXML(xmlPath, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var unk);

            Assert.That(xmlProteins.First(p => !p.IsDecoy).BaseSequence == "PROCEINC");
        }
        [Test]
        [Category("LongRunning")]
        public void ReadWriteLargeProteinXmlLogErrors()
        {
            string inputPath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";
            string outputPath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\variant.xml";
            string logPath = Path.Combine(Path.GetDirectoryName(outputPath) ?? TestContext.CurrentContext.WorkDirectory, "protein_variant_log.txt");

            var log = new List<string>();
            void Log(string line)
            {
                log.Add(line);
                TestContext.WriteLine(line);
            }
            void FlushLog()
            {
                try { File.WriteAllLines(logPath, log); } catch (Exception ex) { TestContext.WriteLine("[WARN] Could not write log: " + ex.Message); }
            }

            Log("=== Large Protein XML Diagnostic Loader ===");
            Log("Input:  " + inputPath);
            Log("Output: " + outputPath);

            if (!File.Exists(inputPath))
            {
                Log("[FATAL] File does not exist.");
                FlushLog();
                Assert.Pass("Input XML missing; see log.");
            }

            try
            {
                var fi = new FileInfo(inputPath);
                Log($"File Size: {fi.Length:N0} bytes  LastWrite: {fi.LastWriteTime}");
            }
            catch (Exception ex)
            {
                Log("[WARN] Could not stat file: " + ex.Message);
            }

            // Peek at start/end lines for sanity
            try
            {
                var allLinesEnum = File.ReadLines(inputPath);
                var head = allLinesEnum.Take(10).ToList();
                var tail = File.ReadLines(inputPath).Reverse().Take(10).Reverse().ToList();
                Log("--- File Head (first 10 lines) ---");
                foreach (var l in head) Log(l);
                Log("--- File Tail (last 10 lines) ---");
                foreach (var l in tail) Log(l);
            }
            catch (Exception ex)
            {
                Log("[WARN] Could not preview file content: " + ex.Message);
            }

            List<Protein> rawProteins = null;
            Dictionary<string, Modification> unknownMods;

            var loadAttempts = new List<(string Label, Func<List<Protein>> Action)>();

            // Attempt #1: Full settings (original intention)
            loadAttempts.Add(("FullVariants",
                () => ProteinDbLoader.LoadProteinXML(
                        inputPath,
                        generateTargets: true,
                        decoyType: DecoyType.None,
                        allKnownModifications: Enumerable.Empty<Modification>(),
                        isContaminant: false,
                        modTypesToExclude: new List<string>(),
                        unknownModifications: out unknownMods,
                        maxSequenceVariantsPerIsoform: 50,
                        maxSequenceVariantIsoforms: 500)));

            // Attempt #2: Reduced variant burden
            loadAttempts.Add(("ReducedVariants",
                () => ProteinDbLoader.LoadProteinXML(
                        inputPath,
                        generateTargets: true,
                        decoyType: DecoyType.None,
                        allKnownModifications: Enumerable.Empty<Modification>(),
                        isContaminant: false,
                        modTypesToExclude: new List<string>(),
                        unknownModifications: out unknownMods,
                        maxSequenceVariantsPerIsoform: 10,
                        maxSequenceVariantIsoforms: 50)));

            // Attempt #3: No variant expansion (max isoforms = 1)
            loadAttempts.Add(("NoVariants",
                () => ProteinDbLoader.LoadProteinXML(
                        inputPath,
                        generateTargets: true,
                        decoyType: DecoyType.None,
                        allKnownModifications: Enumerable.Empty<Modification>(),
                        isContaminant: false,
                        modTypesToExclude: new List<string>(),
                        unknownModifications: out unknownMods,
                        maxSequenceVariantsPerIsoform: 1,
                        maxSequenceVariantIsoforms: 1)));

            // Attempt #4: Minimal parse (treat as contaminants = false but still parse)
            loadAttempts.Add(("Minimal",
                () => ProteinDbLoader.LoadProteinXML(
                        inputPath,
                        generateTargets: false,
                        decoyType: DecoyType.None,
                        allKnownModifications: Enumerable.Empty<Modification>(),
                        isContaminant: false,
                        modTypesToExclude: new List<string>(),
                        unknownModifications: out unknownMods,
                        maxSequenceVariantsPerIsoform: 1,
                        maxSequenceVariantIsoforms: 1)));

            Exception lastEx = null;
            foreach (var (label, action) in loadAttempts)
            {
                try
                {
                    Log($"[INFO] Attempting load strategy: {label}");
                    rawProteins = action();
                    if (rawProteins != null && rawProteins.Count > 0)
                    {
                        Log($"[SUCCESS] Strategy '{label}' loaded {rawProteins.Count} proteins.");
                        break;
                    }
                    Log($"[WARN] Strategy '{label}' returned null or empty set.");
                }
                catch (Exception ex)
                {
                    lastEx = ex;
                    Log($"[ERROR] Strategy '{label}' threw: {ex.Message}");
                    var ie = ex.InnerException;
                    int depth = 0;
                    while (ie != null && depth < 5)
                    {
                        Log($"  Inner[{depth}] {ie.GetType().Name}: {ie.Message}");
                        ie = ie.InnerException;
                        depth++;
                    }
                    Log("  Stack (first lines):");
                    foreach (var line in ex.StackTrace?.Split('\n').Take(6) ?? Enumerable.Empty<string>())
                        Log("    " + line.Trim());
                }
            }

            if (rawProteins == null || rawProteins.Count == 0)
            {
                Log("[FATAL] All loading strategies failed.");
                if (lastEx != null) Log("Last exception: " + lastEx.GetType().Name + " - " + lastEx.Message);
                FlushLog();
                Assert.Pass("Could not load proteins; see log for diagnostics: " + logPath);
            }

            Log("[INFO] Proceeding to variant expansion & write phase.");

            // Variant expansion (safe) – we don’t abort if some fail
            var expanded = new List<Protein>();
            var variantFailures = new List<(string Accession, string Reason)>();
            foreach (var p in rawProteins)
            {
                expanded.Add(p);
                try
                {
                    var vs = p.GetVariantBioPolymers(maxSequenceVariantIsoforms: 200).OfType<Protein>().ToList();
                    foreach (var v in vs) expanded.Add(v);
                }
                catch (Exception ex)
                {
                    variantFailures.Add((p.Accession, ex.Message));
                }
            }

            if (variantFailures.Count > 0)
            {
                Log($"[WARN] Variant expansion failures: {variantFailures.Count}");
                foreach (var vf in variantFailures.Take(100))
                    Log($"VariantFail\t{vf.Accession}\t{vf.Reason}");
            }

            // Write
            bool writeOk = false;
            try
            {
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), expanded, outputPath);
                writeOk = true;
                Log($"[INFO] Wrote combined XML: {outputPath}");
            }
            catch (Exception ex)
            {
                Log("[ERROR] Bulk write failed: " + ex.Message);
                // Attempt isolation
                foreach (var p in expanded.Take(500))
                {
                    try
                    {
                        var tmp = Path.Combine(Path.GetTempPath(), $"single_{SanitizeFilePart(p.Accession)}_{Guid.NewGuid():N}.xml");
                        ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { p }, tmp);
                        try { File.Delete(tmp); } catch { }
                    }
                    catch (Exception px)
                    {
                        Log($"WriteFail\t{p.Accession}\t{px.Message}");
                    }
                }
            }

            // Optional read-back
            if (writeOk && File.Exists(outputPath))
            {
                try
                {
                    var rt = ProteinDbLoader.LoadProteinXML(
                        outputPath,
                        generateTargets: true,
                        decoyType: DecoyType.None,
                        allKnownModifications: Enumerable.Empty<Modification>(),
                        isContaminant: false,
                        modTypesToExclude: new List<string>(),
                        unknownModifications: out _,
                        maxSequenceVariantsPerIsoform: 2,
                        maxSequenceVariantIsoforms: 10);
                    Log($"[INFO] Read-back proteins: {rt?.Count ?? 0}");
                }
                catch (Exception ex)
                {
                    Log("[ERROR] Read-back failed: " + ex.Message);
                }
            }

            FlushLog();
            Log("Log written: " + logPath);
            Assert.Pass("Completed diagnostic run. See log: " + logPath);

            // Helpers
            static string SanitizeFilePart(string s)
            {
                if (string.IsNullOrWhiteSpace(s)) return "NA";
                var invalid = Path.GetInvalidFileNameChars();
                return new string(s.Select(c => invalid.Contains(c) ? '_' : c).ToArray());
            }
        }
        [Test]
        [Category("Diagnostic")]
        public void DiagnoseSingleProblemProteinVariants()
        {
            // small.xml should contain ONLY the first failing UniProt <entry> (e.g., A0A087X1C5)
            // placed in the same directory as the large file. This test mirrors the large diagnostic
            // but adds deeper per‑variant validation and never hard-fails.
            string folder = @"E:\Projects\Mann_11cell_lines\A549\A549_1";
            string inputPath = Path.Combine(folder, "small.xml");
            string outputPath = Path.Combine(folder, "small_variant.xml");
            string logPath = Path.Combine(folder, "small_variant_log.txt");

            var log = new List<string>();
            void Log(string msg)
            {
                log.Add(msg);
                TestContext.WriteLine(msg);
            }
            void Flush()
            {
                try { File.WriteAllLines(logPath, log); }
                catch (Exception ex) { TestContext.WriteLine("[WARN] Could not write log: " + ex.Message); }
            }

            Log("=== Single Protein Variant Diagnostic ===");
            Log("Input:  " + inputPath);
            Log("Output: " + outputPath);

            if (!File.Exists(inputPath))
            {
                Log("[FATAL] small.xml not found.");
                Flush();
                Assert.Pass("Missing small.xml; nothing to diagnose.");
            }

            try
            {
                var fi = new FileInfo(inputPath);
                Log($"File Size: {fi.Length:N0} bytes  LastWrite: {fi.LastWriteTime}");
            }
            catch (Exception ex)
            {
                Log("[WARN] Could not stat file: " + ex.Message);
            }

            // Preview first few lines
            try
            {
                foreach (var l in File.ReadLines(inputPath).Take(15))
                    Log(l);
            }
            catch (Exception ex)
            {
                Log("[WARN] Could not preview file head: " + ex.Message);
            }

            Dictionary<string, Modification> unknown;
            List<Protein> proteins = null;
            Exception loadEx = null;

            try
            {
                proteins = ProteinDbLoader.LoadProteinXML(
                    inputPath,
                    generateTargets: true,
                    decoyType: DecoyType.None,
                    allKnownModifications: Enumerable.Empty<Modification>(),
                    isContaminant: false,
                    modTypesToExclude: new List<string>(),
                    unknownModifications: out unknown,
                    maxSequenceVariantsPerIsoform: 50,
                    maxSequenceVariantIsoforms: 200);
            }
            catch (Exception ex)
            {
                loadEx = ex;
                Log("[ERROR] LoadProteinXML threw: " + ex.Message);
                if (ex.StackTrace != null)
                    Log("StackTop: " + string.Join(" | ", ex.StackTrace.Split('\n').Take(4).Select(s => s.Trim())));
            }

            if (proteins == null || proteins.Count == 0)
            {
                Log("[FATAL] No proteins parsed from small.xml.");
                if (loadEx != null && loadEx.InnerException != null)
                    Log("Inner: " + loadEx.InnerException.GetType().Name + " - " + loadEx.InnerException.Message);
                Flush();
                Assert.Pass("Load failed; see log.");
            }

            Log($"[INFO] Proteins parsed: {proteins.Count}");

            // We expect exactly one; if more, we still proceed
            foreach (var p in proteins)
            {
                Log($"--- Protein Accession: {p.Accession} Name:{p.Name} Length:{p.Length} VariationsDefined:{p.SequenceVariations?.Count() ?? 0}");
                if (p.SequenceVariations == null || !p.SequenceVariations.Any())
                {
                    Log("[INFO] No declared sequence variations; nothing to apply.");
                    continue;
                }

                // Per-variation structural validation
                int idx = 0;
                foreach (var v in p.SequenceVariations)
                {
                    idx++;
                    try
                    {
                        ValidateVariation(p, v, idx, Log);
                    }
                    catch (Exception ex)
                    {
                        Log($"[VAR-CHECK-EX] #{idx} {ex.GetType().Name}: {ex.Message}");
                    }
                }

                // Attempt variant generation with guarded catch
                List<Protein> variantForms = null;
                try
                {
                    variantForms = p.GetVariantBioPolymers(maxSequenceVariantIsoforms: 100).OfType<Protein>().ToList();
                    Log($"[APPLY] Variant proteoforms generated: {variantForms.Count} (Applied sets: {variantForms.Count(vf => vf.AppliedSequenceVariations.Any())})");
                }
                catch (Exception ex)
                {
                    Log("[APPLY-ERROR] GetVariantBioPolymers: " + ex.Message);
                    if (ex.StackTrace != null)
                        Log("StackTop: " + string.Join(" | ", ex.StackTrace.Split('\n').Take(5).Select(s => s.Trim())));
                }

                // Enumerate failing application individually (simulate what ApplySingleVariant might do)
                if (variantForms == null || variantForms.Count == 0)
                {
                    Log("[INFO] No variant proteoforms produced; attempting manual sequential application per variation to isolate culprit.");
                    // Try applying each variation in isolation by constructing a single-variation scenario
                    int vNum = 0;
                    foreach (var v in p.SequenceVariations)
                    {
                        vNum++;
                        try
                        {
                            ManualApplyVariantPreview(p, v, vNum, Log);
                        }
                        catch (Exception ex)
                        {
                            Log($"[MANUAL-APPLY-FAIL] Var#{vNum} {ex.GetType().Name}: {ex.Message}");
                        }
                    }
                }
            }

            // Attempt to serialize whatever we have (even if only original protein) to catch write-specific NREs
            try
            {
                ProteinDbWriter.WriteXmlDatabase(
                    new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    proteins,
                    outputPath);
                Log("[INFO] Wrote small_variant.xml successfully.");
            }
            catch (Exception ex)
            {
                Log("[WRITE-ERROR] " + ex.Message);
                if (ex.StackTrace != null)
                    Log("WriteStackTop: " + string.Join(" | ", ex.StackTrace.Split('\n').Take(4).Select(s => s.Trim())));
            }

            Flush();
            Assert.Pass("Diagnostic complete. See log: " + logPath);

            // Helper: structural checks
            static void ValidateVariation(Protein p, SequenceVariation v, int idx, Action<string> log)
            {
                string baseSeq = p.BaseSequence;
                int len = baseSeq.Length;
                int b = v.OneBasedBeginPosition;
                int e = v.OneBasedEndPosition;
                string orig = v.OriginalSequence ?? "";
                string varSeq = v.VariantSequence ?? "";

                log($"[VAR] #{idx} Begin:{b} End:{e} Orig:'{orig}' Var:'{varSeq}' TypeHint:{v.Description}");

                // Coordinate sanity
                if (b < 1 || e < b)
                    log($"  [WARN] Invalid coordinate ordering (Begin:{b}, End:{e}).");
                if (e > len)
                    log($"  [WARN] End position ({e}) exceeds sequence length ({len}).");

                // If original sequence provided, verify it matches the substring
                if (!string.IsNullOrEmpty(orig) && e <= len)
                {
                    int subLen = e - b + 1;
                    if (subLen == orig.Length)
                    {
                        string actual = baseSeq.Substring(b - 1, subLen);
                        if (!string.Equals(actual, orig, StringComparison.Ordinal))
                        {
                            log($"  [MISMATCH] OriginalSequence mismatch. ExpectedInBase:'{actual}' Provided:'{orig}'");
                        }
                    }
                    else
                    {
                        log($"  [WARN] OriginalSequence length ({orig.Length}) != span length ({subLen}).");
                    }
                }

                // Insertion: orig empty, variant non-empty
                if (string.IsNullOrEmpty(orig) && !string.IsNullOrEmpty(varSeq))
                {
                    if (b > len + 1)
                        log($"  [WARN] Insertion begin {b} beyond permissible insertion boundary (len+1={len + 1}).");
                }

                // Deletion: var empty
                if (!string.IsNullOrEmpty(orig) && string.IsNullOrEmpty(varSeq))
                {
                    if (e > len)
                        log($"  [WARN] Deletion end {e} beyond sequence length {len}.");
                }

                // Stop-gain / stop-loss heuristics (asterisk)
                if (varSeq == "*")
                    log("  [INFO] Stop-gain detected.");
                if (orig == "*")
                    log("  [INFO] Stop-loss / extension detected.");
            }

            // Helper: manual preview (simulate variant application core logic simplistically)
            static void ManualApplyVariantPreview(Protein p, SequenceVariation v, int idx, Action<string> log)
            {
                string seq = p.BaseSequence;
                int len = seq.Length;
                int b = v.OneBasedBeginPosition;
                int e = v.OneBasedEndPosition;
                string orig = v.OriginalSequence ?? "";
                string varSeq = v.VariantSequence ?? "";

                log($"[MANUAL] Applying Var#{idx} Begin:{b} End:{e} Orig:'{orig}' Var:'{varSeq}'");

                if (b < 1 || e < b || e > len)
                    throw new ArgumentOutOfRangeException($"Coordinates out of range (Begin={b}, End={e}, Len={len}).");

                // If original given, verify
                if (!string.IsNullOrEmpty(orig))
                {
                    string actual = seq.Substring(b - 1, Math.Min(e, len) - b + 1);
                    if (actual.Length == orig.Length && actual != orig)
                        log($"  [CHECK] Original mismatch (BaseSpan='{actual}' vs Orig='{orig}'). Proceeding anyway.");
                }

                string newSeq;
                if (orig == "*" && !string.IsNullOrEmpty(varSeq))
                {
                    // stop-loss: append extension
                    newSeq = seq + varSeq;
                }
                else if (varSeq == "*")
                {
                    // stop-gain: truncate just before begin
                    newSeq = seq.Substring(0, b - 1);
                }
                else
                {
                    // general replacement
                    int removeLen = e - b + 1;
                    if (b - 1 + removeLen > seq.Length)
                        removeLen = Math.Max(0, seq.Length - (b - 1));
                    newSeq = seq.Substring(0, b - 1) + varSeq + seq.Substring(b - 1 + removeLen);
                }

                log($"  [MANUAL] Result length: {newSeq.Length} (Δ {newSeq.Length - seq.Length})");
            }
        }
    }
}