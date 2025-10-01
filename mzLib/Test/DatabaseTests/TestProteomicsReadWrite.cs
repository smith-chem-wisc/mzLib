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
using NUnit.Framework.Legacy;

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
    }
}