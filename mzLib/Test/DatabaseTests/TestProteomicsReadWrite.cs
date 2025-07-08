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
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(new_mods, ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"));
            Assert.AreEqual(1, newModResEntries.Count);
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(2, ok[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(3, ok2[0].OneBasedPossibleLocalizedModifications.Count);
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
            ProteinDbWriter.WriteXmlDatabase(null, proteinListToWrite, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"differentlyConstuctedProteins.xml"));

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
            Modification mod = new Modification("mod1", null, "modType1", null, null, null, null, null, null, null, null, null, null, null);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification mod2 = new Modification("mod2 on E", null, "modType1", null, motif, "Anywhere.", null, null, null, null, null, null, null, null);
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif3);
            Modification mod3 = new Modification("mod3 on N", null, "modType1", null, motif3, "Anywhere.", null, 10, null, null, null, null, null, null);

            List<Tuple<string, string>> gene_names = new List<Tuple<string, string>> { new Tuple<string, string>("a", "b") };
            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {3, new List<Modification>{mod} },
                {4, new List<Modification>{mod2} },
                {5, new List<Modification>{mod3} }
            };
            List<TruncationProduct> proteolysisProducts = new List<TruncationProduct> { new TruncationProduct(1, 2, "propeptide") };

            string name = "testName";

            string full_name = "testFullName";

            List<DatabaseReference> databaseReferences = new List<DatabaseReference> {
                new DatabaseReference("type1", "id1", new List<Tuple<string, string>> { new Tuple<string, string>("e1", "e2") }) };

            List<SequenceVariation> sequenceVariations = new List<SequenceVariation> { new SequenceVariation(3,"Q", "N", "replace Q by N"),
            new SequenceVariation(3,4,"QE", "NN", "replace QE by NN")};

            List<DisulfideBond> disulfideBonds = new List<DisulfideBond> { new DisulfideBond(1, "ds1"), new DisulfideBond(2, 3, "ds2") };

            Protein originalProtein = new Protein(
                "SEQENCE",
                "a1",
                geneNames: gene_names,
                oneBasedModifications: oneBasedModifications,
                proteolysisProducts: proteolysisProducts,
                name: name,
                fullName: full_name,
                isDecoy: false,
                isContaminant: true,
                databaseReferences: databaseReferences,
                sequenceVariations: sequenceVariations,
                disulfideBonds: disulfideBonds,
                databaseFilePath: Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"));

            // Generate data for files
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { originalProtein },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"));

            IEnumerable<string> modTypesToExclude = new List<string>();
            IEnumerable<Modification> allKnownModifications = new List<Modification>();
            List<Protein> proteinReadFromXml = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"), true, DecoyType.None,
                allKnownModifications, true, modTypesToExclude, out Dictionary<string, Modification> unknownModifications);
            Assert.AreEqual(originalProtein.Accession, proteinReadFromXml[0].Accession);
            Assert.AreEqual(originalProtein.BaseSequence, proteinReadFromXml[0].BaseSequence);
            Assert.AreEqual(originalProtein.DatabaseReferences.First().Id, proteinReadFromXml[0].DatabaseReferences.First().Id);
            Assert.AreEqual(originalProtein.DatabaseReferences.First().Properties.First().Item1, proteinReadFromXml[0].DatabaseReferences.First().Properties.First().Item1);
            Assert.AreEqual(originalProtein.DatabaseReferences.First().Properties.First().Item2, proteinReadFromXml[0].DatabaseReferences.First().Properties.First().Item2);
            Assert.AreEqual(originalProtein.DatabaseReferences.First().Type, proteinReadFromXml[0].DatabaseReferences.First().Type);

            Assert.AreEqual(originalProtein.DisulfideBonds.First().Description, proteinReadFromXml[0].DisulfideBonds.First().Description);
            Assert.AreEqual(originalProtein.DisulfideBonds.First().OneBasedBeginPosition, proteinReadFromXml[0].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(originalProtein.DisulfideBonds.First().OneBasedEndPosition, proteinReadFromXml[0].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(originalProtein.DisulfideBonds.Last().Description, proteinReadFromXml[0].DisulfideBonds.Last().Description);
            Assert.AreEqual(originalProtein.DisulfideBonds.Last().OneBasedBeginPosition, proteinReadFromXml[0].DisulfideBonds.Last().OneBasedBeginPosition);
            Assert.AreEqual(originalProtein.DisulfideBonds.Last().OneBasedEndPosition, proteinReadFromXml[0].DisulfideBonds.Last().OneBasedEndPosition);

            Assert.AreEqual(originalProtein.FullDescription, proteinReadFromXml[0].FullDescription);
            Assert.AreEqual(originalProtein.FullName, proteinReadFromXml[0].FullName);
            Assert.AreEqual(originalProtein.GeneNames, proteinReadFromXml[0].GeneNames);
            Assert.AreEqual(originalProtein.IsContaminant, proteinReadFromXml[0].IsContaminant);
            Assert.AreEqual(originalProtein.IsDecoy, proteinReadFromXml[0].IsDecoy);
            Assert.AreEqual(originalProtein.Length, proteinReadFromXml[0].Length);
            Assert.AreEqual(originalProtein.Name, proteinReadFromXml[0].Name);
            Assert.AreEqual(originalProtein.Organism, proteinReadFromXml[0].Organism);
            Assert.AreEqual(originalProtein.DatabaseFilePath, proteinReadFromXml[0].DatabaseFilePath);
            Assert.AreEqual(1, originalProtein.OneBasedPossibleLocalizedModifications.Keys.Count);
            Assert.AreEqual(1, proteinReadFromXml[0].OneBasedPossibleLocalizedModifications.Keys.Count);
            Assert.AreEqual(originalProtein.OneBasedPossibleLocalizedModifications.Keys.First(), proteinReadFromXml[0].OneBasedPossibleLocalizedModifications.Keys.First());
            Assert.IsTrue(originalProtein.OneBasedPossibleLocalizedModifications[5][0].Equals(proteinReadFromXml[0].OneBasedPossibleLocalizedModifications[5][0]));

            Assert.AreEqual(originalProtein.TruncationProducts.First().OneBasedBeginPosition, proteinReadFromXml[0].TruncationProducts.First().OneBasedBeginPosition);
            Assert.AreEqual(originalProtein.TruncationProducts.First().OneBasedEndPosition, proteinReadFromXml[0].TruncationProducts.First().OneBasedEndPosition);
            Assert.AreEqual(originalProtein.TruncationProducts.First().Type, proteinReadFromXml[0].TruncationProducts.First().Type.Split('(')[0]);

            Assert.AreEqual(originalProtein.SequenceVariations.First().Description, proteinReadFromXml[0].SequenceVariations.First().Description);
            Assert.AreEqual(originalProtein.SequenceVariations.First().OneBasedBeginPosition, proteinReadFromXml[0].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(originalProtein.SequenceVariations.First().OneBasedEndPosition, proteinReadFromXml[0].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreEqual(originalProtein.SequenceVariations.First().OriginalSequence, proteinReadFromXml[0].SequenceVariations.First().OriginalSequence);
            Assert.AreEqual(originalProtein.SequenceVariations.First().VariantSequence, proteinReadFromXml[0].SequenceVariations.First().VariantSequence);
            Assert.AreEqual(originalProtein.SequenceVariations.Last().Description, proteinReadFromXml[0].SequenceVariations.Last().Description);
            Assert.AreEqual(originalProtein.SequenceVariations.Last().OneBasedBeginPosition, proteinReadFromXml[0].SequenceVariations.Last().OneBasedBeginPosition);
            Assert.AreEqual(originalProtein.SequenceVariations.Last().OneBasedEndPosition, proteinReadFromXml[0].SequenceVariations.Last().OneBasedEndPosition);
            Assert.AreEqual(originalProtein.SequenceVariations.Last().OriginalSequence, proteinReadFromXml[0].SequenceVariations.Last().OriginalSequence);
            Assert.AreEqual(originalProtein.SequenceVariations.Last().VariantSequence, proteinReadFromXml[0].SequenceVariations.Last().VariantSequence);
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
            Assert.AreEqual(ok[0].SequenceVariations.First().Description, ok2[0].SequenceVariations.First().Description);
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
            Assert.AreEqual(ok[0].SequenceVariations.First().Description, ok2[0].SequenceVariations.First().Description);
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
    }
}