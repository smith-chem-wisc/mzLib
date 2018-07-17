using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class TestProteomicsReadWrite
    {
        [Test]
        public void ReadXmlNulls()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None,
                null, false, null, out Dictionary<string, ModificationGeneral> un);
        }

        [Test]
        public void Test_read_write_read_xml()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, nice, false, null,
                out Dictionary<string, ModificationGeneral> un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml2.xml"));
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
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void Test_read_Ensembl_pepAllFasta()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out var a);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.xml"), true, DecoyType.None, nice,
                false, null, out Dictionary<string, ModificationGeneral> un);

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

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void Test_read_write_read_fasta()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out var a);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out var b);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void Test_read_xml_write_read_fasta()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, nice, false, null,
                out Dictionary<string, ModificationGeneral> un);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml_test.fasta"), "|");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml_test.fasta"), true, DecoyType.None, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, out var b);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].Name == ok2[i].Name));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].Organism == ok2[i].Organism));
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].GeneNames.First().Item2 == ok2[i].GeneNames.First().Item2));

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void Test_accession_regex_weird()
        {
            FastaHeaderFieldRegex bad = new FastaHeaderFieldRegex("", @"/()/", 0, 1);
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"test_ensembl.pep.all.fasta"), true, DecoyType.None, false,
                bad, bad, bad, bad, bad, out var a);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_test_ensembl.pep.all.fasta"), true, DecoyType.None, false,
                bad, bad, bad, bad, bad, out var b);

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

            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motiff, "Anywhere.", null, null, null, null, null, null, null, null),
                new ModificationGeneral("Phosphoserine", null, "mt", null, m1, "Anywhere.", null, null, null, null, null, null, null, null),
                new ModificationGeneral("Phosphothreonine", null, "mt", null,  m2, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            ModificationGeneral m = new ModificationGeneral("mod", null, "mt", null, motif, "Anywhere.", null, 1, null, null, null, new Dictionary<DissociationType, List<double>>() { { DissociationType.AnyActivationType, new List<double> { -1 } } }, null, null);

            Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>> new_mods = new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>
            {
                {  "P53863", new HashSet<Tuple<int, ModificationGeneral>> {new Tuple<int, ModificationGeneral>(2, m ) } }
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"), true, DecoyType.None, nice, false, new List<string>(),
                out Dictionary<string, ModificationGeneral> un);
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
            List<ModificationGeneral> variableModifications = new List<ModificationGeneral>();
            List<ModificationGeneral> fixedModifications = new List<ModificationGeneral>();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<ModificationGeneral>>(), null,
                "name1", "fullname1", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), new List<DisulfideBond>());

            List<ProteolysisProduct> pp = new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") };
            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<ModificationGeneral>>(), pp,
                "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), new List<DisulfideBond>());

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), new List<Protein> { ParentProtein, proteinWithChain }, Path.Combine(TestContext.CurrentContext.TestDirectory, @"fdsfsd.xml"));
        }

        [Test]
        public void TestEmptyProteins()
        {
            Protein p1 = new Protein("SEQENCE", "p1");
            Assert.AreEqual("p1||", p1.FullDescription);
            Protein p2 = new Protein("SEQENCE", "p2", name: "namep2");

            var proteinListToWrite = new List<Protein> { p1, p2 };

            // Generate data for files
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), proteinListToWrite,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"differentlyConstuctedProteins.xml"));

            IEnumerable<string> modTypesToExclude = new List<string>();
            IEnumerable<ModificationGeneral> allKnownModifications = new List<ModificationGeneral>();
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"differentlyConstuctedProteins.xml"), true, DecoyType.None,
                allKnownModifications, false, modTypesToExclude, out Dictionary<string, ModificationGeneral> un);
            Assert.AreEqual(p1.Accession, ok[0].Accession);
            Assert.AreEqual(p2.Accession, ok[1].Accession);
            Assert.AreEqual(p1.Name, ok[0].Name);
            Assert.AreEqual(p2.Name, ok[1].Name);
        }

        [Test]
        public void TestFullProteinReadWrite()
        {
            ModificationGeneral mod = new ModificationGeneral("mod1", null, "modType1", null, null, null, null, null, null, null, null, null, null, null);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            ModificationGeneral mod2 = new ModificationGeneral("mod2", null, "modType1", null, motif, "Anywhere.", null, null, null, null, null, null, null, null);
            ModificationMotif.TryGetMotif("N", out ModificationMotif motif3);
            ModificationGeneral mod3 = new ModificationGeneral("mod3", null, "modType1", null, motif3, "Anywhere.", null, 10, null, null, null, null, null, null);

            List<Tuple<string, string>> gene_names = new List<Tuple<string, string>> { new Tuple<string, string>("a", "b") };
            IDictionary<int, List<ModificationGeneral>> oneBasedModifications = new Dictionary<int, List<ModificationGeneral>>
            {
                {3, new List<ModificationGeneral>{mod} },
                {4, new List<ModificationGeneral>{mod2} },
                {5, new List<ModificationGeneral>{mod3} }
            };
            List<ProteolysisProduct> proteolysisProducts = new List<ProteolysisProduct> { new ProteolysisProduct(1, 2, "propeptide") };

            string name = "testName";

            string full_name = "testFullName";

            List<DatabaseReference> databaseReferences = new List<DatabaseReference> {
                new DatabaseReference("type1", "id1", new List<Tuple<string, string>> { new Tuple<string, string>("e1", "e2") }) };

            List<SequenceVariation> sequenceVariations = new List<SequenceVariation> { new SequenceVariation(3,"Q", "N", "replace Q by N"),
            new SequenceVariation(3,4,"QE", "NN", "replace QE by NN")};

            List<DisulfideBond> disulfideBonds = new List<DisulfideBond> { new DisulfideBond(1, "ds1"), new DisulfideBond(2, 3, "ds2") };

            Protein p1 = new Protein(
                "SEQENCE",
                "a1",
                gene_names: gene_names,
                oneBasedModifications: oneBasedModifications,
                proteolysisProducts: proteolysisProducts,
                name: name,
                full_name: full_name,
                isDecoy: false,
                isContaminant: true,
                databaseReferences: databaseReferences,
                sequenceVariations: sequenceVariations,
                disulfideBonds: disulfideBonds,
                databaseFilePath: Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"));

            // Generate data for files
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), new List<Protein> { p1 },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"));

            IEnumerable<string> modTypesToExclude = new List<string>();
            IEnumerable<ModificationGeneral> allKnownModifications = new List<ModificationGeneral>();
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bnueiwhf.xml"), true, DecoyType.None,
                allKnownModifications, true, modTypesToExclude, out Dictionary<string, ModificationGeneral> unknownModifications);
            Assert.AreEqual(p1.Accession, ok[0].Accession);
            Assert.AreEqual(p1.BaseSequence, ok[0].BaseSequence);
            Assert.AreEqual(p1.DatabaseReferences.First().Id, ok[0].DatabaseReferences.First().Id);
            Assert.AreEqual(p1.DatabaseReferences.First().Properties.First().Item1, ok[0].DatabaseReferences.First().Properties.First().Item1);
            Assert.AreEqual(p1.DatabaseReferences.First().Properties.First().Item2, ok[0].DatabaseReferences.First().Properties.First().Item2);
            Assert.AreEqual(p1.DatabaseReferences.First().Type, ok[0].DatabaseReferences.First().Type);

            Assert.AreEqual(p1.DisulfideBonds.First().Description, ok[0].DisulfideBonds.First().Description);
            Assert.AreEqual(p1.DisulfideBonds.First().OneBasedBeginPosition, ok[0].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(p1.DisulfideBonds.First().OneBasedEndPosition, ok[0].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(p1.DisulfideBonds.Last().Description, ok[0].DisulfideBonds.Last().Description);
            Assert.AreEqual(p1.DisulfideBonds.Last().OneBasedBeginPosition, ok[0].DisulfideBonds.Last().OneBasedBeginPosition);
            Assert.AreEqual(p1.DisulfideBonds.Last().OneBasedEndPosition, ok[0].DisulfideBonds.Last().OneBasedEndPosition);

            Assert.AreEqual(p1.FullDescription, ok[0].FullDescription);
            Assert.AreEqual(p1.FullName, ok[0].FullName);
            Assert.AreEqual(p1.GeneNames, ok[0].GeneNames);
            Assert.AreEqual(p1.IsContaminant, ok[0].IsContaminant);
            Assert.AreEqual(p1.IsDecoy, ok[0].IsDecoy);
            Assert.AreEqual(p1.Length, ok[0].Length);
            Assert.AreEqual(p1.Name, ok[0].Name);
            Assert.AreEqual(p1.Organism, ok[0].Organism);
            Assert.AreEqual(p1.DatabaseFilePath, ok[0].DatabaseFilePath);
            Assert.AreEqual(p1.OneBasedPossibleLocalizedModifications[3][0], ok[0].OneBasedPossibleLocalizedModifications[3][0]);
            Assert.AreEqual(p1.OneBasedPossibleLocalizedModifications[3][0].Id, ok[0].OneBasedPossibleLocalizedModifications[3][0].Id);
            Assert.AreEqual(p1.OneBasedPossibleLocalizedModifications[3][0].ModificationType, ok[0].OneBasedPossibleLocalizedModifications[3][0].ModificationType);

            Assert.AreEqual(p1.OneBasedPossibleLocalizedModifications[4][0].Id, ok[0].OneBasedPossibleLocalizedModifications[4][0].Id);
            Assert.AreEqual(p1.OneBasedPossibleLocalizedModifications[4][0].ModificationType, ok[0].OneBasedPossibleLocalizedModifications[4][0].ModificationType);
            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).DatabaseReference, (ok[0].OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).DatabaseReference);
            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).Target, (ok[0].OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).Target);
            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).LocationRestriction, (ok[0].OneBasedPossibleLocalizedModifications[4][0] as ModificationGeneral).LocationRestriction);

            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).DiagnosticIons, (ok[0].OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).DiagnosticIons);
            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).NeutralLosses, (ok[0].OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).NeutralLosses);
            Assert.AreEqual((p1.OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).MonoisotopicMass, (ok[0].OneBasedPossibleLocalizedModifications[5][0] as ModificationGeneral).MonoisotopicMass);

            Assert.AreEqual(p1.ProteolysisProducts.First().OneBasedBeginPosition, ok[0].ProteolysisProducts.First().OneBasedBeginPosition);
            Assert.AreEqual(p1.ProteolysisProducts.First().OneBasedEndPosition, ok[0].ProteolysisProducts.First().OneBasedEndPosition);
            Assert.AreEqual(p1.ProteolysisProducts.First().Type, ok[0].ProteolysisProducts.First().Type);

            Assert.AreEqual(p1.SequenceVariations.First().Description, ok[0].SequenceVariations.First().Description);
            Assert.AreEqual(p1.SequenceVariations.First().OneBasedBeginPosition, ok[0].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(p1.SequenceVariations.First().OneBasedEndPosition, ok[0].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreEqual(p1.SequenceVariations.First().OriginalSequence, ok[0].SequenceVariations.First().OriginalSequence);
            Assert.AreEqual(p1.SequenceVariations.First().VariantSequence, ok[0].SequenceVariations.First().VariantSequence);
            Assert.AreEqual(p1.SequenceVariations.Last().Description, ok[0].SequenceVariations.Last().Description);
            Assert.AreEqual(p1.SequenceVariations.Last().OneBasedBeginPosition, ok[0].SequenceVariations.Last().OneBasedBeginPosition);
            Assert.AreEqual(p1.SequenceVariations.Last().OneBasedEndPosition, ok[0].SequenceVariations.Last().OneBasedEndPosition);
            Assert.AreEqual(p1.SequenceVariations.Last().OriginalSequence, ok[0].SequenceVariations.Last().OriginalSequence);
            Assert.AreEqual(p1.SequenceVariations.Last().VariantSequence, ok[0].SequenceVariations.Last().VariantSequence);
        }

        [Test]
        public void TestReadWriteSeqVars()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), true, DecoyType.None,
                nice, false, null, out Dictionary<string, ModificationGeneral> un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_xml.xml"));
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
            var nice = new List<ModificationGeneral>
            {
                new ModificationGeneral("fayk", null, "mt", null, motif, "Anywhere.", null, null, null, null, null, null, null, null)
            };

            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"seqvartests.xml"), true, DecoyType.None,
                nice, false, new List<string>(), out Dictionary<string, ModificationGeneral> un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationGeneral>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"rewrite_seqvartests.xml"));
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
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt")).ToList();
            char[] myChar = { '"' };
            string output = a.First().ToString();
            Assert.AreEqual(output.TrimStart(myChar).TrimEnd(myChar), "ID   4-carboxyglutamate\r\nMT   Biological\r\nTG   E\r\nPP   Anywhere.\r\nCF   CO2\r\nMM   43.989829\r\n");
        }

        [Test]
        public void TestModificationGeneral_Equals()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt")).ToList();
            var b = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "CommonBiological.txt")).ToList();

            Assert.IsTrue(a.First().Equals(b.First()));
        }
    }
}