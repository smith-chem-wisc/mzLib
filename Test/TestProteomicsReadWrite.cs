using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class TestProteomicsReadWrite
    {

        #region Public Methods

        [Test]
        public void test_read_write_read_xml()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, nice, false, new List<string> { "GO" }, null, out un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"), false, nice, false, new List<string> { "GO" }, null, out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(9, ok[0].DatabaseReferences.Count());
            Assert.AreEqual(3, ok[0].DatabaseReferences.First().Properties.Count());
            Assert.AreEqual(3, ok[0].GeneNames.Count());
            Assert.AreEqual("primary", ok[0].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", ok[0].GeneNames.First().Item2);

            Assert.AreEqual(9, ok2[0].DatabaseReferences.Count());
            Assert.AreEqual(3, ok2[0].DatabaseReferences.First().Properties.Count());
            Assert.AreEqual(3, ok2[0].GeneNames.Count());
            Assert.AreEqual("primary", ok2[0].GeneNames.First().Item1);
            Assert.AreEqual("JJJ1", ok2[0].GeneNames.First().Item2);

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void test_read_Ensembl_pepAllFasta()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_gene_expression);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.xml"), false, nice, false, null, null, out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual("ENSP00000381386", ok[0].Accession);
            Assert.AreEqual("ENSP00000215773", ok[1].Accession);
            Assert.AreEqual("ENSG00000099977", ok[0].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok[1].GeneNames.First().Item2);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[0].FullName);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[1].FullName);

            Assert.AreEqual("ENSP00000381386", ok2[0].Accession);
            Assert.AreEqual("ENSP00000215773", ok2[1].Accession);
            Assert.AreEqual("ENSG00000099977", ok2[0].GeneNames.First().Item2);
            Assert.AreEqual("ENSG00000099977", ok2[1].GeneNames.First().Item2);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[0].FullName);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[1].FullName);

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void test_read_write_read_fasta()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_gene_expression);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_gene_expression);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));

            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedBeginPosition == null || prod.OneBasedBeginPosition > 0 && prod.OneBasedBeginPosition <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.All(prod => prod.OneBasedEndPosition == null || prod.OneBasedEndPosition > 0 && prod.OneBasedEndPosition <= p.Length)));
        }

        [Test]
        public void test_accession_regex_weird()
        {
            Regex bad = new Regex(@"/()/");
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_ensembl.pep.all.fasta"), false, false, bad, bad, bad, bad);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), false, false, bad, bad, bad, bad);

            Assert.AreEqual("ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[0].Accession);
            Assert.AreEqual("ENSP00000381386 pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok2[0].Accession);
            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
        }

        [Test]
        public void test_write_with_custom_mods()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            ModificationMotif motif;
            ModificationMotif.TryGetMotif("K", out motif);
            ModificationWithMass m = new ModificationWithMass("mod", new Tuple<string, string>("", ""), motif, ModificationSites.Any, 1, new Dictionary<string, IList<string>>(), -1, new List<double>(), new List<double>(), "");

            Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> new_mods = new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>
            {
                {  "P53863", new HashSet<Tuple<int, ModificationWithMass>> {new Tuple<int, ModificationWithMass>(2, m ) } }
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, nice, false, null, null, out un);
            ProteinDbWriter.WriteXmlDatabase(new_mods, ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"), false, nice, false, null, null, out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual(2, ok[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(3, ok2[0].OneBasedPossibleLocalizedModifications.Count);
        }

        [Test]
        public void AnotherTest()
        {
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], "name1", "fullname1", false, false, new List<DatabaseReference>());

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new int?[] { 4 }, new int?[] { 8 }, new string[] { "chain" }, "name2", "fullname2", false, false, new List<DatabaseReference>());

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { ParentProtein, proteinWithChain }, Path.Combine(TestContext.CurrentContext.TestDirectory, @"fdsfsd.xml"));
        }

        #endregion Public Methods

    }
}