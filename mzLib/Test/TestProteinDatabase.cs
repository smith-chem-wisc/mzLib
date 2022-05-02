﻿using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    [TestFixture]
    public sealed class TestProteinDatabase
    {
        [Test]
        public static void MakeAnewProteinWithAndWithoutTruncations()
        {
            Protein noTruncationProtein1 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addTruncations: false);
            Assert.AreEqual(0, noTruncationProtein1.ProteolysisProducts.Count());

            noTruncationProtein1.AddIntactProteoformToTruncationsProducts(7);
            Assert.AreEqual(1, noTruncationProtein1.ProteolysisProducts.Count());

            Protein noTruncationProtein2 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addTruncations: false);
            Assert.AreEqual(0, noTruncationProtein2.ProteolysisProducts.Count());

            noTruncationProtein2.AddIntactProteoformToTruncationsProducts(7);
            Assert.AreEqual(1, noTruncationProtein2.ProteolysisProducts.Count());

            Protein noTruncationProtein3 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addTruncations: false);
            Assert.AreEqual(0, noTruncationProtein3.ProteolysisProducts.Count());

            noTruncationProtein3.AddIntactProteoformToTruncationsProducts(7);
            Assert.AreEqual(1, noTruncationProtein3.ProteolysisProducts.Count());

            Protein truncationProtein1 = new("PEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addTruncations: true);
            Assert.AreEqual(11, truncationProtein1.ProteolysisProducts.Count());

            Protein truncationProtein2 = new("PEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addTruncations: false);
            truncationProtein2.AddIntactProteoformToTruncationsProducts(7);
            Assert.AreEqual(1, truncationProtein2.ProteolysisProducts.Count());
        }

        [Test]
        public static void AddTruncationsToProteolysisProducts()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            Protein insulinProteinFromXml1
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications1, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml1.ProteolysisProducts.Count());
            insulinProteinFromXml1.AddTruncationsToExistingProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, 7, 5, "truncation");
            Assert.AreEqual(20, insulinProteinFromXml1.ProteolysisProducts.Count());

            Protein insulinProteinFromXml2
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications2, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml2.ProteolysisProducts.Count());
            insulinProteinFromXml2.AddTruncationsToExistingProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, 7, 5, "truncation");
            Assert.AreEqual(20, insulinProteinFromXml2.ProteolysisProducts.Count());

            Protein insulinProteinFromXml3
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications3, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml3.ProteolysisProducts.Count());
            insulinProteinFromXml3.AddTruncationsToExistingProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, 7, 5, "truncation");
            Assert.AreEqual(20, insulinProteinFromXml3.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestRemoveMethionineWhenAppropriate()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");

            Protein insulinProteinFromXml1
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications1, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml1.ProteolysisProducts.Count());

            Protein insulinProteinFromXml2
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications2, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml2.ProteolysisProducts.Count());

            Protein insulinProteinFromXml3
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications3, addTruncations: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml3.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestAddTruncationsIntactAndExistingProteolysisProducts()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addTruncations: true)[0];

            Assert.AreEqual(17, insulinProteinFromFasta.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinProteinFromFasta.ProteolysisProducts.Where(p => p.Type == "full-length proteoform").Count());
            Assert.AreEqual(16, insulinProteinFromFasta.ProteolysisProducts.Where(p => p.Type.Contains("truncation")).Count());

            List<int> expectedBegins = new() { 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1 };
            List<int> expectedEnds = new() { 110, 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 109, 108, 107, 106, 105 };

            CollectionAssert.AreEquivalent(expectedBegins, insulinProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedBeginPosition).ToList());
            CollectionAssert.AreEquivalent(expectedEnds, insulinProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedEndPosition).ToList());

            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            Protein insulinProteinFromXml
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications, addTruncations: true)[0];

            Assert.AreEqual(68, insulinProteinFromXml.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinProteinFromXml.ProteolysisProducts.Where(p => p.Type == "full-length proteoform").Count());
            Assert.AreEqual(62, insulinProteinFromXml.ProteolysisProducts.Where(p => p.Type.Contains("truncation")).Count()); //4 are original proteolysis products

            expectedBegins = new List<int> { 1, 25, 57, 90, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 26, 27, 28, 29, 30, 25, 25, 25, 25, 25, 58, 59, 60, 61, 62, 57, 57, 57, 57, 57, 91, 92, 93, 94, 95, 90, 90, 90, 90, 90, 25 };
            expectedEnds = new List<int> { 24, 54, 87, 110, 110, 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 109, 108, 107, 106, 105, 24, 24, 24, 24, 24, 24, 23, 22, 21, 20, 19, 23, 22, 21, 20, 19, 54, 54, 54, 54, 54, 53, 52, 51, 50, 49, 87, 87, 87, 87, 87, 86, 85, 84, 83, 82, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 110 };

            List<int> reportedBegins = insulinProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedBeginPosition.Value).ToList();
            List<int> reportedEnds = insulinProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedEndPosition.Value).ToList();
            CollectionAssert.AreEquivalent(expectedBegins, reportedBegins);
            CollectionAssert.AreEquivalent(expectedEnds, reportedEnds);
        }

        [Test]
        public static void TestMethionineCleave()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addTruncations: false)[0];

            Assert.AreEqual(0, insulinProteinFromFasta.ProteolysisProducts.Count());
            insulinProteinFromFasta.AddTruncations();
            Assert.AreEqual(17, insulinProteinFromFasta.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestMethionineCleaveNoMethionine()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addTruncations: false)[0];

            Protein noMethionine = new(insulinProteinFromFasta.BaseSequence.Substring(1,insulinProteinFromFasta.BaseSequence.Length-1), insulinProteinFromFasta.Accession);

            Assert.AreEqual(0, noMethionine.ProteolysisProducts.Count());
            noMethionine.AddTruncations();
            Assert.AreEqual(11, noMethionine.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestMethionineVariable()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addTruncations: false)[0];

            Assert.AreEqual(0, insulinProteinFromFasta.ProteolysisProducts.Count());
            insulinProteinFromFasta.AddTruncations();
            Assert.AreEqual(17, insulinProteinFromFasta.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestDoNotWriteTruncationsToXml()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "TestProtein.xml");
            List<Protein> proteins
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, addTruncations: true);

            Assert.AreEqual(16, proteins[0].ProteolysisProducts.Where(p => p.Type.Contains("truncation")).Count());

            string testOutXml = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "testOutXml.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), testOutXml);
            string[] lines = File.ReadAllLines(testOutXml);
            Assert.AreEqual(0, lines.Where(l => l.Contains("truncation")).Count());

            List<Protein> moreProteins
                = ProteinDbLoader.LoadProteinXML(testOutXml, true,
                DecoyType.Reverse, null, false, null, out var moreUnknownModifications, addTruncations: false);
            Assert.AreEqual(0, moreProteins[0].ProteolysisProducts.Where(p => p.Type.Contains("truncation")).Count());

            File.Delete(testOutXml);
        }
    }
}