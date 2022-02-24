using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestProteinDatabase
    {
        [Test]
        public static void MakeAnewProteinWithAndWithoutBiomarkers()
        {
            Protein noBiomarkerProtein1 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addBiomarkers: false);
            Assert.AreEqual(0, noBiomarkerProtein1.ProteolysisProducts.Count());

            noBiomarkerProtein1.AddIntactProteoformToProteolysisProducts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Cleave, 7);
            Assert.AreEqual(1, noBiomarkerProtein1.ProteolysisProducts.Count());

            Protein noBiomarkerProtein2 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addBiomarkers: false);
            Assert.AreEqual(0, noBiomarkerProtein2.ProteolysisProducts.Count());

            noBiomarkerProtein2.AddIntactProteoformToProteolysisProducts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Retain, 7);
            Assert.AreEqual(1, noBiomarkerProtein2.ProteolysisProducts.Count());

            Protein noBiomarkerProtein3 = new("MPEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addBiomarkers: false);
            Assert.AreEqual(0, noBiomarkerProtein3.ProteolysisProducts.Count());

            noBiomarkerProtein3.AddIntactProteoformToProteolysisProducts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Variable, 7);
            Assert.AreEqual(1, noBiomarkerProtein3.ProteolysisProducts.Count());

            Protein biomarkerProtein1 = new("PEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addBiomarkers: true);
            Assert.AreEqual(11, biomarkerProtein1.ProteolysisProducts.Count());

            Protein biomarkerProtein2 = new("PEPTIDEPEPTIDEPEPTIDE", "ACCESSION", addBiomarkers: false);
            biomarkerProtein2.AddIntactProteoformToProteolysisProducts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Cleave, 7);
            Assert.AreEqual(1, biomarkerProtein2.ProteolysisProducts.Count());
        }

        [Test]
        public static void AddBiomarkersToProteolysisProducts()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            Protein insulinProteinFromXml1
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications1, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml1.ProteolysisProducts.Count());
            insulinProteinFromXml1.AddBiomarkersToProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Retain, 7, 5, "biomarker");
            Assert.AreEqual(14, insulinProteinFromXml1.ProteolysisProducts.Count());

            Protein insulinProteinFromXml2
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications2, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml2.ProteolysisProducts.Count());
            insulinProteinFromXml2.AddBiomarkersToProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Cleave, 7, 5, "biomarker");
            Assert.AreEqual(14, insulinProteinFromXml2.ProteolysisProducts.Count());


            Protein insulinProteinFromXml3
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications3, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml3.ProteolysisProducts.Count());
            insulinProteinFromXml3.AddBiomarkersToProteolysisProducts(1, insulinProteinFromXml1.BaseSequence.Length, true, true, Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Variable, 7, 5, "biomarker");
            Assert.AreEqual(15, insulinProteinFromXml3.ProteolysisProducts.Count());
        }

        [Test]
        public static void TestRemoveMethionineWhenAppropriate()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");

            Protein insulinProteinFromXml1
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications1, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml1.ProteolysisProducts.Count());
            insulinProteinFromXml1.RemoveMethionineWhenAppropriateFromExistingProduts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Retain);
            Assert.AreEqual(1, insulinProteinFromXml1.ProteolysisProducts.First().OneBasedBeginPosition.Value);

            Protein insulinProteinFromXml2
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications2, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml2.ProteolysisProducts.Count());
            insulinProteinFromXml2.RemoveMethionineWhenAppropriateFromExistingProduts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Cleave);
            Assert.AreEqual(2, insulinProteinFromXml2.ProteolysisProducts.ToList()[3].OneBasedBeginPosition.Value);


            Protein insulinProteinFromXml3
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications3, addBiomarkers: false)[0];

            Assert.AreEqual(4, insulinProteinFromXml3.ProteolysisProducts.Count());
            insulinProteinFromXml3.RemoveMethionineWhenAppropriateFromExistingProduts(Proteomics.ProteolyticDigestion.InitiatorMethionineBehavior.Variable);
            Assert.AreEqual(1, insulinProteinFromXml3.ProteolysisProducts.ToList()[0].OneBasedBeginPosition.Value);
        }

        [Test]
        public static void TestAddBiomarkersIntactAndExistingProteolysisProducts()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addBiomarkers: true)[0];

            Assert.AreEqual(11, insulinProteinFromFasta.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinProteinFromFasta.ProteolysisProducts.Where(p => p.Type.Contains("intact")).Count());
            Assert.AreEqual(11, insulinProteinFromFasta.ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count());
            List<int> expectedBegins = new List<int> { 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1 };
            List<int> expectedEnds = new List<int> { 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105 };

            CollectionAssert.AreEquivalent(expectedBegins, insulinProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedBeginPosition).ToList());
            CollectionAssert.AreEquivalent(expectedEnds, insulinProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedEndPosition).ToList());

            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            Protein insulinDecoyProteinFromXml
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications, addBiomarkers: true)[0];

            Assert.AreEqual(56, insulinDecoyProteinFromXml.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinDecoyProteinFromXml.ProteolysisProducts.Where(p => p.Type.Contains("intact")).Count());
            Assert.AreEqual(51, insulinDecoyProteinFromXml.ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count()); //4 are original proteolysis products

            expectedBegins = new List<int> { 1, 25, 57, 90, 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 26, 27, 28, 29, 30, 25, 25, 25, 25, 25, 58, 59, 60, 61, 62, 57, 57, 57, 57, 57, 91, 92, 93, 94, 95, 90, 90, 90, 90, 90, 25 };
            expectedEnds = new List<int> { 24, 54, 87, 110, 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 24, 24, 24, 24, 24, 23, 22, 21, 20, 19, 54, 54, 54, 54, 54, 53, 52, 51, 50, 49, 87, 87, 87, 87, 87, 86, 85, 84, 83, 82, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 110 };

            List<int> reportedBegins = insulinDecoyProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedBeginPosition.Value).ToList();
            List<int> reportedEnds = insulinDecoyProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedEndPosition.Value).ToList();
            CollectionAssert.AreEquivalent(expectedBegins, reportedBegins);
            CollectionAssert.AreEquivalent(expectedEnds, reportedEnds);
        }

        [Test]
        public static void TestDecoyBiomarkersAreAdded()
        {
            //Note: existing proteoloysis products are now subjected to additional proteolysis.

            //with fasta (there are no existing proteolysis products. so we rely on the code to deal with that non-factor)
            string fastaDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.fasta");
            Protein insulinDecoyProteinFromFasta = ProteinDbLoader.LoadProteinFasta(fastaDatabase, true, DecoyType.Reverse, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, addBiomarkers: true)[1];

            Assert.AreEqual(11, insulinDecoyProteinFromFasta.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinDecoyProteinFromFasta.ProteolysisProducts.Where(p => p.Type.Contains("intact")).Count());
            Assert.AreEqual(11, insulinDecoyProteinFromFasta.ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count());
            List<int> expectedBegins = new() { 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1 };
            List<int> expectedEnds = new() { 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105 };

            CollectionAssert.AreEquivalent(expectedBegins, insulinDecoyProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedBeginPosition).ToList());
            CollectionAssert.AreEquivalent(expectedEnds, insulinDecoyProteinFromFasta.ProteolysisProducts.Select(p => p.OneBasedEndPosition).ToList());

            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            Protein insulinProteinFromXml
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, addBiomarkers: true)[1];

            Assert.AreEqual(56, insulinProteinFromXml.ProteolysisProducts.Count());
            Assert.AreEqual(1, insulinProteinFromXml.ProteolysisProducts.Where(p => p.Type.Contains("intact")).Count());
            Assert.AreEqual(51, insulinProteinFromXml.ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count()); //4 are original proteolysis products

            expectedBegins = new List<int> { 1, 25, 57, 90, 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 1, 1, 1, 1, 1, 26, 27, 28, 29, 30, 25, 25, 25, 25, 25, 58, 59, 60, 61, 62, 57, 57, 57, 57, 57, 91, 92, 93, 94, 95, 90, 90, 90, 90, 90, 25 };
            expectedEnds = new List<int> { 24, 54, 87, 110, 110, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 24, 24, 24, 24, 24, 23, 22, 21, 20, 19, 54, 54, 54, 54, 54, 53, 52, 51, 50, 49, 87, 87, 87, 87, 87, 86, 85, 84, 83, 82, 110, 110, 110, 110, 110, 109, 108, 107, 106, 105, 110 };

            List<int> reportedBegins = insulinProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedBeginPosition.Value).ToList();
            List<int> reportedEnds = insulinProteinFromXml.ProteolysisProducts.Select(p => p.OneBasedEndPosition.Value).ToList();
            CollectionAssert.AreEquivalent(expectedBegins, reportedBegins);
            CollectionAssert.AreEquivalent(expectedEnds, reportedEnds);
        }

        [Test]
        public static void TestDoNotWriteBiomarkersToXml()
        {
            //with xml, here for this protein, there are existing proteolysis products
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "TestProtein.xml");
            List<Protein> proteins
                = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, addBiomarkers: true);

            Assert.AreEqual(11, proteins[0].ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count());

            string testOutXml = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "testOutXml.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), testOutXml);
            string[] lines = File.ReadAllLines(testOutXml);
            Assert.AreEqual(0, lines.Where(l => l.Contains("biomarker")).Count());

            List<Protein> moreProteins
                = ProteinDbLoader.LoadProteinXML(testOutXml, true,
                DecoyType.Reverse, null, false, null, out var moreUnknownModifications, addBiomarkers: false);
            Assert.AreEqual(0, moreProteins[0].ProteolysisProducts.Where(p => p.Type.Contains("biomarker")).Count());

            File.Delete(testOutXml);
        }
    }
}