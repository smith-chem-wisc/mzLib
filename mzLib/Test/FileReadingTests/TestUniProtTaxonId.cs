using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The NCBI taxonomy identifier must survive loading, from EITHER database format.
    ///
    /// The XML path already kept it, and not by design for taxonomy specifically: ProteinXmlEntry
    /// stores every dbReference generically, so the one inside the &lt;organism&gt; block comes
    /// along for free. The FASTA path did not -- OX= appeared in two header regexes but only ever
    /// as a terminator for the preceding field, so the same protein loaded from a FASTA carried no
    /// taxonomy at all. These tests pin both paths to the same representation.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestUniProtTaxonId
    {
        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        [Test]
        public void Xml_RetainsTheTaxonomyId()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Data("DatabaseTests", "disulfidetests.xml"),
                true, DecoyType.None, null, false, null, out _);

            var taxon = proteins.First().DatabaseReferences
                .FirstOrDefault(r => r.Type == ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType);

            Assert.That(proteins.First().Organism, Is.EqualTo("Mus musculus"));
            Assert.That(taxon, Is.Not.Null);
            Assert.That(taxon!.Id, Is.EqualTo("10090"));
        }

        [Test]
        public void Fasta_RetainsTheTaxonomyIdFromOx()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"taxon_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(path,
                    ">sp|P12345|TEST_HUMAN Test protein OS=Homo sapiens OX=9606 GN=TST PE=1 SV=2\n" +
                    "PEPTIDEK\n");

                var proteins = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out _);

                var protein = proteins.Single();
                Assert.That(protein.Organism, Is.EqualTo("Homo sapiens"));

                var taxon = protein.DatabaseReferences
                    .FirstOrDefault(r => r.Type == ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType);
                Assert.That(taxon, Is.Not.Null, "OX= was previously only a delimiter, never captured");
                Assert.That(taxon!.Id, Is.EqualTo("9606"));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void Fasta_WithoutOx_HasNoTaxonomyReferenceAndStillLoads()
        {
            // Plenty of FASTA headers predate OX= or are not UniProt at all. Absence must be
            // absence, not an empty or zero taxon.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"notaxon_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(path,
                    ">sp|P12345|TEST_HUMAN Test protein OS=Homo sapiens GN=TST PE=1 SV=2\n" +
                    "PEPTIDEK\n");

                var protein = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out _).Single();

                Assert.That(protein.Organism, Is.EqualTo("Homo sapiens"));
                Assert.That(protein.DatabaseReferences
                    .Any(r => r.Type == ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType), Is.False);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void BothFormatsAgreeOnTheRepresentation()
        {
            // The point of the change: a consumer must not have to know which format a protein
            // came from in order to find its taxonomy.
            string fasta = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"agree_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(fasta,
                    ">sp|P99999|MOUSE_TEST Test OS=Mus musculus OX=10090 GN=TST PE=1 SV=1\nPEPTIDEK\n");

                var fromFasta = ProteinDbLoader.LoadProteinFasta(fasta, true, DecoyType.None, false, out _).Single();
                var fromXml = ProteinDbLoader.LoadProteinXML(Data("DatabaseTests", "disulfidetests.xml"),
                    true, DecoyType.None, null, false, null, out _).First();

                string TaxonOf(Proteomics.Protein p) => p.DatabaseReferences
                    .First(r => r.Type == ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType).Id;

                Assert.That(TaxonOf(fromFasta), Is.EqualTo("10090"));
                Assert.That(TaxonOf(fromXml), Is.EqualTo("10090"));
            }
            finally
            {
                if (File.Exists(fasta)) File.Delete(fasta);
            }
        }
    }
}
