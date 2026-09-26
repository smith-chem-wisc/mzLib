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

        /// <summary>
        /// Regression: decoy proteins did not inherit the taxonomy their targets carry, so half of
        /// a decoy-containing database answered the organism question and half did not.
        /// </summary>
        [Test]
        public void Decoys_InheritTheTaxonomyOfTheirTarget()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"decoytaxon_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(path,
                    ">sp|P12345|TEST_HUMAN Test OS=Homo sapiens OX=9606 GN=TST PE=1 SV=2\nPEPTIDEKPEPTIDER\n");

                var proteins = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.Reverse, false, out _);

                var target = proteins.Single(p => !p.IsDecoy);
                var decoy = proteins.Single(p => p.IsDecoy);

                Assert.That(target.NcbiTaxonomyId, Is.EqualTo("9606"));
                Assert.That(decoy.NcbiTaxonomyId, Is.EqualTo("9606"),
                    "a decoy is a reversed sequence from this organism's database");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void Decoys_DoNotInheritUnrelatedDatabaseReferences()
        {
            // Only the taxonomy travels. A decoy is not annotated in EMBL, GO or anywhere else, and
            // claiming otherwise would assert things about it that are false.
            var proteins = ProteinDbLoader.LoadProteinXML(Data("DatabaseTests", "disulfidetests.xml"),
                true, DecoyType.Reverse, null, false, null, out _);

            var decoy = proteins.First(p => p.IsDecoy);
            Assert.That(decoy.NcbiTaxonomyId, Is.EqualTo("10090"));
            Assert.That(decoy.DatabaseReferences.Select(r => r.Type).Distinct(),
                Is.EquivalentTo(new[] { ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType }));
        }

        /// <summary>
        /// Regression: LoadProteinFasta retained the taxonomy but GetUniProtFastaHeader did not
        /// write it, so a FASTA-to-FASTA round trip silently dropped it again.
        /// </summary>
        [Test]
        public void FastaHeader_RoundTripsTheTaxonomyId()
        {
            string input = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"rt_in_{System.Guid.NewGuid():N}.fasta");
            string output = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"rt_out_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(input,
                    ">sp|P12345|TEST_HUMAN Test OS=Homo sapiens OX=9606 GN=TST PE=1 SV=2\nPEPTIDEK\n");

                var loaded = ProteinDbLoader.LoadProteinFasta(input, true, DecoyType.None, false, out _);
                ProteinDbWriter.WriteFastaDatabase(loaded.ToList(), output, "|");

                string header = File.ReadAllLines(output).First(l => l.StartsWith(">"));
                Assert.That(header, Does.Contain("OX=9606"));

                var reloaded = ProteinDbLoader.LoadProteinFasta(output, true, DecoyType.None, false, out _);
                Assert.That(reloaded.Single().NcbiTaxonomyId, Is.EqualTo("9606"));
                Assert.That(reloaded.Single().Organism, Is.EqualTo("Homo sapiens"));
            }
            finally
            {
                foreach (var f in new[] { input, output })
                    if (File.Exists(f)) File.Delete(f);
            }
        }

        [Test]
        public void FastaHeader_OmitsOxWhenTheTaxonomyIsUnknown()
        {
            var protein = new Proteomics.Protein("PEPTIDEK", "P00001", organism: "Homo sapiens");
            Assert.That(protein.NcbiTaxonomyId, Is.Null);
            Assert.That(protein.GetUniProtFastaHeader(), Does.Not.Contain("OX="));
        }

        [Test]
        public void OxRegex_IsAnchoredAndDoesNotMatchATrailingToken()
        {
            // Unanchored, a description token ending in "OX=" followed by digits captured the wrong
            // number. A taxonomy id is a join key: a silently wrong value is worse than none.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"anchor_{System.Guid.NewGuid():N}.fasta");
            try
            {
                File.WriteAllText(path,
                    ">sp|P12345|TEST_HUMAN Protein SOX=5 family OS=Homo sapiens OX=9606 GN=TST PE=1 SV=2\nPEPTIDEK\n");

                var protein = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out _).Single();
                Assert.That(protein.NcbiTaxonomyId, Is.EqualTo("9606"), "SOX=5 must not be read as the taxonomy");
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
