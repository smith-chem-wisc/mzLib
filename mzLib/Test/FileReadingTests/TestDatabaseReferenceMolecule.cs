using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// UniProt scopes some cross-references to one isoform with a &lt;molecule id="P04406-2"/&gt; child of the
    /// &lt;dbReference&gt;. In the UniProt 2026_03 gene-centric XMLs that is 610 of 47,539 human Reactome references,
    /// and tens of thousands of Ensembl, RefSeq and CCDS ones. Read without the molecule, an isoform-only pathway or
    /// transcript reads as if it belonged to the whole entry.
    ///
    /// DatabaseReference.Molecule keeps it: null means the reference covers the whole entry. It is read only inside a
    /// dbReference (a comment's text-form &lt;molecule&gt;Isoform 3&lt;/molecule&gt; is not one), it survives the
    /// writer and the duplicate-entry merge, and Equals compares it, so two references that differ only by isoform
    /// stay two.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDatabaseReferenceMolecule
    {
        private string _tempDirectory;

        [SetUp]
        public void SetUp()
        {
            _tempDirectory = Path.Combine(Path.GetTempPath(), "TestDatabaseReferenceMolecule_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_tempDirectory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_tempDirectory))
            {
                Directory.Delete(_tempDirectory, true);
            }
        }

        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        private static List<Protein> Load(string path) =>
            ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None, null, false, null, out _);

        private static List<Protein> LoadIsoformScope() => Load(Data("DatabaseTests", "isoformScope.xml"));

        private static DatabaseReference Reference(Protein protein, string type, string id) =>
            protein.DatabaseReferences.Single(r => r.Type == type && r.Id == id);

        private static List<Tuple<string, string>> PathwayName(string name) =>
            new List<Tuple<string, string>> { new Tuple<string, string>("pathway name", name) };

        private string WriteAndReload(List<Protein> proteins, out List<Protein> reloaded)
        {
            string path = Path.Combine(_tempDirectory, "written.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins, path);
            reloaded = Load(path);
            return File.ReadAllText(path);
        }

        [Test]
        public void Molecule_ReadFromDbReference()
        {
            var protein = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");

            Assert.That(Reference(protein, "Reactome", "R-HSA-77289").Molecule, Is.EqualTo("Q9ZZZ1-1"));
            Assert.That(Reference(protein, "Reactome", "R-HSA-556833").Molecule, Is.EqualTo("Q9ZZZ1-2"));
        }

        [Test]
        public void Molecule_NullWhenDbReferenceHasNone()
        {
            var proteins = LoadIsoformScope();
            var first = proteins.Single(p => p.Accession == "Q9ZZZ1");
            var second = proteins.Single(p => p.Accession == "Q9ZZZ2");

            Assert.That(Reference(first, "KEGG", "hsa:9000001").Molecule, Is.Null, "null means the reference covers the whole entry");
            Assert.That(Reference(first, "NCBI Taxonomy", "9606").Molecule, Is.Null);
            Assert.That(Reference(second, "Reactome", "R-HSA-70171").Molecule, Is.Null);
        }

        [Test]
        public void Molecule_DoesNotLeakIntoTheNextReference()
        {
            var protein = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");

            // R-HSA-70171 follows a reference scoped to Q9ZZZ1-1, and the KEGG reference follows one scoped to Q9ZZZ1-2.
            Assert.That(Reference(protein, "Reactome", "R-HSA-70171").Molecule, Is.Null);
            Assert.That(Reference(protein, "KEGG", "hsa:9000001").Molecule, Is.Null);
        }

        [Test]
        public void Molecule_TextMoleculeInACommentIsIgnored()
        {
            var protein = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");

            var scoped = protein.DatabaseReferences.Where(r => r.Molecule != null)
                .Select(r => (r.Type, r.Id, r.Molecule)).OrderBy(r => r.Id, StringComparer.Ordinal).ToList();

            Assert.That(scoped, Is.EqualTo(new[]
            {
                ("Ensembl", "ENST00009000001.1", "Q9ZZZ1-2"),
                ("Reactome", "R-HSA-556833", "Q9ZZZ1-2"),
                ("Reactome", "R-HSA-77289", "Q9ZZZ1-1"),
            }), "the subcellular-location comment's <molecule>Isoform 3</molecule> belongs to no dbReference");
        }

        [Test]
        public void Molecule_ReadOnEnsemblAndCcdsReferencesToo()
        {
            var fixture = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");
            Assert.That(Reference(fixture, "Ensembl", "ENST00009000001.1").Molecule, Is.EqualTo("Q9ZZZ1-2"));

            var gapdh = Load(Data("DatabaseTests", "humanGAPDH.xml")).Single(p => p.Accession == "P04406");
            Assert.That(Reference(gapdh, "CCDS", "CCDS58201.1").Molecule, Is.EqualTo("P04406-2"));
            Assert.That(Reference(gapdh, "CCDS", "CCDS8549.1").Molecule, Is.EqualTo("P04406-1"));
        }

        [Test]
        public void WriteThenRead_KeepsMolecule()
        {
            WriteAndReload(LoadIsoformScope(), out var reloaded);
            var protein = reloaded.Single(p => p.Accession == "Q9ZZZ1");

            Assert.That(Reference(protein, "Reactome", "R-HSA-77289").Molecule, Is.EqualTo("Q9ZZZ1-1"));
            Assert.That(Reference(protein, "Reactome", "R-HSA-556833").Molecule, Is.EqualTo("Q9ZZZ1-2"));
            Assert.That(Reference(protein, "Ensembl", "ENST00009000001.1").Molecule, Is.EqualTo("Q9ZZZ1-2"));
            Assert.That(Reference(protein, "Reactome", "R-HSA-70171").Molecule, Is.Null);
        }

        [Test]
        public void Writer_ReferenceWithoutMolecule_WritesNoMoleculeElement()
        {
            var unscoped = LoadIsoformScope().Where(p => p.Accession == "Q9ZZZ2").ToList();

            string written = WriteAndReload(unscoped, out _);

            Assert.That(written, Does.Not.Contain("<molecule"), "output for references without a molecule is unchanged");
        }

        [Test]
        public void Merge_KeepsReferencesThatDifferOnlyByMolecule()
        {
            var isoformTwo = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"), "P00001-2");
            var isoformThree = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"), "P00001-3");
            var first = new Protein("PEPTIDEK", "P00001", databaseReferences: new List<DatabaseReference> { isoformTwo });
            var second = new Protein("PEPTIDEK", "P00001", databaseReferences: new List<DatabaseReference> { isoformThree });

            var merged = ProteinDbLoader.Merge(new[] { first, second }).Single();

            Assert.That(merged.DatabaseReferences.Select(r => r.Molecule).OrderBy(m => m, StringComparer.Ordinal),
                Is.EqualTo(new[] { "P00001-2", "P00001-3" }), "two isoform claims must not collapse into one");
        }

        [Test]
        public void ThreeArgumentConstructor_LeavesMoleculeNull()
        {
            var reference = new DatabaseReference("Reactome", "R-HSA-70171", PathwayName("Glycolysis"));

            Assert.That(reference.Molecule, Is.Null);
        }

        [Test]
        public void Equals_ComparesMolecule()
        {
            var wholeEntry = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"));
            var isoformTwo = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"), "P00001-2");
            var isoformTwoAgain = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"), "P00001-2");
            var isoformThree = new DatabaseReference("Reactome", "R-HSA-77289", PathwayName("Mitochondrial Fatty Acid Beta-Oxidation"), "P00001-3");

            Assert.That(isoformTwo.Equals(isoformTwoAgain), Is.True);
            Assert.That(isoformTwo.GetHashCode(), Is.EqualTo(isoformTwoAgain.GetHashCode()));
            Assert.That(isoformTwo.Equals(isoformThree), Is.False);
            Assert.That(isoformTwo.Equals(wholeEntry), Is.False);
            Assert.That(wholeEntry.Equals(isoformTwo), Is.False);
        }

        [Test]
        public void Equals_NonDatabaseReference_ReturnsFalse()
        {
            var reference = new DatabaseReference("Reactome", "R-HSA-70171", PathwayName("Glycolysis"));

            Assert.That(reference.Equals("R-HSA-70171"), Is.False);
            Assert.That(reference.Equals(null), Is.False);
        }
    }
}
