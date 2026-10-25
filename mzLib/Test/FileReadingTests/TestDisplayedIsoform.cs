using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// An entry's &lt;sequence&gt; is the sequence of one isoform, the one UniProt displays. Which isoform that is
    /// matters as soon as a reference is scoped with &lt;molecule&gt;: a Reactome pathway listed for P12345-2 applies
    /// to the entry's own sequence only if P12345-2 is the displayed isoform. UniProt says so in
    /// &lt;comment type="alternative products"&gt;, on the &lt;isoform&gt; whose &lt;sequence type="displayed"/&gt;,
    /// and the displayed isoform can carry any number (ProteinAccession's remarks), so "-1" is never assumed.
    ///
    /// UniProtSequenceAttributes.DisplayedIsoformId keeps it beside the sequence's other attributes, so it travels
    /// through every clone, variant and decoy path without a Protein constructor change. Two parser traps are pinned:
    /// an interaction comment's &lt;interactant&gt;&lt;id&gt; sits at the same depth as an isoform's &lt;id&gt;, and an
    /// isoform's own &lt;sequence type="described"/&gt; must not be read as the entry's sequence.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDisplayedIsoform
    {
        private string _tempDirectory;

        [SetUp]
        public void SetUp()
        {
            _tempDirectory = Path.Combine(Path.GetTempPath(), "TestDisplayedIsoform_" + Guid.NewGuid().ToString("N"));
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

        private static List<Protein> Load(string path, DecoyType decoyType = DecoyType.None) =>
            ProteinDbLoader.LoadProteinXML(path, true, decoyType, null, false, null, out _);

        private static List<Protein> LoadIsoformScope(DecoyType decoyType = DecoyType.None) =>
            Load(Data("DatabaseTests", "isoformScope.xml"), decoyType);

        private static UniProtSequenceAttributes AttributesDisplaying(string displayedIsoformId) =>
            new UniProtSequenceAttributes(8, 900, "0000000000000001", new DateTime(2026, 9, 2), 1, null,
                UniProtSequenceAttributes.FragmentType.unspecified, displayedIsoformId);

        private string WriteAndReload(List<Protein> proteins, out List<Protein> reloaded)
        {
            string path = Path.Combine(_tempDirectory, "written.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins, path);
            reloaded = Load(path);
            return File.ReadAllText(path);
        }

        [Test]
        public void DisplayedIsoformId_ReadFromAlternativeProducts()
        {
            var gapdh = Load(Data("DatabaseTests", "humanGAPDH.xml")).Single(p => p.Accession == "P04406");

            Assert.That(gapdh.UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("P04406-1"));
        }

        [Test]
        public void DisplayedIsoformId_IsNotAssumedToBeIsoformOne()
        {
            var protein = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");

            Assert.That(protein.UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("Q9ZZZ1-2"),
                "Q9ZZZ1 lists isoforms -1, -2 and -3 and displays -2");
        }

        [Test]
        public void DisplayedIsoformId_IgnoresInteractantIds()
        {
            var proteins = LoadIsoformScope();

            // Both entries carry an interaction comment whose <interactant><id> is Q9ZZZ1-3; only Q9ZZZ1 lists isoforms.
            Assert.That(proteins.Single(p => p.Accession == "Q9ZZZ1").UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("Q9ZZZ1-2"));
            Assert.That(proteins.Single(p => p.Accession == "Q9ZZZ2").UniProtSequenceAttributes.DisplayedIsoformId, Is.Null);
        }

        [Test]
        public void DisplayedIsoformId_NullWhenEntryListsNoIsoforms()
        {
            string fastaPath = Path.Combine(_tempDirectory, "one.fasta");
            File.WriteAllText(fastaPath, ">sp|Q9ZZZ3|ZZZ3_HUMAN Isoform scope test protein 3 OS=Homo sapiens OX=9606 GN=ZZZ3 PE=1 SV=1\nMSEQWENCER\n");

            var fromFasta = ProteinDbLoader.LoadProteinFasta(fastaPath, true, DecoyType.None, false, out _).Single();
            var inMemory = new Protein("PEPTIDEK", "P00001");

            Assert.That(fromFasta.UniProtSequenceAttributes.DisplayedIsoformId, Is.Null);
            Assert.That(inMemory.UniProtSequenceAttributes.DisplayedIsoformId, Is.Null);
        }

        [Test]
        public void IsoformSequenceElement_DoesNotReplaceEntrySequence()
        {
            var protein = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");
            var gapdh = Load(Data("DatabaseTests", "humanGAPDH.xml")).Single(p => p.Accession == "P04406");

            Assert.That(protein.BaseSequence, Is.EqualTo("MPEPTIDEKAAR"));
            Assert.That(protein.UniProtSequenceAttributes.Checksum, Is.EqualTo("00000000000000A1"));
            Assert.That(protein.UniProtSequenceAttributes.SequenceVersion, Is.EqualTo(2));
            Assert.That(gapdh.UniProtSequenceAttributes.Length, Is.EqualTo(335));
        }

        [Test]
        public void AppliedVariant_CarriesParentDisplayedIsoformId()
        {
            var parent = LoadIsoformScope().Single(p => p.Accession == "Q9ZZZ1");

            var variant = new Protein("MPEPTIDEKAAK", parent, new[] { new SequenceVariation(12, "R", "K", "test variant", null) }, null, null, null);

            Assert.That(variant.Accession, Is.Not.EqualTo(parent.Accession));
            Assert.That(variant.UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("Q9ZZZ1-2"),
                "a variant applied to the displayed sequence is a variant of the displayed isoform");
        }

        [Test]
        public void Decoys_CarryNoIsoformScopedReferences()
        {
            var decoys = LoadIsoformScope(DecoyType.Reverse).Where(p => p.IsDecoy).ToList();

            Assert.That(decoys, Is.Not.Empty);
            Assert.That(decoys.SelectMany(d => d.DatabaseReferences).Select(r => r.Type).Distinct(), Is.EqualTo(new[] { "NCBI Taxonomy" }));
            Assert.That(decoys.SelectMany(d => d.DatabaseReferences).Any(r => r.Molecule != null), Is.False);
        }

        [Test]
        public void Merge_KeepsDisplayedIsoformId()
        {
            var first = new Protein("PEPTIDEK", "P00001", uniProtSequenceAttributes: AttributesDisplaying("P00001-2"));
            var second = new Protein("PEPTIDEK", "P00001", uniProtSequenceAttributes: AttributesDisplaying("P00001-2"));

            var merged = ProteinDbLoader.Merge(new[] { first, second }).Single();

            Assert.That(merged.UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("P00001-2"));
        }

        [Test]
        public void WriteThenRead_KeepsDisplayedIsoformId()
        {
            WriteAndReload(LoadIsoformScope(), out var reloaded);

            Assert.That(reloaded.Single(p => p.Accession == "Q9ZZZ1").UniProtSequenceAttributes.DisplayedIsoformId, Is.EqualTo("Q9ZZZ1-2"));
            Assert.That(reloaded.Single(p => p.Accession == "Q9ZZZ2").UniProtSequenceAttributes.DisplayedIsoformId, Is.Null);
        }

        [Test]
        public void Writer_WritesNoCommentWithoutDisplayedIsoform()
        {
            var withoutIsoforms = LoadIsoformScope().Where(p => p.Accession == "Q9ZZZ2").ToList();

            string written = WriteAndReload(withoutIsoforms, out _);

            Assert.That(written, Does.Not.Contain("alternative products"), "output for entries without isoforms is unchanged");
        }
    }
}
