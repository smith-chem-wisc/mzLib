using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using NUnit.Framework;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Compara's gene-tree peptide alignments: which residue of one protein shares a column with a
    /// residue of another protein of the same tree.
    ///
    /// The fixture has two trees. Tree T1 aligns a human, a mouse and a fish protein; the human protein
    /// has an insertion the mouse lacks, so the human residues there face a gap. Tree T2 holds a second
    /// human protein. The fish protein is not a canonical protein of any gene kept, so it is checked but
    /// not held. Residues are mapped through columns, never by position, so a protein that starts with
    /// gaps still maps correctly.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestComparaGeneTreeAlignment
    {
        private const string AlignmentName = "Compara.116.protein_default.aa.fasta.gz";

        private const string TwoTrees =
            ">ENSP00000000001\n" +
            "MKTAYIAK\n" +
            "QR\n" +
            ">ENSMUSP00000000001\n" +
            "MKT--IAKQR\n" +
            ">ENSDARP00000000001\n" +
            "--TAYIAKQ-\n" +
            "\n//\n\n" +
            ">ENSP00000000002\n" +
            "PEPTIDE\n" +
            "\n//\n";

        private string _dir;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "GeneTreeAlignment_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private string WriteGz(string name, string text)
        {
            string path = Path.Combine(_dir, name);
            using var file = File.Create(path);
            using var gz = new GZipStream(file, CompressionLevel.Optimal);
            using var writer = new StreamWriter(gz);
            writer.Write(text);
            return path;
        }

        private ComparaGeneTreeContent Trees(string extra = "") => ComparaGeneTreeContent.Load(WriteGz(
            "vertebrates.GeneTree_content.default.e116.txt.gz",
            "ENSGT1\tENSP00000000001\tENSG00000000001\tY\n" +
            "ENSGT1\tENSMUSP00000000001\tENSMUSG00000000001\tY\n" +
            "ENSGT1\tENSMUSP00000000009\tENSMUSG00000000001\tN\n" +
            "ENSGT2\tENSP00000000002\tENSG00000000002\tY\n" + extra));

        private ComparaGeneTreeAlignment Load(string text = TwoTrees, ComparaGeneTreeContent trees = null) =>
            ComparaGeneTreeAlignment.Load(WriteGz(AlignmentName, text), trees ?? Trees());

        [Test]
        public void Load_KeepsTheCanonicalProteinsOfTheTrees_AndRecordsTheSource()
        {
            var trees = Trees();
            var a = Load(trees: trees);
            Assert.That(a.Count, Is.EqualTo(3));
            Assert.That(a.Proteins.Select(p => p.ProteinId),
                Is.EqualTo(new[] { "ENSMUSP00000000001", "ENSP00000000001", "ENSP00000000002" }));
            Assert.That(a.AlignmentCount, Is.EqualTo(2));
            Assert.That(a.Release, Is.EqualTo("116"));
            Assert.That(a.Collection, Is.EqualTo("protein_default"));
            Assert.That(a.SourceFileName, Is.EqualTo(AlignmentName));
            Assert.That(a.SourceSha256, Has.Length.EqualTo(64));
            Assert.That(a.RestrictedToGeneTreeContent, Is.EqualTo(trees.SourceSha256));
            Assert.That(a.TryGetProtein("ENSDARP00000000001", out _), Is.False);
        }

        [Test]
        public void Load_JoinsWrappedLines_AndRemovesGapsFromTheSequence()
        {
            var a = Load();
            Assert.That(a.TryGetProtein("ENSP00000000001", out var human), Is.True);
            Assert.That(human.Sequence, Is.EqualTo("MKTAYIAKQR"));
            Assert.That(human.TreeId, Is.EqualTo("ENSGT1"));
            Assert.That(human.AlignmentIndex, Is.EqualTo(0));
            Assert.That(a.TryGetProtein("ENSMUSP00000000001", out var mouse), Is.True);
            Assert.That(mouse.Sequence, Is.EqualTo("MKTIAKQR"));
            Assert.That(a.TryGetProtein("ENSP00000000002", out var second), Is.True);
            Assert.That(second.AlignmentIndex, Is.EqualTo(1));
            Assert.That(second.TreeId, Is.EqualTo("ENSGT2"));
            Assert.That(a.TryGetProtein(null, out _), Is.False);
        }

        [Test]
        public void MapResidue_FollowsTheColumn_AcrossAGapInTheSource()
        {
            var a = Load();
            // Human I is residue 6, column 6; mouse has I at residue 4, column 6.
            var m = a.MapResidue("ENSP00000000001", 6, "ENSMUSP00000000001");
            Assert.That(m, Is.EqualTo(new ComparaColumnMapping(ComparaColumnOutcome.Aligned, 6, 4, 'I')));
            // And back.
            var back = a.MapResidue("ENSMUSP00000000001", 4, "ENSP00000000001");
            Assert.That(back, Is.EqualTo(new ComparaColumnMapping(ComparaColumnOutcome.Aligned, 6, 6, 'I')));
        }

        [Test]
        public void MapResidue_AGapInTheTarget_IsAnOutcome_NotTheNearestResidue()
        {
            var a = Load();
            // Human A (residue 4) and Y (residue 5) sit in columns 4 and 5, where mouse has gaps.
            Assert.That(a.MapResidue("ENSP00000000001", 4, "ENSMUSP00000000001"),
                Is.EqualTo(new ComparaColumnMapping(ComparaColumnOutcome.GapInTarget, 4, null, null)));
            Assert.That(a.MapResidue("ENSP00000000001", 5, "ENSMUSP00000000001").Outcome,
                Is.EqualTo(ComparaColumnOutcome.GapInTarget));
        }

        [Test]
        public void MapResidue_SaysWhyThereIsNoColumn()
        {
            var a = Load();
            Assert.That(a.MapResidue("ENSP00000000001", 1, "ENSP00000000002").Outcome,
                Is.EqualTo(ComparaColumnOutcome.DifferentAlignments));
            Assert.That(a.MapResidue("ENSDARP00000000001", 1, "ENSP00000000001").Outcome,
                Is.EqualTo(ComparaColumnOutcome.SourceNotInAlignment));
            Assert.That(a.MapResidue("ENSP00000000001", 1, "ENSDARP00000000001").Outcome,
                Is.EqualTo(ComparaColumnOutcome.TargetNotInAlignment));
            Assert.Throws<ArgumentOutOfRangeException>(() => a.MapResidue("ENSP00000000001", 11, "ENSMUSP00000000001"));
            Assert.Throws<ArgumentOutOfRangeException>(() => a.MapResidue("ENSP00000000001", 0, "ENSMUSP00000000001"));
            Assert.Throws<ArgumentNullException>(() => a.MapResidue(null, 1, "ENSMUSP00000000001"));
        }

        [Test]
        public void Load_AcceptsAFileThatDoesNotEndWithTheSeparator()
        {
            var a = Load(">ENSP00000000002\nPEPTIDE\n");
            Assert.That(a.AlignmentCount, Is.EqualTo(1));
            Assert.That(a.TryGetProtein("ENSP00000000002", out var p), Is.True);
            Assert.That(p.Sequence, Is.EqualTo("PEPTIDE"));
        }

        private static IEnumerable<TestCaseData> Malformed()
        {
            yield return new TestCaseData("PEPTIDE\n>ENSP00000000002\nPEPTIDE\n", "line 1: a sequence line before any header")
                .SetName("Load_RefusesASequenceBeforeAnyHeader");
            yield return new TestCaseData(">\nPEPTIDE\n", "line 1: a header with no protein id")
                .SetName("Load_RefusesAnEmptyHeader");
            yield return new TestCaseData(">ENSP00000000002\nPEP1IDE\n", "line 2: character '1' in ENSP00000000002")
                .SetName("Load_RefusesACharacterThatIsNotAResidue");
            yield return new TestCaseData(">ENSDARP00000000001\nPEPT\n>ENSDARP00000000002\nPEP\n//\n", "ENSDARP00000000002 is 3 columns wide, the alignment 4")
                .SetName("Load_RefusesRaggedRows_EvenOnRowsThatAreNotKept");
            yield return new TestCaseData(">ENSP00000000002\nPEPTIDE\n//\n//\n", "line 4: an alignment with no rows")
                .SetName("Load_RefusesAnEmptyAlignment");
            yield return new TestCaseData(">ENSP00000000002\nPEPTIDE\n//\n>ENSP00000000002\nPEPTIDE\n//\n", "ENSP00000000002 appears twice")
                .SetName("Load_RefusesAKeptProteinThatAppearsTwice");
            yield return new TestCaseData(">ENSP00000000001\nPEPTIDE\n>ENSP00000000002\nPEPTIDE\n//\n", "alignment 0 holds ENSP00000000001 of ENSGT1 and ENSP00000000002 of ENSGT2")
                .SetName("Load_RefusesAnAlignmentHoldingTwoTrees");
        }

        [TestCaseSource(nameof(Malformed))]
        public void Load_RefusesRatherThanGuesses(string text, string message)
        {
            var ex = Assert.Throws<InvalidDataException>(() => Load(text));
            Assert.That(ex!.Message, Does.Contain(message));
            Assert.That(ex.Message, Does.StartWith(AlignmentName));
        }

        [Test]
        public void Load_MissingFileThrows()
        {
            Assert.Throws<FileNotFoundException>(() =>
                ComparaGeneTreeAlignment.Load(Path.Combine(_dir, "absent.aa.fasta.gz"), Trees()));
        }

        [Test]
        public void Load_RequiresTheTreesToKeep()
        {
            Assert.Throws<ArgumentNullException>(() => ComparaGeneTreeAlignment.Load(WriteGz(AlignmentName, TwoTrees), null));
        }

        /// <summary>
        /// Release 116's real file (908 MB). Set MZLIB_COMPARA_ALIGNMENT to its path, MZLIB_COMPARA_TREES to
        /// vertebrates.GeneTree_content.default.e116.txt.gz, and MZLIB_GENE_SETS to the gene tables to keep
        /// (';'-separated). Keeping every gene of the collection would hold millions of proteins.
        /// </summary>
        [Test]
        [Explicit("Requires Ensembl 116's Compara.116.protein_default.aa.fasta.gz; set MZLIB_COMPARA_ALIGNMENT, MZLIB_COMPARA_TREES and MZLIB_GENE_SETS.")]
        public void Load_RealRelease116_ReadsEveryAlignment()
        {
            string alignment = Environment.GetEnvironmentVariable("MZLIB_COMPARA_ALIGNMENT");
            string trees = Environment.GetEnvironmentVariable("MZLIB_COMPARA_TREES");
            var sets = Environment.GetEnvironmentVariable("MZLIB_GENE_SETS").Split(';').Select(EnsemblGeneSetReader.Load).ToList();
            var content = ComparaGeneTreeContent.Load(trees, sets);
            var a = ComparaGeneTreeAlignment.Load(alignment, content);
            Assert.That(a.Count, Is.EqualTo(content.Count));
            Assert.That(a.Release, Is.EqualTo("116"));
            TestContext.WriteLine($"{a.AlignmentCount} alignments, {a.Count} proteins kept");
        }
    }
}
