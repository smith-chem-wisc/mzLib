using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using UsefulProteomicsDatabases.Ensembl;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Ensembl Compara's homologies among a chosen set of species, and the two dumps they are read from.
    ///
    /// The fixture is a three-species world built to hold the two traps the builder exists for. Its
    /// human-mouse relationships sit only in the MOUSE dump, as every one does in release 116, so a
    /// builder that read only the dump named after each species would lose them. And H2~M2 and M2~R2 are
    /// orthologs while H2~R2 is not called, so a builder that chained calls would invent a triple.
    ///
    /// Checked on real data: human, mouse and rat release 116 give 23,764 human-mouse rows, human->mouse
    /// statuses 17,826 / 531 / 1,333 / 853 over protein-coding and IG/TR genes, and 15,508 human genes
    /// in an all-one-to-one triple, each equal to an independent implementation.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestOrthologySnapshot
    {
        private const string Human = "homo_sapiens", Mouse = "mus_musculus", Rat = "rattus_norvegicus";
        private const string H1 = "ENSG00000000001", H2 = "ENSG00000000002", H3 = "ENSG00000000003", H4 = "ENSG00000000004";
        private const string M1 = "ENSMUSG00000000001", M2 = "ENSMUSG00000000002", M3 = "ENSMUSG00000000003";
        private const string R1 = "ENSRNOG00000000001", R2 = "ENSRNOG00000000002";

        private static readonly string Header = string.Join("\t", ComparaHomologyDump.Columns) + "\n";

        private static string Row(string id, string type, string sp, string gene, double identity,
            string otherSp, string other, double otherIdentity, string goc = "100", string high = "1") =>
            $"{gene}\t{gene}P\t{sp}\t{identity}\t{type}\t{other}\t{other}P\t{otherSp}\t{otherIdentity}\tNULL\tNULL\t{goc}\tNULL\t{high}\t{id}\n";

        private string _dir;
        private Dictionary<string, EnsemblGeneSet> _sets;

        [SetUp]
        public void SetUp()
        {
            _dir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "Orthology_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_dir);
            _sets = new Dictionary<string, EnsemblGeneSet>
            {
                [Human] = GeneSet("Homo_sapiens.GRCh38.116.gtf", H1, H2, H3, H4),
                [Mouse] = GeneSet("Mus_musculus.GRCm39.116.gtf", M1, M2, M3),
                [Rat] = GeneSet("Rattus_norvegicus.GRCr8.116.gtf", R1, R2),
            };
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_dir)) Directory.Delete(_dir, true);
        }

        private EnsemblGeneSet GeneSet(string name, params string[] genes)
        {
            string path = Path.Combine(_dir, name);
            File.WriteAllText(path, string.Concat(genes.Select(g =>
                $"1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"{g}\"; gene_version \"1\"; gene_biotype \"protein_coding\";\n")));
            return EnsemblGeneSet.LoadGtf(path);
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

        private static readonly string[] All = { Human, Mouse, Rat };

        private ComparaHomologyDump HumanDump(IEnumerable<string> species = null) => ComparaHomologyDump.Load(
            WriteGz("human.homologies.tsv.gz", Header +
                Row("10", "ortholog_one2one", Human, H1, 80, Rat, R1, 81)),
            species ?? All);

        private ComparaHomologyDump MouseDump(IEnumerable<string> species = null) => ComparaHomologyDump.Load(
            WriteGz("Compara.116.protein_default.homologies.tsv.gz", Header +
                Row("11", "ortholog_one2one", Mouse, M1, 90, Human, H1, 91) +
                Row("12", "ortholog_one2one", Mouse, M1, 95, Rat, R1, 96) +
                Row("13", "ortholog_one2many", Mouse, M2, 70, Human, H2, 71, goc: "NULL", high: "NULL") +
                Row("14", "ortholog_one2many", Mouse, M3, 60, Human, H2, 61, high: "0") +
                Row("15", "within_species_paralog", Mouse, M2, 50, Mouse, M3, 51)),
            species ?? All);

        private ComparaHomologyDump RatDump(string extra = "") => ComparaHomologyDump.Load(
            WriteGz("rat.homologies.tsv.gz", Header +
                Row("16", "ortholog_one2one", Rat, R2, 88, Mouse, M2, 89) + extra),
            All);

        private ComparaGeneTreeContent Trees(IEnumerable<EnsemblGeneSet> restrictTo = null) => ComparaGeneTreeContent.Load(
            WriteGz("vertebrates.GeneTree_content.default.e116.txt.gz",
                $"ENSGT1\t{H1}P\t{H1}\tY\nENSGT1\t{M1}P\t{M1}\tY\nENSGT1\t{R1}P\t{R1}\tY\nENSGT1\t{R1}Q\t{R1}\tN\n" +
                $"ENSGT2\t{H2}P\t{H2}\tY\nENSGT2\t{M2}P\t{M2}\tY\nENSGT2\t{M3}P\t{M3}\tY\nENSGT2\t{R2}P\t{R2}\tY\n" +
                $"ENSGT3\t{H3}P\t{H3}\tY\nENSGT9\tENSDARP1\tENSDARG00000000001\tY\n"),
            restrictTo);

        private OrthologySnapshot Standard() =>
            OrthologySnapshot.Build("116", _sets, Trees(_sets.Values), new[] { HumanDump(), MouseDump(), RatDump() });

        // ---- the homology dump ----

        [Test]
        public void Dump_KeepsRowsOfTheRequestedSpecies_ValuesVerbatim()
        {
            var dump = MouseDump(new[] { Human, Mouse });

            Assert.That(dump.Genome, Is.EqualTo(Mouse));
            Assert.That(dump.Rows.Select(r => r.HomologyId), Is.EqualTo(new[] { "11", "13", "14", "15" }),
                "the mouse-rat row is checked but not kept");
            Assert.That(dump.RowsRead, Is.EqualTo(5));
            Assert.That(dump.Species, Is.EqualTo(new[] { Human, Mouse }));

            var r13 = dump.Rows.Single(r => r.HomologyId == "13");
            Assert.That((r13.GeneA, r13.IdentityA, r13.GeneB, r13.IdentityB), Is.EqualTo((M2, 70.0, H2, 71.0)));
            Assert.That((r13.GocScore, r13.IsHighConfidence, r13.Dn), Is.EqualTo(((int?)null, (bool?)null, (double?)null)),
                "NULL is the source declining to say; it is not 0 or false");
            Assert.That(dump.Rows.Single(r => r.HomologyId == "14").IsHighConfidence, Is.False);
            Assert.That(dump.Rows.Single(r => r.HomologyId == "15").Class, Is.EqualTo(HomologyClass.Paralog));
        }

        [Test]
        public void Dump_ReleaseFromEnsemblsFileName_AndSha256OfTheCompressedBytes()
        {
            var dump = MouseDump();
            string path = Path.Combine(_dir, "Compara.116.protein_default.homologies.tsv.gz");

            Assert.That(dump.Release, Is.EqualTo("116"));
            Assert.That(dump.SourceSha256, Is.EqualTo(Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant()));
            Assert.That(HumanDump().Release, Is.Null, "a name without the release gives null, not a guess");
        }

        private static readonly object[] MalformedDumps =
        {
            new object[] { "gene_stable_id\tspecies\n", "unexpected columns" },
            new object[] { Header + "ENSG1\tP\thomo_sapiens\n", "line 2: 3 cells, expected 15" },
            new object[] { Header + Row("1", "ortholog_one2few", Human, H1, 1, Mouse, M1, 1), "unknown homology_type 'ortholog_one2few'" },
            new object[] { Header + Row("1", "other_paralog", Human, H1, 1, Mouse, M1, 1), "is a other_paralog between" },
            new object[] { Header + Row("1", "ortholog_one2one", Human, H1, 1, Human, H2, 1), "is a ortholog_one2one within" },
            new object[] { Header + Row("1", "ortholog_one2one", Human, H1, 1, Mouse, M1, 1) + Row("2", "ortholog_one2one", Mouse, M1, 1, Rat, R1, 1),
                "line 3: species mus_musculus in a dump of homo_sapiens" },
            new object[] { Header + Row("1", "ortholog_one2one", Human, H1, 1, Mouse, M1, 1).Replace("\t1\tortholog", "\tx\tortholog"), "identity 'x' is not a number" },
            new object[] { Header + Row("1", "ortholog_one2one", Human, H1, 1, Mouse, M1, 1, goc: "12.5"), "goc_score '12.5' is not a whole number" },
            new object[] { Header + Row("1", "ortholog_one2one", Human, H1, 1, Mouse, M1, 1, high: "yes"), "is_high_confidence 'yes'" },
            new object[] { Header + Row("1", "ortholog_one2few", "danio_rerio", "D1", 1, "danio_rerio", "D2", 1), "unknown homology_type" },
        };

        [TestCaseSource(nameof(MalformedDumps))]
        public void Dump_RefusesRatherThanGuesses(string text, string message)
        {
            var ex = Assert.Throws<InvalidDataException>(() => ComparaHomologyDump.Load(WriteGz("bad.tsv.gz", text), All));
            Assert.That(ex.Message, Does.Contain(message));
        }

        [Test]
        public void Dump_MissingFileOrNoSpeciesThrows()
        {
            Assert.Throws<FileNotFoundException>(() => ComparaHomologyDump.Load(Path.Combine(_dir, "absent.tsv.gz"), All));
            Assert.Throws<ArgumentException>(() => ComparaHomologyDump.Load(WriteGz("x.tsv.gz", Header), Array.Empty<string>()));
        }

        [Test]
        public void Dump_ClassNames()
        {
            Assert.That(Enum.GetValues<HomologyClass>().Select(ComparaHomologyDump.ClassName),
                Is.EqualTo(new[] { "ortholog", "paralog", "homoeolog" }));
            Assert.That(ComparaHomologyDump.ClassOf("gene_split"), Is.EqualTo(HomologyClass.Paralog));
            Assert.Throws<ArgumentOutOfRangeException>(() => ComparaHomologyDump.ClassOf("ortholog"));
            Assert.That(ComparaHomologyDump.KnownHomologyTypes, Has.Count.EqualTo(9));
        }

        // ---- the gene-tree dump ----

        [Test]
        public void Trees_OneTreeAndOneCanonicalProteinPerGene()
        {
            var trees = Trees();

            Assert.That(trees.TryGetMember(R1, out var r1), Is.True);
            Assert.That(r1, Is.EqualTo(new ComparaGeneTreeMember(R1, "ENSGT1", R1 + "P")), "the N row is not canonical");
            Assert.That(trees.Count, Is.EqualTo(9));
            Assert.That(trees.RestrictedToGeneSets, Is.Null);
            Assert.That((trees.Release, trees.Collection), Is.EqualTo(("116", "vertebrates")));
            Assert.That(trees.Members.First().GeneId, Is.EqualTo("ENSDARG00000000001"), "ordinal order");
        }

        [Test]
        public void Trees_RestrictedToGeneSets_RecordsWhichOnes()
        {
            var trees = Trees(new[] { _sets[Human] });

            Assert.That(trees.Count, Is.EqualTo(3));
            Assert.That(trees.TryGetMember(M1, out _), Is.False);
            Assert.That(trees.RestrictedToGeneSets, Is.EqualTo(new[] { _sets[Human].SourceSha256 }));
        }

        private static readonly object[] MalformedTrees =
        {
            new object[] { "ENSGT1\tP1\tG1\n", "line 1: 3 cells, expected 4" },
            new object[] { "ENSGT1\tP1\tG1\tyes\n", "canonical flag 'yes'" },
            new object[] { "ENSGT1\tP1\tG1\tY\nENSGT2\tP2\tG1\tN\n", "line 2: G1 is in two gene trees, ENSGT1 and ENSGT2" },
            new object[] { "ENSGT1\tP1\tG1\tY\nENSGT1\tP2\tG1\tY\n", "G1 has 2 canonical proteins" },
            new object[] { "ENSGT1\tP1\tG1\tN\n", "G1 has 0 canonical proteins" },
        };

        [TestCaseSource(nameof(MalformedTrees))]
        public void Trees_RefuseRatherThanGuess(string text, string message)
        {
            var ex = Assert.Throws<InvalidDataException>(() => ComparaGeneTreeContent.Load(WriteGz("t.txt.gz", text)));
            Assert.That(ex.Message, Does.Contain(message));
        }

        // ---- the snapshot ----

        [Test]
        public void Build_ReadsEveryDump_BecauseAPairsRowsSitInEitherOne()
        {
            var snap = Standard();

            Assert.That(snap.Pair(Human, Mouse).Select(r => r.HomologyId), Is.EqualTo(new[] { "11", "13", "14" }),
                "every human-mouse row came from the mouse dump");
            Assert.That(snap.Pair(Human, Rat).Select(r => r.HomologyId), Is.EqualTo(new[] { "10" }));
            Assert.That(snap.Pair(Mouse, Rat).Select(r => r.HomologyId), Is.EqualTo(new[] { "12", "16" }));
            Assert.That(snap.Pair(Mouse, Mouse).Single().HomologyType, Is.EqualTo("within_species_paralog"));
            Assert.That(snap.Pair(Human, Human), Is.Empty);
            Assert.That(snap.Pairs, Is.EqualTo(new[]
            {
                (Human, Human), (Human, Mouse), (Human, Rat), (Mouse, Mouse), (Mouse, Rat), (Rat, Rat)
            }));
            Assert.That(snap.IdenticalDuplicatesDropped, Is.Zero);
        }

        [Test]
        public void Build_OrientsSideAToTheSpeciesThatSortsFirst_IdentityMovesWithItsGene()
        {
            var snap = Standard();

            var r11 = snap.Pair(Mouse, Human).Single(r => r.HomologyId == "11");
            Assert.That((r11.SpeciesA, r11.GeneA, r11.IdentityA, r11.GeneB, r11.IdentityB),
                Is.EqualTo((Human, H1, 91.0, M1, 90.0)), "written mouse-first; stored human-first");

            var fromMouse = snap.Orthologs(Mouse, Human).Single(r => r.HomologyId == "11");
            Assert.That((fromMouse.GeneA, fromMouse.IdentityA), Is.EqualTo((M1, 90.0)));
            Assert.That(fromMouse.SourceFile, Is.EqualTo("Compara.116.protein_default.homologies.tsv.gz"));
        }

        [Test]
        public void PairStatus_SaysHowFarTheInferenceGot()
        {
            var snap = Standard();

            Assert.That(snap.PairStatus(Human, Mouse).Select(x => (x.Gene.GeneId, x.Status)), Is.EqualTo(new[]
            {
                (H1, OrthologyStatus.HasOrtholog),
                (H2, OrthologyStatus.HasOrtholog),
                (H3, OrthologyStatus.TreeLacksTargetSpecies),
                (H4, OrthologyStatus.NotInAnyTree),
            }));
            Assert.That(snap.StatusOf(H2, Rat), Is.EqualTo(OrthologyStatus.NoEdgeInSharedTree),
                "ENSGT2 holds R2, and no H2-R2 ortholog was called");
            Assert.That(snap.StatusOf("ENSG00000099999", Mouse), Is.EqualTo(OrthologyStatus.NotInGeneSet));
            Assert.That(snap.StatusOf(null, Mouse), Is.EqualTo(OrthologyStatus.NotInGeneSet));
            Assert.That(snap.SpeciesOf(M3), Is.EqualTo(Mouse));
        }

        [Test]
        public void SpeciesSet_RequiresEveryPair_NeverChains()
        {
            var snap = Standard();

            var triples = snap.SpeciesSet(Human, Mouse, Rat).ToList();
            Assert.That(triples, Has.Count.EqualTo(1), "H2~M2 and M2~R2 without H2~R2 is not a triple");
            Assert.That(triples[0].GeneIds, Is.EqualTo(new[] { H1, M1, R1 }));
            Assert.That(triples[0].AllOneToOne, Is.True);

            var pairs = snap.SpeciesSet(Human, Mouse).ToList();
            Assert.That(pairs.Select(t => (t.GeneIds[0], t.GeneIds[1], t.AllOneToOne)), Is.EqualTo(new[]
            {
                (H1, M1, true), (H2, M2, false), (H2, M3, false)
            }), "one-to-many is kept as two rows, never collapsed or picked");
            Assert.That(snap.SpeciesSet(Rat, Human, Mouse).Single().GeneIds, Is.EqualTo(new[] { R1, H1, M1 }),
                "genes come in the order the species were given");
        }

        [Test]
        public void Views_RefuseSpeciesTheyCannotAnswerFor()
        {
            var snap = Standard();

            Assert.Throws<ArgumentException>(() => snap.PairStatus(Human, Human));
            Assert.Throws<ArgumentException>(() => snap.Orthologs(Human, "danio_rerio"));
            Assert.Throws<ArgumentException>(() => snap.SpeciesSet(Human));
            Assert.Throws<ArgumentException>(() => snap.SpeciesSet(Human, Human));
            Assert.Throws<ArgumentException>(() => snap.StatusOf(H1, Human));
        }

        [Test]
        public void Build_WithoutEverySpeciesDump_IsRefused()
        {
            var ex = Assert.Throws<ArgumentException>(() =>
                OrthologySnapshot.Build("116", _sets, Trees(), new[] { HumanDump(), RatDump() }));
            Assert.That(ex.Message, Does.Contain("No homology dump for mus_musculus"));
        }

        [Test]
        public void Build_IdenticalRowInBothDumps_IsCountedOnce_ADifferentOneIsRefused()
        {
            var same = OrthologySnapshot.Build("116", _sets, Trees(),
                new[] { HumanDump(), MouseDump(), RatDump(Row("12", "ortholog_one2one", Rat, R1, 96, Mouse, M1, 95)) });
            Assert.That(same.IdenticalDuplicatesDropped, Is.EqualTo(1), "the same relationship, written from the other side");
            Assert.That(same.Pair(Mouse, Rat).Count(r => r.HomologyId == "12"), Is.EqualTo(1));

            var ex = Assert.Throws<InvalidDataException>(() => OrthologySnapshot.Build("116", _sets, Trees(),
                new[] { HumanDump(), MouseDump(), RatDump(Row("12", "ortholog_one2many", Rat, R1, 96, Mouse, M1, 95)) }));
            Assert.That(ex.Message, Does.Contain("Homology 12 has different rows"));
        }

        [Test]
        public void Build_AGeneOutsideItsSpeciesGeneSet_IsRefused()
        {
            var ex = Assert.Throws<InvalidDataException>(() => OrthologySnapshot.Build("116", _sets, Trees(),
                new[] { HumanDump(), MouseDump(), RatDump(Row("17", "ortholog_one2one", Rat, "ENSRNOG00000000099", 1, Mouse, M1, 1)) }));
            Assert.That(ex.Message, Does.Contain("names ENSRNOG00000000099, which is not in the rattus_norvegicus gene set"));
        }

        [Test]
        public void Build_InputsThatDoNotFitTogether_AreRefused()
        {
            var trees = Trees();
            Assert.That(Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("117", _sets, trees,
                new[] { HumanDump(), MouseDump(), RatDump() })).Message, Does.Contain("is from release 116, not 117"));

            Assert.That(Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("116", _sets, trees,
                new[] { HumanDump(new[] { Human, Rat }), MouseDump(), RatDump() })).Message, Does.Contain("was loaded without mus_musculus"));

            Assert.That(Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("116", _sets, Trees(new[] { _sets[Human] }),
                new[] { HumanDump(), MouseDump(), RatDump() })).Message, Does.Contain("restricted to other gene sets"));

            Assert.That(Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("116", _sets, trees,
                new[] { HumanDump(), MouseDump(), MouseDump(), RatDump() })).Message, Does.Contain("Two dumps of mus_musculus"));

            var twoSpecies = new Dictionary<string, EnsemblGeneSet> { [Human] = _sets[Human], [Mouse] = _sets[Mouse] };
            Assert.That(Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("116", twoSpecies, trees,
                new[] { HumanDump(), MouseDump(), RatDump() })).Message, Does.Contain("not one of the species"));

            Assert.Throws<ArgumentException>(() => OrthologySnapshot.Build("116", new Dictionary<string, EnsemblGeneSet>(), trees, Array.Empty<ComparaHomologyDump>()));
        }

        [Test]
        public void Build_ASubsetOfTheLoadedSpecies_DropsTheRest()
        {
            var twoSpecies = new Dictionary<string, EnsemblGeneSet> { [Human] = _sets[Human], [Mouse] = _sets[Mouse] };
            var snap = OrthologySnapshot.Build("116", twoSpecies, Trees(), new[] { HumanDump(), MouseDump() });

            Assert.That(snap.Species, Is.EqualTo(new[] { Human, Mouse }));
            Assert.That(snap.Pair(Human, Mouse), Has.Count.EqualTo(3));
            Assert.That(snap.StatusOf(H1, Mouse), Is.EqualTo(OrthologyStatus.HasOrtholog));
            Assert.That(snap.Dumps.Select(d => d.Genome), Is.EqualTo(new[] { Human, Mouse }));
        }
    }
}
