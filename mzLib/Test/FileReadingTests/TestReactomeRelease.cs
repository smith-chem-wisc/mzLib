using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using NUnit.Framework;
using UsefulProteomicsDatabases.Pathways;

namespace Test.FileReadingTests
{
    /// <summary>
    /// A Reactome release, read from Reactome's two flat files: ReactomePathways.txt (id, name, species; every
    /// species) and ReactomePathwaysRelation.txt (parent, child). UniProt cross-references only the lowest-level
    /// pathways a protein takes part in, so rolling a membership up to "Metabolism" or "Disease" needs this hierarchy,
    /// pinned by sha256 like go.obo.
    ///
    /// What these pin:
    ///  - the release number comes from the caller, because neither file carries one;
    ///  - the hierarchy is a DAG (388 pathways in release 96 have two or more parents), so ancestors are a closure,
    ///    counted once each, and a cycle or an edge to a pathway the file lacks is refused rather than walked;
    ///  - the human equivalent of an inferred pathway is R-HSA with the same number AND the same name, and nothing
    ///    else: a chicken pathway Reactome curated separately has none.
    /// The fixture is release 96 cut to the ancestor closure of six seeds (DataFiles/Pathways/PROVENANCE.md).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestReactomeRelease
    {
        private string _tempDirectory;

        [SetUp]
        public void SetUp()
        {
            _tempDirectory = Path.Combine(Path.GetTempPath(), "TestReactomeRelease_" + Guid.NewGuid().ToString("N"));
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

        internal static string PathwaysFixture =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Pathways", "ReactomePathways-96-trimmed.txt");

        internal static string RelationFixture =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Pathways", "ReactomePathwaysRelation-96-trimmed.txt");

        internal static ReactomeRelease LoadFixture() => ReactomeRelease.Load(PathwaysFixture, RelationFixture, "96");

        private (string Pathways, string Relation) Write(string pathways, string relation)
        {
            string pathwaysPath = Path.Combine(_tempDirectory, "pathways.txt");
            string relationPath = Path.Combine(_tempDirectory, "relation.txt");
            File.WriteAllText(pathwaysPath, pathways);
            File.WriteAllText(relationPath, relation);
            return (pathwaysPath, relationPath);
        }

        private static string Sha256Of(string path) =>
            Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

        [Test]
        public void Load_ReadsTheTrimmedRelease96Slice()
        {
            var release = LoadFixture();

            Assert.That(release.Release, Is.EqualTo("96"));
            Assert.That(release.Pathways.Count, Is.EqualTo(25));
            Assert.That(release.Hierarchy.Count, Is.EqualTo(19));
            var integrin = release.Pathways.Single(p => p.PathwayId == "R-HSA-354192");
            Assert.That(integrin.Name, Is.EqualTo("Integrin signaling"));
            Assert.That(integrin.SpeciesName, Is.EqualTo("Homo sapiens"));
            Assert.That(release.Pathways.Select(p => p.PathwayId), Is.Ordered.Using((IComparer<string>)StringComparer.Ordinal));
        }

        [Test]
        public void Load_HashesBothFilesAsLowercaseHex()
        {
            var release = LoadFixture();

            Assert.That(release.PathwaysSha256, Is.EqualTo(Sha256Of(PathwaysFixture)));
            Assert.That(release.RelationSha256, Is.EqualTo(Sha256Of(RelationFixture)));
            Assert.That(release.PathwaysSha256, Does.Match("^[0-9a-f]{64}$"));
        }

        [TestCase(null)]
        [TestCase("")]
        public void Load_RequiresTheRelease(string release)
        {
            Assert.Throws<ArgumentException>(() => ReactomeRelease.Load(PathwaysFixture, RelationFixture, release));
        }

        [Test]
        public void Load_MissingFile_Throws()
        {
            Assert.Throws<FileNotFoundException>(() => ReactomeRelease.Load(Path.Combine(_tempDirectory, "absent.txt"), RelationFixture, "96"));
            Assert.Throws<FileNotFoundException>(() => ReactomeRelease.Load(PathwaysFixture, Path.Combine(_tempDirectory, "absent.txt"), "96"));
        }

        [TestCase("R-HSA-1\tA\tHomo sapiens\nR-HSA-2\tB\tHomo sapiens\n", "R-HSA-1\tR-HSA-2\nR-HSA-2\tR-HSA-1\n", "cycle")]
        [TestCase("R-HSA-1\tA\tHomo sapiens\n", "R-HSA-1\tR-HSA-1\n", "cycle")]
        [TestCase("R-HSA-1\tA\tHomo sapiens\n", "R-HSA-1\tR-HSA-9\n", "R-HSA-9")]
        [TestCase("R-HSA-1\tA\tHomo sapiens\n", "R-HSA-9\tR-HSA-1\n", "R-HSA-9")]
        [TestCase("R-HSA-1\tA\tHomo sapiens\nR-HSA-1\tB\tHomo sapiens\n", "", "twice")]
        [TestCase("R-HSA-1\tA\n", "", "3 cells")]
        [TestCase("\tA\tHomo sapiens\n", "", "empty")]
        [TestCase("R-HSA-1\tA\tHomo sapiens\nR-HSA-2\tB\tHomo sapiens\n", "R-HSA-1\tR-HSA-2\tx\n", "2 cells")]
        public void Load_RefusesAMalformedRelease(string pathways, string relation, string expected)
        {
            var (pathwaysPath, relationPath) = Write(pathways, relation);

            var error = Assert.Throws<InvalidDataException>(() => ReactomeRelease.Load(pathwaysPath, relationPath, "1"));

            Assert.That(error.Message, Does.Contain(expected));
        }

        [Test]
        public void Ancestors_OfMitochondrialBetaOxidation_AreItsChain()
        {
            var release = LoadFixture();

            Assert.That(release.Ancestors("R-HSA-77289"), Is.EquivalentTo(new[] { "R-HSA-8978868", "R-HSA-556833", "R-HSA-1430728" }));
            Assert.That(release.Ancestors("R-MMU-77289"), Is.EquivalentTo(new[] { "R-MMU-8978868", "R-MMU-556833", "R-MMU-1430728" }));
            Assert.That(release.Ancestors("R-HSA-1430728"), Is.Empty, "a top-level pathway has no ancestors");
        }

        [Test]
        public void Ancestors_OfAPathwayWithTwoParents_ReachBothRoots()
        {
            var release = LoadFixture();

            Assert.That(release.Ancestors("R-HSA-354192"), Is.EquivalentTo(new[] { "R-HSA-162582", "R-HSA-76009", "R-HSA-76002", "R-HSA-109582" }));
        }

        [Test]
        public void Ancestors_OfADiamond_CountEachAncestorOnce()
        {
            var (pathwaysPath, relationPath) = Write(
                "R-HSA-1\tTop\tHomo sapiens\nR-HSA-2\tLeft\tHomo sapiens\nR-HSA-3\tRight\tHomo sapiens\nR-HSA-4\tBottom\tHomo sapiens\n",
                "R-HSA-1\tR-HSA-2\nR-HSA-1\tR-HSA-3\nR-HSA-2\tR-HSA-4\nR-HSA-3\tR-HSA-4\n");

            var release = ReactomeRelease.Load(pathwaysPath, relationPath, "1");

            Assert.That(release.Ancestors("R-HSA-4"), Is.EquivalentTo(new[] { "R-HSA-1", "R-HSA-2", "R-HSA-3" }));
        }

        [Test]
        public void Ancestors_UnknownId_Throws()
        {
            var release = LoadFixture();

            Assert.Throws<ArgumentException>(() => release.Ancestors("R-HSA-0"));
            Assert.That(release.Contains("R-HSA-0"), Is.False);
            Assert.That(release.Contains("R-HSA-77289"), Is.True);
        }

        [Test]
        public void Roots_AreTheTopLevelPathways()
        {
            var release = LoadFixture();

            Assert.That(release.Roots, Is.EqualTo(new[]
            {
                "R-GGA-1670466", "R-HSA-109582", "R-HSA-1430728", "R-HSA-162582", "R-HSA-1643685", "R-MMU-1430728", "R-RNO-1430728",
            }));
        }

        [Test]
        public void HumanEquivalent_HumanIdIsItself()
        {
            var release = LoadFixture();

            Assert.That(release.HumanEquivalentId("R-HSA-5602358"), Is.EqualTo("R-HSA-5602358"));
        }

        [Test]
        public void HumanEquivalent_SameNumberAndSameName_IsTheHumanId()
        {
            var release = LoadFixture();

            Assert.That(release.HumanEquivalentId("R-MMU-77289"), Is.EqualTo("R-HSA-77289"));
            Assert.That(release.HumanEquivalentId("R-RNO-1430728"), Is.EqualTo("R-HSA-1430728"));
            Assert.That(release.Pathways.Single(p => p.PathwayId == "R-RNO-556833").HumanEquivalentId, Is.EqualTo("R-HSA-556833"));
        }

        [Test]
        public void HumanEquivalent_SameNumberDifferentName_IsNull()
        {
            var (pathwaysPath, relationPath) = Write("R-HSA-1\tAlpha\tHomo sapiens\nR-MMU-1\tBeta\tMus musculus\n", "");

            var release = ReactomeRelease.Load(pathwaysPath, relationPath, "1");

            Assert.That(release.HumanEquivalentId("R-MMU-1"), Is.Null, "a number shared by two differently named events is not the same pathway");
        }

        [Test]
        public void HumanEquivalent_NoHumanId_IsNull()
        {
            var release = LoadFixture();

            Assert.That(release.HumanEquivalentId("R-GGA-351470"), Is.Null, "a curated chicken pathway with no R-HSA twin");
            Assert.Throws<ArgumentException>(() => release.HumanEquivalentId("R-HSA-0"));
        }
    }
}
