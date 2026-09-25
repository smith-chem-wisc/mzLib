using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Test.Omics;
using Test.Omics.BioPolymerGroupTests;

namespace Test.FileReadingTests.InternalFileReading
{
    /// <summary>
    /// The fixture is six rows of the AllQuantifiedProteinGroups.tsv MetaMorpheus 1.1.11 wrote for the
    /// public PRIDE dataset PXD036557 (18 label-free files), chosen for the cases a reader gets wrong:
    /// a plain group (P68363), a modification name containing commas (P05141), a two-member group
    /// whose occupancy cells hold two entities (P0C0S5|Q71UI9), a contaminant (P02769), a protein
    /// N-terminal site at position 0 (P63104) and a decoy (DECOY_P62750). The header is the full
    /// 94-column one, in MetaMorpheus 1.1.x's column vocabulary.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestProteinGroupFromTsv
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ExternalFileTypes\MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv");

        private string _outputDirectory = "";

        [OneTimeSetUp]
        public void SetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinGroupFromTsv");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            if (Directory.Exists(_outputDirectory))
                Directory.Delete(_outputDirectory, true);
        }

        private static ProteinGroupFromTsv Row(string accession) =>
            new ProteinGroupFromTsvFile(FixturePath).Single(r => r.ProteinGroupName == accession);

        [Test]
        public void ReadsEveryRowThroughBothEntryPoints()
        {
            var direct = new ProteinGroupFromTsvFile(FixturePath);
            var factory = FileReader.ReadFile<ProteinGroupFromTsvFile>(FixturePath);
            Assert.That(direct.Count(), Is.EqualTo(6));
            Assert.That(factory.Count(), Is.EqualTo(6));
            Assert.That(direct.CanRead(FixturePath));
            Assert.That(direct.FileType, Is.EqualTo(SupportedFileType.MetaMorpheusQuantifiedProteinGroups));
            Assert.That(direct.Software, Is.EqualTo(Software.MetaMorpheus));
        }

        /// <summary>MetaMorpheus 1.1.x wrote "Protein Accession", "Protein QValue", "Best Peptide PEP"; mzLib's
        /// schema writes "BioPolymer Accession" and has no PEP column. Both must read.</summary>
        [Test]
        public void ReadsTheFixedColumnsInMetaMorpheus11Vocabulary()
        {
            var row = Row("P68363");
            Assert.Multiple(() =>
            {
                Assert.That(row.Gene, Is.EqualTo("TUBA1B"));
                Assert.That(row.NumberOfPsms, Is.EqualTo(268));
                Assert.That(row.NumberOfMembers, Is.EqualTo(1));
                Assert.That(row.CumulativeTarget, Is.EqualTo(1));
                Assert.That(row.QValue, Is.EqualTo(0));
                Assert.That(row.BestPep, Is.EqualTo(0));
                Assert.That(row.DecoyContaminantTarget, Is.EqualTo("T"));
            });
        }

        /// <summary>Labels are kept verbatim, "-calib" included, and a blank intensity is null, not 0.</summary>
        [Test]
        public void SampleGroupsAreKeyedByVerbatimLabelAndBlankIsNull()
        {
            var row = Row("P68363");
            Assert.That(row.SampleGroups, Has.Count.EqualTo(18));
            var a = row.SampleGroups["QE-002106_GM1_a-calib"];
            Assert.That(a.SpectralCount, Is.EqualTo(10));
            Assert.That(a.Intensity, Is.EqualTo(1838317.1674502683));
            Assert.That(row.SampleGroups["QE-002112_GM2_c-calib"].Intensity, Is.Null);
        }

        /// <summary>"N6,N6,N6-trimethyllysine on K" is one modification. Splitting on the first comma
        /// would read four fields from a three-field record.</summary>
        [Test]
        public void AModificationNameWithCommasIsOneName()
        {
            var site = Row("P05141").SampleGroups["QE-002106_GM1_a-calib"].CountOccupancy.Sites.Single(s => s.Position == 52);
            Assert.That(site.ModificationIdWithMotif, Is.EqualTo("N6,N6,N6-trimethyllysine on K"));
            Assert.That((site.Numerator, site.Denominator, site.Fraction), Is.EqualTo((1.0, 1.0, 1.0)));
        }

        /// <summary>Entities are "|"-joined within one cell, and the accessions are "|"-joined in their own column.
        /// The two are not zipped: the cell does not name its entities.</summary>
        [Test]
        public void AMultiMemberGroupKeepsItsEntitiesApart()
        {
            var row = Row("P0C0S5|Q71UI9");
            Assert.That(row.Accessions, Is.EqualTo(new[] { "P0C0S5", "Q71UI9" }));
            var cell = row.SampleGroups["QE-002108_GM1_c-calib"].CountOccupancy;
            Assert.That(cell.Entities, Has.Count.EqualTo(2));
            Assert.That(cell.Entities.All(e => e.Single().Position == 72));
        }

        [Test]
        public void ContaminantsAndDecoysAreRowsThatSayWhatTheyAre()
        {
            var contaminant = Row("P02769");
            var decoy = Row("DECOY_P62750");
            Assert.That(contaminant.IsContaminant && !contaminant.IsDecoy);
            Assert.That(decoy.IsDecoy && !decoy.IsContaminant);
            Assert.That(decoy.QValue, Is.EqualTo(0.0011299435028248588));
            Assert.That(decoy.SampleGroups.Values.All(s => s.Intensity == null));
        }

        /// <summary>Two sites in one entity, and intensity pairs written at four significant digits in
        /// exponent form ("3.772E+06/3.024E+07").</summary>
        [Test]
        public void IntensityOccupancyReadsExponentPairs()
        {
            var cell = Row("P02769").SampleGroups["QE-002106_GM1_a-calib"].IntensityOccupancy;
            Assert.That(cell.Entities.Single(), Has.Count.EqualTo(2));
            var k = cell.Sites.First();
            Assert.That((k.Position, k.ModificationIdWithMotif), Is.EqualTo((304, "Hydroxylation on K")));
            Assert.That((k.Fraction, k.Numerator, k.Denominator), Is.EqualTo((0.1247, 3.772e6, 3.024e7)));
        }

        [Test]
        public void PositionZeroIsTheNTerminus()
        {
            var site = Row("P63104").SampleGroups["QE-002119_GM6_b-calib"].CountOccupancy.Sites.Single();
            Assert.That(site.Position, Is.EqualTo(0));
            Assert.That(site.IsNTerminus);
            Assert.That(site.ModificationIdWithMotif, Is.EqualTo("N-acetylmethionine on M"));
        }

        /// <summary>Derived values are computed once per row, not on every read.</summary>
        [Test]
        public void ParsedValuesAreKeptAcrossReads()
        {
            var row = Row("P0C0S5|Q71UI9");
            var group = row.SampleGroups["QE-002108_GM1_c-calib"];
            Assert.That(group.CountOccupancy, Is.SameAs(group.CountOccupancy));
            Assert.That(group.IntensityOccupancy, Is.SameAs(group.IntensityOccupancy));
            Assert.That(row.Accessions, Is.SameAs(row.Accessions));

            row.ProteinGroupName = "P1|P2|P3";
            Assert.That(row.Accessions, Is.EqualTo(new[] { "P1", "P2", "P3" }), "a new name is split again");
        }

        [Test]
        public void WritingIsRefused()
        {
            var file = new ProteinGroupFromTsvFile(FixturePath);
            Assert.Throws<NotSupportedException>(() => file.WriteResults(Path.Combine(_outputDirectory, "x")));
        }

        [Test]
        public void AFileThatIsNotAProteinGroupTableIsReportedWithItsPath()
        {
            string path = Path.Combine(_outputDirectory, "Broken_AllQuantifiedProteinGroups.tsv");
            File.WriteAllText(path, "Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\nP1\tT\tnot-a-number\n");
            var ex = Assert.Throws<MzLibException>(() => new ProteinGroupFromTsvFile(path).LoadResults());
            Assert.That(ex!.Message, Does.Contain(path));
        }

        [TestCase("T", false, false)]
        [TestCase("D", true, false)]
        [TestCase("ET", false, true)]
        [TestCase("ED", true, true)]
        public void EntrapmentGroupLabelsAreRead(string label, bool isDecoy, bool isEntrapment)
        {
            string path = Path.Combine(_outputDirectory, $"Entrapment{label}_AllProteinGroups.tsv");
            File.WriteAllText(path, $"Protein Accession\tProtein Decoy/Contaminant/Target\tProtein QValue\nRandom_P1_f0\t{label}\t0.001\n");

            var row = new ProteinGroupFromTsvFile(path).Single();
            Assert.That((row.IsDecoy, row.IsEntrapment, row.IsContaminant), Is.EqualTo((isDecoy, isEntrapment, false)));
        }

        /// <summary>
        /// The reader is the inverse of mzLib's own writer. Groups are rendered by
        /// <see cref="BioPolymerGroupTsvSchema"/> (the "BioPolymer ..." vocabulary) and read back.
        /// </summary>
        [Test]
        public void RoundTripsTheSchemaWriterLabelFree()
        {
            var group = BuildGroup([FileA, FileA, FileB], [100.0, 300.0, 50.0]);
            var fileA = new SpectraFileInfo(FileA, "Control", 0, 0, 0);
            var fileB = new SpectraFileInfo(FileB, "Treatment", 0, 0, 0);
            group.SamplesForQuantification = [fileA, fileB];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { fileA, 1000.0 }, { fileB, 2000.0 } };

            var row = WriteAndRead(group, "LabelFree");
            Assert.That(row.ProteinGroupName, Is.EqualTo("P00001"));
            Assert.That((row.QValue, row.BestScore, row.BestNotchQValue), Is.EqualTo((0.01, 10.0, 0.005)));
            Assert.That(row.BestPep, Is.Null, "the schema writes no PEP column");
            Assert.That(row.SampleGroups.Keys, Is.EquivalentTo(new[] { "goldenA", "goldenB" }));

            var a = row.SampleGroups["goldenA"];
            Assert.That((a.SpectralCount, a.Intensity), Is.EqualTo(((int?)2, (double?)1000.0)));
            var count = a.CountOccupancy.Sites.Single();
            Assert.That((count.Position, count.ModificationIdWithMotif, count.Numerator, count.Denominator),
                Is.EqualTo((3, "Phosphorylation on D", 1.0, 2.0)));
            Assert.That(a.IntensityOccupancy.Sites.Single().Fraction, Is.EqualTo(0.25));

            var b = row.SampleGroups["goldenB"];
            Assert.That(b.CountOccupancyText, Is.Null, "an empty occupancy cell is not reported");
            Assert.That(b.CountOccupancy.Sites, Is.Empty);
        }

        /// <summary>Isobaric: one counting pair per file, one intensity pair per channel. The counting label
        /// and the channel labels differ, so entries carry one half or the other, never a guessed pairing.</summary>
        [Test]
        public void RoundTripsTheSchemaWriterIsobaric()
        {
            var group = BuildGroup([FileA, FileA, FileA], [100.0, 300.0, 50.0]);
            var c126 = new IsobaricQuantSampleInfo(FileA, "Control", 1, 1, 0, 1, "126", 126.0, false);
            var c127 = new IsobaricQuantSampleInfo(FileA, "Control", 1, 1, 0, 2, "127N", 127.0, false);
            group.SamplesForQuantification = [c126, c127];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { c126, 500.0 }, { c127, 750.0 } };

            var row = WriteAndRead(group, "Isobaric");
            Assert.That(row.SampleGroups.Keys, Is.EquivalentTo(new[] { "goldenA", "goldenA_126", "goldenA_127N" }));
            Assert.That(row.SampleGroups["goldenA"].SpectralCount, Is.EqualTo(3));
            Assert.That(row.SampleGroups["goldenA"].Intensity, Is.Null);
            Assert.That(row.SampleGroups["goldenA_127N"].Intensity, Is.EqualTo(750.0));
            Assert.That(row.SampleGroups["goldenA_127N"].SpectralCount, Is.Null);
        }

        private const string FileA = @"C:\goldenA.raw";
        private const string FileB = @"C:\goldenB.raw";

        private ProteinGroupFromTsv WriteAndRead(BioPolymerGroup group, string name)
        {
            string path = Path.Combine(_outputDirectory, $"{name}_AllQuantifiedProteinGroups.tsv");
            File.WriteAllLines(path, [GroupTsv.Header(group), GroupTsv.Row(group)]);
            return FileReader.ReadFile<ProteinGroupFromTsvFile>(path).Single();
        }

        /// <summary>The group <see cref="BioPolymerGroupTsvGoldenTests"/> locks: one protein, a phosphorylated
        /// and an unmodified ACDEF, and a unique GHIK.</summary>
        private static BioPolymerGroup BuildGroup(string[] psmFilePaths, double[] psmIntensities)
        {
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001", organism: "Homo sapiens", name: "TestName",
                fullName: "Test Protein", geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });
            ModificationMotif.TryGetMotif("D", out var motif);
            var phospho = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);
            var modifiedForm = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5,
                new Dictionary<int, Modification> { { 4, phospho } });
            var unmodifiedForm = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
            var uniqueForm = new MockBioPolymerWithSetMods("GHIK", "GHIK", protein, 6, 9);
            var group = new BioPolymerGroup(new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { modifiedForm, uniqueForm },
                new HashSet<IBioPolymerWithSetMods> { uniqueForm });
            var psms = new List<global::Omics.SpectralMatch.ISpectralMatch>
            {
                new MockSpectralMatch(psmFilePaths[0], "ACD[Phosphorylation]EF", "ACDEF", 10.0, 1, new[] { modifiedForm }),
                new MockSpectralMatch(psmFilePaths[1], "ACDEF", "ACDEF", 8.0, 2, new[] { unmodifiedForm }),
                new MockSpectralMatch(psmFilePaths[2], "GHIK", "GHIK", 6.0, 3, new[] { uniqueForm })
            };
            for (int i = 0; i < psms.Count; i++)
                ((MockSpectralMatch)psms[i]).Intensities = [psmIntensities[i]];
            group.AllPsmsBelowOnePercentFDR = [.. psms];
            group.CumulativeTarget = 7;
            group.CumulativeDecoy = 2;
            group.QValue = 0.01;
            group.BestBioPolymerWithSetModsScore = 10.0;
            group.BestBioPolymerWithSetModsQValue = 0.005;
            group.CalculateSequenceCoverage();
            return group;
        }
    }
}
