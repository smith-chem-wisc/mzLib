using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Locks the exact tab-separated header and row emitted by <see cref="BioPolymerGroup"/>.
    /// Other tests in this folder assert with substrings; these compare every field, so a change
    /// to column order, count, naming, or value formatting fails here rather than silently
    /// reshaping output files.
    ///
    /// To regenerate after an intentional column change, run the Explicit DumpGoldenStrings test.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class BioPolymerGroupTsvGoldenTests
    {
        private const string FileA = @"C:\goldenA.raw";
        private const string FileB = @"C:\goldenB.raw";

        private static readonly string[] StaticHeader =
        [
            "BioPolymer Accession", "Gene", "Organism", "BioPolymer Full Name",
            "BioPolymer Unmodified Mass", "Number of BioPolymers in Group",
            "Unique Sequences", "Shared Sequences", "Number of Sequences",
            "Number of Unique Sequences", "Sequence Coverage Fraction", "Sequence Coverage",
            "Sequence Coverage with Mods", "Fragment Sequence Coverage"
        ];

        private static readonly string[] TrailingHeader =
        [
            "Number of PSMs", "BioPolymer Decoy/Contaminant/Target",
            "BioPolymer Cumulative Target", "BioPolymer Cumulative Decoy",
            "BioPolymer QValue", "Best Sequence Score", "Best Sequence Notch QValue"
        ];

        private static readonly string[] StaticRow =
        [
            "P00001", "GENE1", "Homo sapiens", "Test Protein", "NaN", "1",
            "GHIK", "ACDEF", "2", "1", "1", "ACDEFGHIK",
            "ACD[Phosphorylation on D]EFGHIK", "acdefghik"
        ];

        private static readonly string[] TrailingRow =
        [
            "3", "T", "7", "2", "0.01", "10", "0.005"
        ];

        private const string CountOccupancy = "pos3[Phosphorylation on D,info:fraction=0.50(1/2)]";
        private const string IntensityOccupancy = "pos3[Phosphorylation on D,info:fraction=0.2500(100/400)]";

        /// <summary>
        /// Builds a deterministic group: one protein, one shared modified sequence, one unique
        /// sequence, and three PSMs (modified ACDEF, unmodified ACDEF, and GHIK).
        /// Each sequence set yields a single element per pipe-joined field, so HashSet iteration
        /// order cannot affect the output being locked.
        /// </summary>
        private static BioPolymerGroup BuildGroup(string[] psmFilePaths, double[] psmIntensities = null)
        {
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001",
                organism: "Homo sapiens",
                name: "TestName",
                fullName: "Test Protein",
                geneNames: new List<Tuple<string, string>> { new("primary", "GENE1") });

            ModificationMotif.TryGetMotif("D", out var motif);
            var phospho = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);

            var modifiedForm = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", protein, 1, 5,
                new Dictionary<int, Modification> { { 4, phospho } });
            var unmodifiedForm = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
            var uniqueForm = new MockBioPolymerWithSetMods("GHIK", "GHIK", protein, 6, 9);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { modifiedForm, uniqueForm },
                new HashSet<IBioPolymerWithSetMods> { uniqueForm });

            var psms = new List<ISpectralMatch>
            {
                new MockSpectralMatch(psmFilePaths[0], "ACD[Phosphorylation]EF", "ACDEF", 10.0, 1, new[] { modifiedForm }),
                new MockSpectralMatch(psmFilePaths[1], "ACDEF", "ACDEF", 8.0, 2, new[] { unmodifiedForm }),
                new MockSpectralMatch(psmFilePaths[2], "GHIK", "GHIK", 6.0, 3, new[] { uniqueForm })
            };

            if (psmIntensities != null)
            {
                for (int i = 0; i < psms.Count; i++)
                {
                    ((MockSpectralMatch)psms[i]).Intensities = [psmIntensities[i]];
                }
            }

            group.AllPsmsBelowOnePercentFDR = [.. psms];
            group.CumulativeTarget = 7;
            group.CumulativeDecoy = 2;
            group.QValue = 0.01;
            group.BestBioPolymerWithSetModsScore = 10.0;
            group.BestBioPolymerWithSetModsQValue = 0.005;
            group.CalculateSequenceCoverage();

            return group;
        }

        private static BioPolymerGroup BuildNoDesignGroup() =>
            BuildGroup([FileA, FileA, FileA]);

        private static BioPolymerGroup BuildLabelFreeGroup()
        {
            var group = BuildGroup([FileA, FileA, FileB], [100.0, 300.0, 50.0]);
            var fileA = new SpectraFileInfo(FileA, "Control", 0, 0, 0);
            var fileB = new SpectraFileInfo(FileB, "Treatment", 0, 0, 0);

            group.SamplesForQuantification = [fileA, fileB];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { fileA, 1000.0 },
                { fileB, 2000.0 }
            };
            return group;
        }

        private static BioPolymerGroup BuildIsobaricGroup()
        {
            var group = BuildGroup([FileA, FileA, FileA], [100.0, 300.0, 50.0]);
            var channel126 = new IsobaricQuantSampleInfo(FileA, "Control", 1, 1, 0, 1, "126", 126.0, false);
            var channel127 = new IsobaricQuantSampleInfo(FileA, "Control", 1, 1, 0, 2, "127N", 127.0, false);

            group.SamplesForQuantification = [channel126, channel127];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { channel126, 500.0 },
                { channel127, 750.0 }
            };
            return group;
        }

        private static string[] Compose(string[] middle) =>
            [.. StaticHeader, .. middle, .. TrailingHeader];

        private static string[] ComposeRow(string[] middle) =>
            [.. StaticRow, .. middle, .. TrailingRow];

        private static void AssertTsv(BioPolymerGroup group, string[] expectedHeader, string[] expectedRow)
        {
            Assert.Multiple(() =>
            {
                Assert.That(GroupTsv.Header(group).Split('\t'), Is.EqualTo(expectedHeader), "header");
                Assert.That(GroupTsv.Row(group).Split('\t'), Is.EqualTo(expectedRow), "row");
            });
        }

        /// <summary>
        /// No experimental design: PSMs are grouped by their own file path and only spectral
        /// count and count-based occupancy are emitted (no intensity columns).
        /// </summary>
        [Test]
        public void NoExperimentalDesign_HeaderAndRowAreUnchanged()
        {
            AssertTsv(BuildNoDesignGroup(),
                Compose(["SpectralCount_goldenA", "CountOccupancy_goldenA"]),
                ComposeRow(["3", CountOccupancy]));
        }

        /// <summary>
        /// Label-free with two conditions and intensities: four columns per sample group.
        /// goldenB's occupancy fields are empty because its only PSM carries no modification.
        /// </summary>
        [Test]
        public void LabelFreeWithIntensities_HeaderAndRowAreUnchanged()
        {
            AssertTsv(BuildLabelFreeGroup(),
                Compose([
                    "SpectralCount_goldenA", "Intensity_goldenA", "CountOccupancy_goldenA", "IntensityOccupancy_goldenA",
                    "SpectralCount_goldenB", "Intensity_goldenB", "CountOccupancy_goldenB", "IntensityOccupancy_goldenB"
                ]),
                ComposeRow([
                    "2", "1000", CountOccupancy, IntensityOccupancy,
                    "1", "2000", "", ""
                ]));
        }

        /// <summary>
        /// Isobaric (TMT/iTRAQ): one four-column block per channel, ordered by file then channel label.
        /// Spectral count is per file, so both channels report all three PSMs.
        /// </summary>
        [Test]
        public void Isobaric_HeaderAndRowAreUnchanged()
        {
            AssertTsv(BuildIsobaricGroup(),
                Compose([
                    "SpectralCount_goldenA_126", "Intensity_goldenA_126", "CountOccupancy_goldenA_126", "IntensityOccupancy_goldenA_126",
                    "SpectralCount_goldenA_127N", "Intensity_goldenA_127N", "CountOccupancy_goldenA_127N", "IntensityOccupancy_goldenA_127N"
                ]),
                ComposeRow([
                    "3", "500", CountOccupancy, IntensityOccupancy,
                    "3", "750", CountOccupancy, IntensityOccupancy
                ]));
        }

        /// <summary>
        /// The invariant the column-schema refactor exists to guarantee: header and row always
        /// carry the same number of fields, for every experimental design.
        /// </summary>
        [Test]
        public void HeaderAndRowFieldCountsAgree(
            [Values("none", "labelfree", "isobaric")] string design)
        {
            var group = design switch
            {
                "labelfree" => BuildLabelFreeGroup(),
                "isobaric" => BuildIsobaricGroup(),
                _ => BuildNoDesignGroup()
            };

            Assert.That(GroupTsv.Row(group).Split('\t'), Has.Length.EqualTo(GroupTsv.Header(group).Split('\t').Length),
                $"Header/row column count mismatch for '{design}' design.");
        }

        /// <summary>
        /// A group quantified in only some conditions must still fill the full width of the header.
        ///
        /// When each group rendered its own columns, a group with no measured intensity in a
        /// condition emitted two columns where the header had four, shifting every later field on
        /// that row. Deriving one schema from the whole dataset makes the widths agree by
        /// construction: the unquantified condition contributes empty fields instead of vanishing.
        /// </summary>
        [Test]
        public void GroupMissingIntensityInOneCondition_StillMatchesHeaderWidth()
        {
            var fileA = new SpectraFileInfo(FileA, "Control", 0, 0, 0);
            var fileB = new SpectraFileInfo(FileB, "Treatment", 0, 0, 0);
            var samples = new List<ISampleInfo> { fileA, fileB };

            // Quantified in both conditions.
            var quantifiedInBoth = BuildGroup([FileA, FileA, FileB], [100.0, 300.0, 50.0]);
            quantifiedInBoth.SamplesForQuantification = samples;
            quantifiedInBoth.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { fileA, 1000.0 },
                { fileB, 2000.0 }
            };

            // Quantified in Control only — no feature found in Treatment.
            var quantifiedInOne = BuildGroup([FileA, FileA, FileB], [100.0, 300.0, 50.0]);
            quantifiedInOne.SamplesForQuantification = samples;
            quantifiedInOne.IntensitiesBySample = new Dictionary<ISampleInfo, double>
            {
                { fileA, 500.0 }
            };

            BioPolymerGroup[] dataset = [quantifiedInBoth, quantifiedInOne];
            int headerWidth = GroupTsv.Header(dataset).Split('\t').Length;

            Assert.Multiple(() =>
            {
                Assert.That(GroupTsv.RowInDataset(quantifiedInBoth, dataset).Split('\t'),
                    Has.Length.EqualTo(headerWidth), "fully quantified group");
                Assert.That(GroupTsv.RowInDataset(quantifiedInOne, dataset).Split('\t'),
                    Has.Length.EqualTo(headerWidth), "partially quantified group");
            });

            // The Treatment intensity field is present but empty rather than absent.
            var fields = GroupTsv.RowInDataset(quantifiedInOne, dataset).Split('\t');
            var headers = GroupTsv.Header(dataset).Split('\t');
            int treatmentIntensity = Array.IndexOf(headers, "Intensity_goldenB");

            Assert.That(treatmentIntensity, Is.GreaterThan(-1), "Intensity_goldenB column should exist");
            Assert.That(fields[treatmentIntensity], Is.Empty);
        }

        /// <summary>
        /// Maintenance tool, not a check: prints current output so the constants above can be
        /// regenerated after a deliberate column change.
        /// </summary>
        [Test]
        [Explicit("Harvests current output to re-author the golden constants.")]
        public void DumpGoldenStrings()
        {
            foreach (var (name, group) in new (string, BioPolymerGroup)[]
                     {
                         ("NO_DESIGN", BuildNoDesignGroup()),
                         ("LABEL_FREE", BuildLabelFreeGroup()),
                         ("ISOBARIC", BuildIsobaricGroup())
                     })
            {
                TestContext.Out.WriteLine($"===={name}_HEADER====");
                TestContext.Out.WriteLine(GroupTsv.Header(group));
                TestContext.Out.WriteLine($"===={name}_ROW====");
                TestContext.Out.WriteLine(GroupTsv.Row(group));
            }
        }
    }
}
