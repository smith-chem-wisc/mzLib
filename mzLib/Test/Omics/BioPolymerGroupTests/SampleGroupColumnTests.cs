using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Tests for how the dataset's sample groups become quantification columns.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SampleGroupColumnTests
    {
        private static BioPolymerGroup Group(string accession, params string[] psmFiles)
        {
            var protein = new MockBioPolymer("ACDEFGHIK", accession);
            var form = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                new HashSet<IBioPolymerWithSetMods> { form },
                new HashSet<IBioPolymerWithSetMods> { form });

            group.AllPsmsBelowOnePercentFDR = [.. psmFiles.Select((f, i) =>
                (ISpectralMatch)new MockSpectralMatch(f, "ACDEF", "ACDEF", 10.0, i + 1, new[] { form }))];

            return group;
        }

        /// <summary>
        /// Sample group labels are not unique. SampleGroupBuilder names a group after its first file
        /// whenever conditions are undefined or the design looks like SILAC, so two files sharing a
        /// filename in different directories produce two sample groups with the same label.
        ///
        /// Each must still get its own column, carrying its own values. Collapsing them by label
        /// would drop the second group's counts and intensities without any sign in the output.
        /// </summary>
        [Test]
        public void SampleGroupsSharingALabelEachKeepTheirOwnColumns()
        {
            var run1 = new SpectraFileInfo(@"C:\run1\sample.raw", "Control", 0, 0, 0);
            var run2 = new SpectraFileInfo(@"C:\run2\sample.raw", "Treatment", 0, 0, 0);

            var group = Group("P00001", @"C:\run1\sample.raw", @"C:\run2\sample.raw", @"C:\run2\sample.raw");
            group.SamplesForQuantification = [run1, run2];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { run1, 111.0 }, { run2, 999.0 } };

            var headers = GroupTsv.Header(group).Split('\t');
            var fields = GroupTsv.Row(group).Split('\t');

            string Field(string columnName) => fields[Array.IndexOf(headers, columnName)];

            Assert.Multiple(() =>
            {
                Assert.That(fields, Has.Length.EqualTo(headers.Length));

                Assert.That(Field("SpectralCount_run1_sample"), Is.EqualTo("1"));
                Assert.That(Field("SpectralCount_run2_sample"), Is.EqualTo("2"));
                Assert.That(Field("Intensity_run1_sample"), Is.EqualTo("111"));
                Assert.That(Field("Intensity_run2_sample"), Is.EqualTo("999"));
            });
        }

        /// <summary>
        /// The case a single-group test cannot reach: two records observed in DIFFERENT subsets of
        /// same-named files.
        ///
        /// Matching sample groups to columns by position within a record is wrong, because a record
        /// only has sample groups for the files it was observed in. A record seen only in runB has
        /// runB at index 0, where a record seen in both has runA — so runB's counts get filed under
        /// runA's column. Matching on file identity is what makes this correct.
        /// </summary>
        [Test]
        public void RecordsObservedInDifferentFilesAreNotCrossAttributed()
        {
            const string RunA = @"C:\runA\sample.raw";
            const string RunB = @"C:\runB\sample.raw";

            // Seen in both files; seen only in runB.
            var inBoth = Group("P00001", RunA, RunA, RunB);
            var runBOnly = Group("P00002", RunB, RunB, RunB, RunB);

            BioPolymerGroup[] dataset = [inBoth, runBOnly];
            var headers = GroupTsv.Header(dataset).Split('\t');

            string CountFor(BioPolymerGroup g, string columnName)
            {
                var fields = GroupTsv.RowInDataset(g, dataset).Split('\t');
                return fields[Array.IndexOf(headers, columnName)];
            }

            Assert.Multiple(() =>
            {
                // Colliding labels are widened by parent directory, so the columns are tellable apart.
                Assert.That(headers, Does.Contain("SpectralCount_runA_sample"));
                Assert.That(headers, Does.Contain("SpectralCount_runB_sample"));

                Assert.That(CountFor(inBoth, "SpectralCount_runA_sample"), Is.EqualTo("2"));
                Assert.That(CountFor(inBoth, "SpectralCount_runB_sample"), Is.EqualTo("1"));

                // The regression: these were "4" and "" — runB's spectra filed under runA.
                Assert.That(CountFor(runBOnly, "SpectralCount_runA_sample"), Is.Empty,
                    "not observed in runA, so its runA cell must be blank");
                Assert.That(CountFor(runBOnly, "SpectralCount_runB_sample"), Is.EqualTo("4"));
            });
        }

        /// <summary>
        /// Column names must be unique within a file, or a consumer reading by name cannot tell two
        /// samples apart — and different tools disagree about how to mangle duplicates.
        /// </summary>
        [Test]
        public void CollidingLabelsProduceDistinctColumnNames()
        {
            var dataset = new[]
            {
                Group("P00001", @"C:\runA\sample.raw", @"C:\runB\sample.raw")
            };

            var headers = GroupTsv.Header(dataset).Split('\t');
            var quantHeaders = headers.Where(h => h.StartsWith("SpectralCount_")).ToArray();

            Assert.That(quantHeaders, Is.Unique);
            Assert.That(quantHeaders, Is.EqualTo(new[] { "SpectralCount_runA_sample", "SpectralCount_runB_sample" }));
        }

        /// <summary>
        /// A sample group spans fractions and technical replicates, which are routinely stored one
        /// per directory under the same file name. Keying the group's per-file collections by file
        /// name therefore collides: the intensity dictionary silently kept one of the two, and
        /// FilesInGroup threw outright, taking down report generation.
        /// </summary>
        [Test]
        public void FractionsSharingAFileNameDoNotCollide()
        {
            var protein = new MockBioPolymer("ACDEFGHIK", "P00001");
            var form = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", protein, 1, 5);
            var group = new BioPolymerGroup([protein], [form], [form]);

            // (path, condition, biorep, techrep, fraction) — one sample, two fractions, same name.
            var fraction1 = new SpectraFileInfo(@"C:\frac1\sample.raw", "", 0, 0, 0);
            var fraction2 = new SpectraFileInfo(@"C:\frac2\sample.raw", "", 0, 0, 1);

            group.AllPsmsBelowOnePercentFDR =
            [
                new MockSpectralMatch(fraction1.FullFilePathWithExtension, "ACDEF", "ACDEF", 10.0, 1, new[] { form }),
                new MockSpectralMatch(fraction2.FullFilePathWithExtension, "ACDEF", "ACDEF", 10.0, 2, new[] { form })
            ];
            group.SamplesForQuantification = [fraction1, fraction2];
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { fraction1, 111.0 }, { fraction2, 999.0 } };

            Assert.DoesNotThrow(() => group.PopulateSampleGroupResults());

            var result = group.SampleGroupResults.Single();
            Assert.Multiple(() =>
            {
                Assert.That(result.FilesInGroup, Has.Count.EqualTo(2), "both fractions belong to the group");
                Assert.That(result.Intensity, Is.EqualTo(1110.0), "fraction intensities sum rather than overwrite");
            });
        }

        /// <summary>
        /// Intensity-based stoichiometry needs measured intensity behind it. A site with none has a
        /// zero denominator, and printing it produced "fraction=0.0000(0/0)" — a measured-looking
        /// zero — in rows whose Intensity cell was blank, so the two columns contradicted each other.
        /// Count-based occupancy is unaffected: its denominator is the observation count.
        /// </summary>
        [Test]
        public void IntensityOccupancyIsBlankWithoutMeasuredIntensity()
        {
            ModificationMotif.TryGetMotif("D", out var motif);
            var phospho = new Modification("Phospho", null, "Biological", null, motif, "Anywhere.", null, 79.966);
            var file = new SpectraFileInfo(@"C:\d\a.raw", "Control", 0, 0, 0);

            BioPolymerGroup Modified(string accession, double? psmIntensity, double? sampleIntensity)
            {
                var protein = new MockBioPolymer("ACDEFGHIK", accession);
                var form = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phospho]EF", protein, 1, 5,
                    new Dictionary<int, Modification> { { 4, phospho } });

                var group = new BioPolymerGroup([protein], [form], [form]);
                var psm = new MockSpectralMatch(file.FullFilePathWithExtension, "ACD[Phospho]EF", "ACDEF", 10.0, 1, new[] { form });
                if (psmIntensity.HasValue) psm.Intensities = [psmIntensity.Value];

                group.AllPsmsBelowOnePercentFDR = [psm];
                group.SamplesForQuantification = [file];
                if (sampleIntensity.HasValue)
                    group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { file, sampleIntensity.Value } };

                return group;
            }

            // One group carries intensity, so the dataset gets intensity columns; the other has none.
            var quantified = Modified("P00001", 500.0, 1000.0);
            var unquantified = Modified("P00002", null, null);

            BioPolymerGroup[] dataset = [quantified, unquantified];
            var headers = GroupTsv.Header(dataset).Split('\t');

            string Field(BioPolymerGroup g, string column) =>
                GroupTsv.RowInDataset(g, dataset).Split('\t')[Array.IndexOf(headers, column)];

            Assert.Multiple(() =>
            {
                Assert.That(Field(unquantified, "Intensity_a"), Is.Empty);
                Assert.That(Field(unquantified, "IntensityOccupancy_a"), Is.Empty,
                    "no measured intensity means no stoichiometry to report");

                // Count-based occupancy is still reported — the observation is real.
                Assert.That(Field(unquantified, "CountOccupancy_a"), Does.Contain("1/1"));

                // The quantified group is untouched.
                Assert.That(Field(quantified, "IntensityOccupancy_a"), Does.Contain("Phospho"));
            });
        }

        /// <summary>
        /// Widening is itself a source of collisions: the name a widened column takes can be the
        /// plain label of a different file that was never part of that collision. Resolving each
        /// label group in isolation therefore reintroduces the duplicate it set out to remove.
        /// </summary>
        [Test]
        public void WidenedNameDoesNotCollideWithAnotherFilesPlainLabel()
        {
            // run1\sample.raw and run2\sample.raw collide and widen to run1_sample / run2_sample.
            // x\run1_sample.raw is not part of that collision and would otherwise stay run1_sample.
            var dataset = new[]
            {
                Group("P00001", @"C:\run1\sample.raw", @"C:\run2\sample.raw", @"C:\x\run1_sample.raw")
            };

            var quantHeaders = GroupTsv.Header(dataset).Split('\t')
                .Where(h => h.StartsWith("SpectralCount_"))
                .ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(quantHeaders, Is.Unique);
                Assert.That(quantHeaders, Has.Length.EqualTo(3));
            });
        }

        /// <summary>
        /// Widening applies only where it is needed: labels that were already unique are untouched,
        /// so output does not change for datasets that were never ambiguous.
        /// </summary>
        [Test]
        public void NonCollidingLabelsAreLeftAlone()
        {
            var dataset = new[]
            {
                Group("P00001", @"C:\runA\alpha.raw", @"C:\runB\beta.raw")
            };

            var headers = GroupTsv.Header(dataset).Split('\t');

            Assert.Multiple(() =>
            {
                Assert.That(headers, Does.Contain("SpectralCount_alpha"));
                Assert.That(headers, Does.Contain("SpectralCount_beta"));
            });
        }

        /// <summary>
        /// The schema is a snapshot of the dataset's experimental design, so it has to be built
        /// after quantification state is final. Changing samples afterwards cannot corrupt row
        /// width — rows stay aligned to the schema they are rendered against — but the new sample
        /// simply has no column, so build the schema last.
        /// </summary>
        [Test]
        public void RowStaysAlignedWhenQuantificationChangesAfterSchemaIsBuilt()
        {
            var fileA = new SpectraFileInfo(@"C:\a.raw", "Control", 0, 0, 0);
            var fileB = new SpectraFileInfo(@"C:\b.raw", "Treatment", 0, 0, 0);

            var group = Group("P00001", @"C:\a.raw");
            group.SamplesForQuantification = [fileA];

            var schema = BioPolymerGroupTsvSchema.For([group]);
            int headerWidth = TsvWriter.HeaderLine(schema).Split('\t').Length;

            group.SamplesForQuantification = [fileA, fileB];

            Assert.That(TsvWriter.RowLine(schema, group).Split('\t'), Has.Length.EqualTo(headerWidth));
        }
    }
}
