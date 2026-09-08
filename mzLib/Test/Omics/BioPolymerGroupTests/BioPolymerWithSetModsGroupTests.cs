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
    /// Tests for the digestion-product-level group: base-sequence identity, peptide-local
    /// occupancy, and its TSV schema.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class BioPolymerWithSetModsGroupTests
    {
        private const string File1 = @"C:\pepgroupA.raw";
        private const string Sequence = "ACDEFK";

        private static MockBioPolymer Parent(string accession = "P00001") =>
            new("MMMACDEFKGGG", accession);

        private static Modification Mod(string residue, string id, string location)
        {
            ModificationMotif.TryGetMotif(residue, out var motif);
            return new Modification(id, null, "Biological", null, motif, location, null, 79.966);
        }

        /// <summary>
        /// ACDEFK sits at residues 4-9 of the parent. Peptide-local mod keys follow the
        /// AllModsOneIsNterminus convention: 1 is the N-terminus, 2 is residue A, and so on.
        /// </summary>
        private static MockBioPolymerWithSetMods Form(
            IBioPolymer parent, string fullSequence, Dictionary<int, Modification> mods = null) =>
            new(Sequence, fullSequence, parent, 4, 9, mods);

        private static MockSpectralMatch Psm(
            string fullSequence, IBioPolymerWithSetMods form, int scan, double? intensity = null)
        {
            var psm = new MockSpectralMatch(File1, fullSequence, Sequence, 10.0, scan, new[] { form });
            if (intensity.HasValue)
                psm.Intensities = [intensity.Value];
            return psm;
        }

        #region Identity and construction

        [Test]
        public void DistinctPeptidoformsOfOneSequenceAreKeptTogether()
        {
            var parent = Parent();
            var unmodified = Form(parent, "ACDEFK");
            var phosphorylated = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            var group = new BioPolymerWithSetModsGroup(Sequence, [unmodified, phosphorylated]);

            Assert.Multiple(() =>
            {
                Assert.That(group.BaseSequence, Is.EqualTo(Sequence));
                Assert.That(group.Peptidoforms, Has.Count.EqualTo(2),
                    "both forms should survive; dedup is by FullSequence, not by the element's Equals");
                Assert.That(group.ParentBioPolymers, Has.Count.EqualTo(1));
            });
        }

        [Test]
        public void RepeatedPeptidoformIsRecordedOnce()
        {
            var parent = Parent();
            var group = new BioPolymerWithSetModsGroup(Sequence,
                [Form(parent, "ACDEFK"), Form(parent, "ACDEFK")]);

            Assert.That(group.Peptidoforms, Has.Count.EqualTo(1));
        }

        [Test]
        public void SequenceSharedAcrossParentsGivesOneGroupWithBothParents()
        {
            var group = new BioPolymerWithSetModsGroup(Sequence,
                [Form(Parent("P00002"), "ACDEFK"), Form(Parent("P00001"), "ACDEFK")]);

            Assert.That(group.ListOfParentsOrderedByAccession.Select(p => p.Accession).ToArray(),
                Is.EqualTo(new[] { "P00001", "P00002" }), "parents should be accession-ordered");
        }

        [Test]
        public void FormWithDifferentBaseSequenceIsRejected()
        {
            var parent = Parent();
            var wrongForm = new MockBioPolymerWithSetMods("MMMMM", "MMMMM", parent, 1, 5);

            var ex = Assert.Throws<ArgumentException>(
                () => new BioPolymerWithSetModsGroup(Sequence, [wrongForm]));

            Assert.That(ex.Message, Does.Contain(Sequence).And.Contain("MMMMM"),
                "message should name both the group's sequence and the offending one");
        }

        [Test]
        public void PsmForDifferentSequenceIsRejected()
        {
            var parent = Parent();
            var group = new BioPolymerWithSetModsGroup(Sequence, [Form(parent, "ACDEFK")]);
            var foreignForm = new MockBioPolymerWithSetMods("MMMMM", "MMMMM", parent, 1, 5);
            var foreignPsm = new MockSpectralMatch(File1, "MMMMM", "MMMMM", 5.0, 9, new[] { foreignForm });

            Assert.Throws<ArgumentException>(
                () => group.AllPsmsBelowOnePercentFDR = [foreignPsm]);
        }

        [Test]
        public void AmbiguousPsmIsAccepted()
        {
            var parent = Parent();
            var group = new BioPolymerWithSetModsGroup(Sequence, [Form(parent, "ACDEFK")]);
            var ambiguous = new MockSpectralMatch(File1, null, null, 5.0, 9, Array.Empty<IBioPolymerWithSetMods>());

            Assert.DoesNotThrow(() => group.AllPsmsBelowOnePercentFDR = [ambiguous],
                "a PSM with no resolved base sequence should not be treated as a mismatch");
        }

        [Test]
        public void GroupsAreEqualOnBaseSequence()
        {
            var a = new BioPolymerWithSetModsGroup(Sequence, [Form(Parent("P00001"), "ACDEFK")]);
            var b = new BioPolymerWithSetModsGroup(Sequence, [Form(Parent("P00002"), "ACD[Phospho]EFK")]);

            Assert.Multiple(() =>
            {
                Assert.That(a, Is.EqualTo(b));
                Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
            });
        }

        [Test]
        public void CreateGroupsBucketsPsmsByBaseSequence()
        {
            var parent = Parent();
            var acdefk = Form(parent, "ACDEFK");
            var other = new MockBioPolymerWithSetMods("GGG", "GGG", parent, 10, 12);

            var groups = BioPolymerWithSetModsGroup.CreateGroups([
                Psm("ACDEFK", acdefk, 1),
                Psm("ACDEFK", acdefk, 2),
                new MockSpectralMatch(File1, "GGG", "GGG", 4.0, 3, new[] { other })
            ]);

            Assert.Multiple(() =>
            {
                Assert.That(groups, Has.Count.EqualTo(2));
                Assert.That(groups.Single(g => g.BaseSequence == Sequence).AllPsmsBelowOnePercentFDR,
                    Has.Count.EqualTo(2));
            });
        }

        /// <summary>
        /// An unresolved PSM has an empty base sequence rather than null, so it must be filtered
        /// out explicitly — otherwise it forms a group keyed on the empty string.
        /// </summary>
        [Test]
        public void CreateGroupsSkipsPsmsWithNoResolvedSequence()
        {
            var parent = Parent();
            var acdefk = Form(parent, "ACDEFK");

            var groups = BioPolymerWithSetModsGroup.CreateGroups([
                Psm("ACDEFK", acdefk, 1),
                new MockSpectralMatch(File1, null, null, 3.0, 2, Array.Empty<IBioPolymerWithSetMods>())
            ]);

            Assert.That(groups.Select(g => g.BaseSequence).ToArray(), Is.EqualTo(new[] { Sequence }));
        }

        /// <summary>
        /// A group holding both a confident and an unresolved PSM must still compute occupancy;
        /// the unresolved one is not evidence for or against a modification at any position.
        /// </summary>
        [Test]
        public void UnresolvedPsmDoesNotBreakOccupancy()
        {
            var parent = Parent();
            var modified = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            var group = new BioPolymerWithSetModsGroup(Sequence, [modified])
            {
                AllPsmsBelowOnePercentFDR =
                [
                    Psm("ACD[Phospho]EFK", modified, 1),
                    new MockSpectralMatch(File1, null, null, 3.0, 2, Array.Empty<IBioPolymerWithSetMods>())
                ]
            };

            group.PopulateSampleGroupResults();
            var occupancy = group.SampleGroupResults.Single().DigestionProductOccupancy[Sequence];

            Assert.That(occupancy[4].Single().ToModInfoString(), Does.Contain("1/1"),
                "the unresolved PSM should not inflate the denominator");
        }

        #endregion

        #region Occupancy in peptide-local coordinates

        /// <summary>
        /// Positions are relative to the peptide, not the parent. ACDEFK starts at residue 4 of the
        /// parent, so a protein-level report would place these mods 3 residues higher.
        /// </summary>
        [Test]
        [TestCase(1, "pos0", TestName = "Occupancy_NTerminalModReportsPos0")]
        [TestCase(4, "pos3", TestName = "Occupancy_InternalModReportsPeptideResidue")]
        [TestCase(8, "pos7", TestName = "Occupancy_CTerminalModReportsPastLastResidue")]
        public void OccupancyUsesPeptideLocalPositions(int modKey, string expectedPosition)
        {
            var parent = Parent();
            var modified = Form(parent, $"modified{modKey}",
                new Dictionary<int, Modification> { { modKey, Mod("D", "Phospho", "Anywhere.") } });
            var unmodified = Form(parent, "ACDEFK");

            var group = new BioPolymerWithSetModsGroup(Sequence, [modified, unmodified])
            {
                AllPsmsBelowOnePercentFDR =
                [
                    Psm($"modified{modKey}", modified, 1),
                    Psm("ACDEFK", unmodified, 2)
                ]
            };

            group.PopulateSampleGroupResults();
            var occupancy = group.SampleGroupResults.Single().DigestionProductOccupancy[Sequence];

            Assert.Multiple(() =>
            {
                Assert.That(occupancy, Does.ContainKey(modKey));
                Assert.That(occupancy[modKey].Single().ToModInfoString(),
                    Does.StartWith(expectedPosition).And.Contain("1/2"),
                    "one of the two PSMs carries the mod");
            });
        }

        [Test]
        public void AllUnmodifiedYieldsNoOccupancyRatherThanZeroDivision()
        {
            var parent = Parent();
            var unmodified = Form(parent, "ACDEFK");

            var group = new BioPolymerWithSetModsGroup(Sequence, [unmodified])
            {
                AllPsmsBelowOnePercentFDR = [Psm("ACDEFK", unmodified, 1), Psm("ACDEFK", unmodified, 2)]
            };

            group.PopulateSampleGroupResults();

            Assert.That(group.SampleGroupResults.Single().DigestionProductOccupancy, Is.Empty);
        }

        [Test]
        public void PsmWhoseFormCannotBeResolvedStillCountsTowardTheDenominator()
        {
            var parent = Parent();
            var modified = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            // Second PSM reports a full sequence that matches none of the group's forms, so it
            // cannot contribute a modification — but it is still an observation of this sequence.
            var unresolvable = new MockSpectralMatch(File1, "ACD[Unknown]EFK", Sequence, 9.0, 2,
                Array.Empty<IBioPolymerWithSetMods>());

            var group = new BioPolymerWithSetModsGroup(Sequence, [modified])
            {
                AllPsmsBelowOnePercentFDR = [Psm("ACD[Phospho]EFK", modified, 1), unresolvable]
            };

            group.PopulateSampleGroupResults();
            var occupancy = group.SampleGroupResults.Single().DigestionProductOccupancy[Sequence];

            Assert.That(occupancy[4].Single().ToModInfoString(), Does.Contain("1/2"));
        }

        /// <summary>
        /// A sample group containing only unresolved PSMs has nothing to attribute a modification
        /// to, which is not an error. Occupancy runs per sample group, so with no experimental
        /// design a single unresolved PSM alone in its own file forms such a bucket — and throwing
        /// there would take down report generation for the whole file.
        /// </summary>
        [Test]
        public void SampleGroupOfOnlyUnresolvedPsmsYieldsNoOccupancy()
        {
            var parent = Parent();
            var modified = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            var group = new BioPolymerWithSetModsGroup(Sequence, [modified])
            {
                AllPsmsBelowOnePercentFDR =
                [
                    Psm("ACD[Phospho]EFK", modified, 1),
                    // Its own source file, so it lands in a bucket by itself.
                    new MockSpectralMatch(@"C:\other.raw", null, null, 3.0, 2, Array.Empty<IBioPolymerWithSetMods>())
                ]
            };

            Assert.DoesNotThrow(() => group.PopulateSampleGroupResults());

            var buckets = group.SampleGroupResults;
            Assert.Multiple(() =>
            {
                Assert.That(buckets, Has.Count.EqualTo(2), "one bucket per source file");
                Assert.That(buckets.Single(b => b.Label == "other").DigestionProductOccupancy, Is.Empty);
                Assert.That(buckets.Single(b => b.Label == "pepgroupA").DigestionProductOccupancy, Is.Not.Empty);
            });
        }

        /// <summary>
        /// Directly: an all-unresolved list is an empty result, not an exception.
        /// </summary>
        [Test]
        public void CalculatorReturnsEmptyForUnresolvedOnlyInput()
        {
            var unresolved = new MockSpectralMatch(@"C:\x.raw", null, null, 1.0, 1, Array.Empty<IBioPolymerWithSetMods>());

            Assert.Multiple(() =>
            {
                Assert.That(ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([unresolved]), Is.Empty);
                Assert.That(ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy([]), Is.Empty);
            });
        }

        /// <summary>
        /// The setter validates its input, so it must not then hold a reference the caller can keep
        /// mutating — otherwise a foreign PSM added afterwards is filed under this group's sequence
        /// and reports a modification at a residue that sequence does not have.
        /// </summary>
        [Test]
        public void PsmsAreCopiedSoLaterCallerMutationCannotBypassValidation()
        {
            var parent = Parent();
            var form = Form(parent, "ACDEFK");
            var callersSet = new HashSet<ISpectralMatch> { Psm("ACDEFK", form, 1) };

            var group = new BioPolymerWithSetModsGroup(Sequence, [form])
            {
                AllPsmsBelowOnePercentFDR = callersSet
            };

            var foreignForm = new MockBioPolymerWithSetMods("MMMMM", "MMMMM", parent, 1, 5);
            callersSet.Add(new MockSpectralMatch(File1, "MMMMM", "MMMMM", 5.0, 9, new[] { foreignForm }));

            Assert.Multiple(() =>
            {
                Assert.That(group.AllPsmsBelowOnePercentFDR, Has.Count.EqualTo(1),
                    "the group must not see the PSM added to the caller's set after assignment");
                Assert.That(group.AllPsmsBelowOnePercentFDR.Select(p => p.BaseSequence), Is.All.EqualTo(Sequence));
            });
        }

        [Test]
        public void EmptyPsmSetProducesNoSampleGroups()
        {
            var group = new BioPolymerWithSetModsGroup(Sequence, [Form(Parent(), "ACDEFK")]);

            group.PopulateSampleGroupResults();

            Assert.Multiple(() =>
            {
                Assert.That(group.SampleGroupResults, Is.Empty);
                Assert.That(() => GroupTsv.PeptideRow(group), Throws.Nothing);
            });
        }

        [Test]
        public void SinglePsmGivesFullOccupancy()
        {
            var parent = Parent();
            var modified = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            var group = new BioPolymerWithSetModsGroup(Sequence, [modified])
            {
                AllPsmsBelowOnePercentFDR = [Psm("ACD[Phospho]EFK", modified, 1)]
            };

            group.PopulateSampleGroupResults();
            var occupancy = group.SampleGroupResults.Single().DigestionProductOccupancy[Sequence];

            Assert.That(occupancy[4].Single().CountBasedOccupancy, Is.EqualTo(1.0));
        }

        [Test]
        public void ScoreTakesTheBestPsm()
        {
            var parent = Parent();
            var form = Form(parent, "ACDEFK");
            var group = new BioPolymerWithSetModsGroup(Sequence, [form]);

            group.AllPsmsBelowOnePercentFDR =
            [
                new MockSpectralMatch(File1, "ACDEFK", Sequence, 4.0, 1, new[] { form }),
                new MockSpectralMatch(File1, "ACDEFK", Sequence, 12.0, 2, new[] { form })
            ];
            group.Score();

            Assert.That(group.BestPsmScore, Is.EqualTo(12.0));
        }

        [Test]
        public void ScoreOnEmptyGroupIsZero()
        {
            var group = new BioPolymerWithSetModsGroup(Sequence, [Form(Parent(), "ACDEFK")]);

            group.Score();

            Assert.That(group.BestPsmScore, Is.Zero);
        }

        #endregion

        #region Schema

        private static BioPolymerWithSetModsGroup QuantifiedGroup(
            string accession, IReadOnlyDictionary<ISampleInfo, double> intensities, List<ISampleInfo> samples)
        {
            var parent = Parent(accession);
            var modified = Form(parent, "ACD[Phospho]EFK",
                new Dictionary<int, Modification> { { 4, Mod("D", "Phospho", "Anywhere.") } });

            return new BioPolymerWithSetModsGroup(Sequence, [modified])
            {
                AllPsmsBelowOnePercentFDR = [Psm("ACD[Phospho]EFK", modified, 1, 100.0)],
                SamplesForQuantification = samples,
                IntensitiesBySample = intensities.ToDictionary(kvp => kvp.Key, kvp => kvp.Value)
            };
        }

        [Test]
        public void SchemaOmitsParentLevelColumns()
        {
            var group = new BioPolymerWithSetModsGroup(Sequence, [Form(Parent(), "ACDEFK")]);
            var headers = GroupTsv.PeptideHeader(group).Split('\t');

            Assert.Multiple(() =>
            {
                Assert.That(headers, Does.Contain("Base Sequence"));
                Assert.That(headers, Does.Contain("Peptidoforms"));
                Assert.That(headers, Does.Not.Contain("Gene"));
                Assert.That(headers, Does.Not.Contain("Organism"));
                Assert.That(headers, Does.Not.Contain("Sequence Coverage"));
            });
        }

        [Test]
        public void ParentAccessionsAndResidueColumnsLineUp()
        {
            var group = new BioPolymerWithSetModsGroup(Sequence,
                [Form(Parent("P00002"), "ACDEFK"), Form(Parent("P00001"), "ACDEFK")]);

            var headers = GroupTsv.PeptideHeader(group).Split('\t');
            var fields = GroupTsv.PeptideRow(group).Split('\t');

            string Field(string header) => fields[Array.IndexOf(headers, header)];

            Assert.Multiple(() =>
            {
                Assert.That(Field("Parent Accessions"), Is.EqualTo("P00001|P00002"));
                Assert.That(Field("Start Residue in Parent"), Is.EqualTo("4|4"));
                Assert.That(Field("End Residue in Parent"), Is.EqualTo("9|9"));
            });
        }

        /// <summary>
        /// A sequence can occur more than once in the same protein — repeat domains in histones,
        /// collagens and mucins are ordinary, and repeats are exactly where site-level occupancy
        /// interpretation matters most. Deduplicating on (form, parent) alone kept whichever
        /// occurrence was seen first and dropped the rest from the provenance columns silently.
        /// </summary>
        [Test]
        public void SequenceRepeatedWithinOneParentReportsEveryOccurrence()
        {
            // ACDEFK sits at residues 4-9 and again at 13-18.
            var parent = new MockBioPolymer("MMMACDEFKGGGACDEFKGGG", "P00001");
            var atFirst = new MockBioPolymerWithSetMods(Sequence, Sequence, parent, 4, 9);
            var atSecond = new MockBioPolymerWithSetMods(Sequence, Sequence, parent, 13, 18);

            var group = new BioPolymerWithSetModsGroup(Sequence, [atFirst, atSecond]);

            var headers = GroupTsv.PeptideHeader(group).Split('\t');
            var fields = GroupTsv.PeptideRow(group).Split('\t');
            string Field(string header) => fields[Array.IndexOf(headers, header)];

            Assert.Multiple(() =>
            {
                Assert.That(group.Peptidoforms, Has.Count.EqualTo(2), "both occurrences are distinct observations");

                // The accession repeats so the three lists stay readable entry-for-entry.
                Assert.That(Field("Parent Accessions"), Is.EqualTo("P00001|P00001"));
                Assert.That(Field("Start Residue in Parent"), Is.EqualTo("4|13"));
                Assert.That(Field("End Residue in Parent"), Is.EqualTo("9|18"));

                // Still one distinct parent.
                Assert.That(Field("Number of Parents"), Is.EqualTo("1"));
            });
        }

        [Test]
        public void HeaderAndRowFieldCountsAgree()
        {
            var file = new SpectraFileInfo(File1, "Control", 0, 0, 0);
            var samples = new List<ISampleInfo> { file };
            var group = QuantifiedGroup("P00001",
                new Dictionary<ISampleInfo, double> { { file, 1000.0 } }, samples);

            Assert.That(GroupTsv.PeptideRow(group).Split('\t'),
                Has.Length.EqualTo(GroupTsv.PeptideHeader(group).Split('\t').Length));
        }

        /// <summary>
        /// The same dataset-level guarantee the parent-level report has: a group quantified in only
        /// some conditions still fills the header's width.
        /// </summary>
        [Test]
        public void GroupMissingIntensityInOneConditionStillMatchesHeaderWidth()
        {
            var fileA = new SpectraFileInfo(File1, "Control", 0, 0, 0);
            var fileB = new SpectraFileInfo(@"C:\pepgroupB.raw", "Treatment", 0, 0, 0);
            var samples = new List<ISampleInfo> { fileA, fileB };

            var inBoth = QuantifiedGroup("P00001",
                new Dictionary<ISampleInfo, double> { { fileA, 1000.0 }, { fileB, 2000.0 } }, samples);
            var inOne = QuantifiedGroup("P00002",
                new Dictionary<ISampleInfo, double> { { fileA, 500.0 } }, samples);

            BioPolymerWithSetModsGroup[] dataset = [inBoth, inOne];
            int headerWidth = GroupTsv.PeptideHeader(dataset).Split('\t').Length;

            Assert.Multiple(() =>
            {
                Assert.That(GroupTsv.PeptideRowInDataset(inBoth, dataset).Split('\t'), Has.Length.EqualTo(headerWidth));
                Assert.That(GroupTsv.PeptideRowInDataset(inOne, dataset).Split('\t'), Has.Length.EqualTo(headerWidth));
            });
        }

        [Test]
        public void OccupancyColumnCarriesPeptideLocalPosition()
        {
            var file = new SpectraFileInfo(File1, "Control", 0, 0, 0);
            var samples = new List<ISampleInfo> { file };
            var group = QuantifiedGroup("P00001",
                new Dictionary<ISampleInfo, double> { { file, 1000.0 } }, samples);

            var headers = GroupTsv.PeptideHeader(group).Split('\t');
            var fields = GroupTsv.PeptideRow(group).Split('\t');
            var occupancy = fields[Array.IndexOf(headers, "CountOccupancy_pepgroupA")];

            // Peptide residue 3, not parent residue 6.
            Assert.That(occupancy, Does.StartWith("pos3").And.Contain("Phospho"));
        }

        #endregion
    }
}
