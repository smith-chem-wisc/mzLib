using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymerGroup;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Tests for the column-name collision rule.
    ///
    /// A sample group's label is a display name and is not unique: <see cref="SampleGroupBuilder"/>
    /// names a group after its first file whenever there is no experimental design to name it by, so
    /// two files sharing a name in different directories produce two sample groups labelled the same.
    /// Emitting two columns with one name leaves a reader unable to tell which file a value came from.
    ///
    /// The widening half of the rule is reachable through the schema. What is pinned here is the half
    /// that widening cannot reach -- a label with no path to widen with, which is exactly what an
    /// experimental design produces, since its labels are built from condition and replicate rather
    /// than from a file. That path is the difference between two distinguishable columns and two
    /// columns a reader must guess between.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SampleGroupLabelsTests
    {
        private static (IReadOnlyList<string> Identities, IReadOnlyDictionary<string, (string, string?)> Labels)
            Input(params (string Identity, string Label, string? SourcePath)[] entries)
            => (entries.Select(e => e.Identity).ToList(),
                entries.ToDictionary(e => e.Identity, e => (e.Label, e.SourcePath)));

        /// <summary>
        /// The design names a sample "Control_1" whatever file it came from, so there is no path to
        /// widen with. Two such samples must still get two distinct column names.
        /// </summary>
        [Test]
        public void LabelsWithNoPathToWidenAreSeparatedByOrdinal()
        {
            var (identities, labels) = Input(
                ("Control|0", "Control_1", null),
                ("Control|1", "Control_1", null));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result["Control|0"], Is.EqualTo("Control_1"), "the first occurrence keeps its label");
                Assert.That(result["Control|1"], Is.EqualTo("Control_1_2"));
            });
        }

        /// <summary>
        /// Numbering starts at 2 and continues, so an n-way collision yields n distinct names rather
        /// than one name reused.
        /// </summary>
        [Test]
        public void EveryDuplicateAfterTheFirstIsNumbered()
        {
            var (identities, labels) = Input(
                ("a", "sample", null),
                ("b", "sample", null),
                ("c", "sample", null));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result.Values, Is.EquivalentTo(new[] { "sample", "sample_2", "sample_3" }));
                Assert.That(result.Values.Distinct().Count(), Is.EqualTo(3), "a column name may not repeat");
            });
        }

        /// <summary>
        /// The ordinal it would have taken can already be somebody else's real label. Suffixing
        /// blindly would produce the collision the rule exists to remove, so the ordinal advances
        /// past what is taken and the sample that owns the name keeps it.
        /// </summary>
        [Test]
        public void AnOrdinalAlreadyInUseIsSkippedRatherThanReused()
        {
            var (identities, labels) = Input(
                ("a", "sample", null),
                ("b", "sample_2", null),
                ("c", "sample", null));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result["b"], Is.EqualTo("sample_2"), "the sample whose real label this is keeps it");
                Assert.That(result["c"], Is.EqualTo("sample_3"), "so the duplicate advances past it");
                Assert.That(result.Values.Distinct().Count(), Is.EqualTo(3));
            });
        }

        /// <summary>
        /// A path with no directory component has no ancestor to widen with either, so it takes the
        /// same last resort as a design-supplied label.
        /// </summary>
        [Test]
        public void ABareFilenameHasNoParentToWidenWith()
        {
            var (identities, labels) = Input(
                ("a", "sample", "sample.raw"),
                ("b", "sample", "sample.raw"));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.That(result.Values, Is.EquivalentTo(new[] { "sample", "sample_2" }));
        }

        /// <summary>
        /// Widening is preferred to numbering whenever there is a path to widen with, because
        /// "Rep2_sample" says which file the column came from and "sample_2" does not.
        /// </summary>
        [Test]
        public void ACollisionWidensBeforeItSuffixes()
        {
            var (identities, labels) = Input(
                ("a", "sample", @"C:\Exp\Rep1\sample.raw"),
                ("b", "sample", @"C:\Exp\Rep2\sample.raw"));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result["a"], Is.EqualTo("Rep1_sample"));
                Assert.That(result["b"], Is.EqualTo("Rep2_sample"));
            });
        }

        /// <summary>
        /// Labels that do not collide are returned unchanged, so the datasets that were already
        /// unambiguous see no change in their column names.
        /// </summary>
        [Test]
        public void LabelsThatDoNotCollideAreLeftAlone()
        {
            var (identities, labels) = Input(
                ("a", "Control_1", @"C:\Exp\Rep1\a.raw"),
                ("b", "Treated_1", @"C:\Exp\Rep2\b.raw"));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result["a"], Is.EqualTo("Control_1"));
                Assert.That(result["b"], Is.EqualTo("Treated_1"));
            });
        }

        /// <summary>
        /// Widening is bounded, and the bound is reachable in principle: two paths that agree for
        /// more ancestors than the rule will walk. The guarantee that survives is the one that
        /// matters -- it terminates, and the names it returns are still distinct.
        /// </summary>
        [Test]
        public void ACollisionDeeperThanTheWideningBoundStillTerminatesWithDistinctNames()
        {
            // 66 shared ancestors, differing only above them: deeper than the 64 rounds of widening.
            string shared = string.Concat(Enumerable.Repeat(@"a\", 66));
            var (identities, labels) = Input(
                ("a", "sample", $@"C:\p\{shared}sample.raw"),
                ("b", "sample", $@"C:\q\{shared}sample.raw"));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(2));
                Assert.That(result.Values.Distinct().Count(), Is.EqualTo(2), "two columns may not share a name");
                Assert.That(result.Values, Has.All.Not.Empty);
            });
        }

        /// <summary>
        /// A file sitting at the root of a drive has a containing directory, but that directory has
        /// no name of its own. Widening has to stop there rather than prefix an empty segment, which
        /// would produce a column named "_sample" and still collide.
        /// </summary>
        [Test]
        public void AFileAtADriveRootHasNoNamedAncestorToWidenWith()
        {
            var (identities, labels) = Input(
                ("a", "sample", @"C:\sample.raw"),
                ("b", "sample", @"D:\sample.raw"));

            var result = SampleGroupLabels.Disambiguate(identities, labels);

            Assert.Multiple(() =>
            {
                Assert.That(result.Values, Is.EquivalentTo(new[] { "sample", "sample_2" }));
                Assert.That(result.Values, Has.None.StartsWith("_"), "no empty path segment is prefixed");
            });
        }
    }
}
