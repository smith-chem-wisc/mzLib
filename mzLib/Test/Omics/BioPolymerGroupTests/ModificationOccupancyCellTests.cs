using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymerGroup;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// <see cref="ModificationOccupancyCell"/> is the inverse of
    /// <see cref="SiteSpecificModificationOccupancy.ToModInfoString"/>; these tests hold it to that.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class ModificationOccupancyCellTests
    {
        [TestCase("Oxidation on M")]
        [TestCase("N6,N6,N6-trimethyllysine on K")]
        [TestCase("HexNAc(1)Hex(1)[glycan] on S")]
        public void RoundTripsTheCountFormatter(string name)
        {
            var site = new SiteSpecificModificationOccupancy(53, name) { ModifiedCount = 3, TotalCount = 7 };
            var parsed = ModificationOccupancyCell.Parse(site.ToModInfoString()).Sites.Single();
            Assert.That(parsed, Is.EqualTo(new OccupancySite(52, name, 0.43, 3, 7)));
        }

        [Test]
        public void RoundTripsTheIntensityFormatterToItsPrintedPrecision()
        {
            var site = new SiteSpecificModificationOccupancy(1, "Acetylation on X")
            { ModifiedIntensity = 1234567.0, TotalIntensity = 98765432.0 };
            var parsed = ModificationOccupancyCell.Parse(site.ToModInfoString(intensityBased: true)).Sites.Single();
            Assert.That(parsed.IsNTerminus);
            Assert.That(parsed.Fraction, Is.EqualTo(0.0125));
            Assert.That((parsed.Numerator, parsed.Denominator), Is.EqualTo((1.235e6, 9.877e7)), "G4: four significant digits");
        }

        [Test]
        public void SemicolonsJoinSitesAndPipesJoinEntities()
        {
            var cell = ModificationOccupancyCell.Parse(
                "pos1[A on K,info:fraction=1.00(1/1)];pos5[B on S,info:fraction=0.50(1/2)]|pos9[C on T,info:fraction=0.25(1/4)]");
            Assert.That(cell.Entities.Select(e => e.Count), Is.EqualTo(new[] { 2, 1 }));
            Assert.That(cell.IsTruncated, Is.False);
        }

        [TestCase(null)]
        [TestCase("")]
        [TestCase("   ")]
        public void BlankIsNoSites(string? text)
        {
            var cell = ModificationOccupancyCell.Parse(text);
            Assert.That(cell.Sites, Is.Empty);
            Assert.That(cell.IsTruncated, Is.False);
        }

        /// <summary>MetaMorpheus 1.1.x replaces an over-long cell wholesale. That is not an empty cell.</summary>
        [Test]
        public void TheExcelPlaceholderIsTruncatedNotEmpty()
        {
            var cell = ModificationOccupancyCell.Parse(ModificationOccupancyCell.ExcelTruncationText);
            Assert.That(cell.IsTruncated);
            Assert.That(cell.Sites, Is.Empty);
        }

        /// <summary>mzLib's writer cuts at <see cref="BioPolymerGroupTsvSchema.MaxStringLength"/>, mid-site.
        /// The complete sites before the cut are kept and the cell says it was cut.</summary>
        [Test]
        public void ACutAtTheLengthLimitKeepsTheCompleteSites()
        {
            string one = "pos12[Phosphorylation on S,info:fraction=0.50(1/2)]";
            string full = string.Join(";", Enumerable.Repeat(one, BioPolymerGroupTsvSchema.MaxStringLength / one.Length + 2));
            var cell = ModificationOccupancyCell.Parse(full[..BioPolymerGroupTsvSchema.MaxStringLength]);
            Assert.That(cell.IsTruncated);
            Assert.That(cell.Sites.Count(), Is.EqualTo(BioPolymerGroupTsvSchema.MaxStringLength / (one.Length + 1)));
        }

        /// <summary>The limit is a writer setting the reader cannot know, so a cut is recognised at any
        /// length and at any character of the unfinished site, including inside a bracketed name.</summary>
        [Test]
        public void ACutIsRecognisedWhateverTheLimitWas()
        {
            string done = "pos1[A on K,info:fraction=1.00(1/1)];";
            string last = "pos7[HexNAc(1)Hex(1)[glycan] on S,info:fraction=0.50(1/2)]";
            for (int cut = 1; cut < last.Length; cut++)
            {
                if (last[..cut].EndsWith(']'))
                    continue; // a cut straight after a "]" in the name cannot be told from a finished site
                var cell = ModificationOccupancyCell.Parse(done + last[..cut]);
                Assert.That(cell.IsTruncated, $"cut after {cut} characters");
                Assert.That(cell.Sites.Single().ModificationIdWithMotif, Is.EqualTo("A on K"));
            }
            Assert.That(ModificationOccupancyCell.Parse("pos12[Phospho").IsTruncated, "a cut inside the first site");
        }

        [TestCase("not an occupancy cell")]
        [TestCase("pos3[Oxidation on M,info:occupancy=0.50(1/2)]")]
        [TestCase("pos1[A on K,info:fraction=1.00(1/1)];pos3[Oxidation on M,info:occupancy=0.50(1/2)]")]
        [TestCase("pos1[A on K,info:fraction=1.00(1/1)];pos3[Oxidation on M,info:fraction=0.50(1/2)]x")]
        [TestCase("pos1[A on K,info:fraction=1.00(1/1)];pos3[B on M];pos5[C on S,info:fraction=1.00(1/1)]")]
        [TestCase("pos1[A on K,info:fraction=1.00(1/1)]pos2[B on K,info:fraction=1.00(1/1)]")]
        [TestCase("pos1[A on K,info:fraction=1.00(1/1)];pos2[B on K,info:fraction=1.00(1/1)]pos3[C on K,info:fraction=1.00(1/1)]")]
        public void AnythingElseIsRefusedNotGuessed(string text)
        {
            Assert.Throws<FormatException>(() => ModificationOccupancyCell.Parse(text));
        }
    }
}
