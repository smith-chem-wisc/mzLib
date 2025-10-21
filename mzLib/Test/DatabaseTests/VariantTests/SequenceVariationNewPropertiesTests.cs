using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    public class SequenceVariationNewPropertiesTests
    {
        private static Modification DummyMod(string id = "Mod1") => new Modification(_originalId: id);

        [Test]
        public void SearchableAnnotation_PrefersVcfLine()
        {
            string vcf = "1\t100\t.\tA\tT\t.\tPASS\tANN=T|missense_variant|MODERATE|GENE1|GENE1|transcript|TX1|protein_coding|1/1|c.100A>T|p.Lys34Asn|100/1000|34/300|34/100|0|\tGT:AD:DP\t0/1:5,4:9\t1/1:0,10:10";
            var sv = new SequenceVariation(10, 10, "A", "T", "free", vcf);
            Assert.That(sv.SearchableAnnotation, Is.EqualTo(vcf));
        }

        [Test]
        public void SearchableAnnotation_FallsBackToDescription()
        {
            var sv = new SequenceVariation(5, 5, "K", "R", "myDesc");
            Assert.That(sv.SearchableAnnotation, Is.EqualTo("myDesc"));
        }

        [Test]
        public void AllelePassthrough_Reference_Alternate()
        {
            string vcf = "1\t100\t.\tA\tT\t.\tPASS\tANN=T|missense_variant|MODERATE|G|G|transcript|TX|protein_coding|1/1|c.100A>T|p.K34N|100/1000|34/300|34/100|0|\tGT:AD:DP\t0/1:5,4:9";
            var sv = new SequenceVariation(10, 10, "A", "T", "desc", vcf);
            Assert.Multiple(() =>
            {
                Assert.That(sv.ReferenceAllele, Is.EqualTo("A"));
                Assert.That(sv.AlternateAllele, Is.EqualTo("T"));
            });
        }

        [Test]
        public void ClassificationPredicates_Work()
        {
            var point = new SequenceVariation(1, 1, "A", "V", "point");
            Assert.Multiple(() =>
            {
                Assert.That(point.IsPointSubstitution, Is.True);
                Assert.That(point.IsMultiResidueSubstitution, Is.False);
                Assert.That(point.IsInsertion, Is.False);
                Assert.That(point.IsDeletion, Is.False);
                Assert.That(point.IsStopGain, Is.False);
                Assert.That(point.IsLikelyFrameshift, Is.False);
            });

            var multi = new SequenceVariation(2, 3, "AA", "VV", "multi");
            Assert.That(multi.IsMultiResidueSubstitution, Is.True);

            var insertion = new SequenceVariation(5, null, "M", "ins");
            Assert.That(insertion.IsInsertion, Is.True);

            var deletion = new SequenceVariation(7, 9, "ABC", "", "del");
            Assert.That(deletion.IsDeletion, Is.True);
                
            var stop = new SequenceVariation(4, 4, "Q", "W*", "stop");
            Assert.That(stop.IsStopGain, Is.True);

            var frameshift = new SequenceVariation(10, 12, "ABC", "AB", "fs");
            Assert.That(frameshift.IsLikelyFrameshift, Is.True);
        }

        [Test]
        public void PointSubstitution_FalseWhenNoChange()
        {
            Assert.That(() => new SequenceVariation(3, 3, "A", "A", "noop"),
                Throws.TypeOf<ArgumentException>());

            // identical but with a variant-specific mod is allowed
            var mods = new Dictionary<int, List<Modification>> { { 3, new List<Modification> { DummyMod() } } };
            var sv = new SequenceVariation(
                3,
                3,
                "A",
                "A",
                "noopWithMod",
                variantCallFormatDataString: null,          // disambiguate (string? overload)
                oneBasedModifications: mods);
            Assert.Multiple(() =>
            {
                Assert.That(sv.IsPointSubstitution, Is.False);
                Assert.That(sv.AreValid(), Is.True);
            });
        }

        [Test]
        public void InvalidModificationPositions_Throw()
        {
            var badMods = new Dictionary<int, List<Modification>> { { 6, new List<Modification> { DummyMod() } } };
            Assert.That(() => new SequenceVariation(
                    5,
                    7,
                    "ABC",
                    "A",
                    "shrink",
                    variantCallFormatDataString: null,   // disambiguate overload (string? param)
                    oneBasedModifications: badMods),
                Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void DeletionModificationInvalid()
        {
            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { DummyMod() } } };
            Assert.That(() => new SequenceVariation(
                    5,
                    7,
                    "ABC",
                    "",
                    "del",
                    variantCallFormatDataString: null,   // disambiguate overload
                    oneBasedModifications: mods),
                Throws.TypeOf<ArgumentException>());
        }
            
        [Test]
        public void SplitPerGenotype_ProducesExpectedVariants()
        {
            string vcf =
                "1\t100\t.\tA\tT\t.\tPASS\t" +
                "ANN=T|missense_variant|MODERATE|GENE1|GENE1|transcript|TX1|protein_coding|1/1|c.100A>T|p.K34N|100/1000|34/300|34/100|0|\t" +
                "GT:AD:DP\t0/1:5,4:9\t1/1:0,10:10";

            var sv = new SequenceVariation(34, 34, "K", "N", "origDesc", vcf);
            var perSample = sv.SplitPerGenotype(includeReferenceForHeterozygous: true);

            // Rationale:
            // The constructor/validation forbids no?op variants (ref->ref with no variant-specific mods).
            // The heterozygous reference copy therefore cannot be materialized and is skipped.
            // Expected:
            //   Sample 0: HeterozygousAlt
            //   Sample 1: HomozygousAlt
            // Total: 2 variants (both with VCF metadata)
            Assert.Multiple(() =>
            {
                Assert.That(perSample, Has.Count.EqualTo(2));
                Assert.That(perSample.Count(v => v.Description.Contains("Sample=0")), Is.EqualTo(1));
                Assert.That(perSample.Count(v => v.Description.Contains("Sample=1")), Is.EqualTo(1));
                Assert.That(perSample.All(v => v.VariantCallFormatData != null), Is.True);
                Assert.That(perSample.Any(v => v.Description.Contains("HeterozygousAlt")), Is.True);
                Assert.That(perSample.Any(v => v.Description.Contains("HomozygousAlt")), Is.True);
                Assert.That(perSample.Any(v => v.Description.Contains("HeterozygousRef")), Is.False);
            });
        }

        [Test]
        public void CombineEquivalent_MergesDescriptionsAndMods()
        {
            var a1 = new SequenceVariation(10, 11, "AA", "VV", "desc1");
            var a2 = new SequenceVariation(
                10,
                11,
                "AA",
                "VV",
                "desc2",
                variantCallFormatDataString: null,          // disambiguate
                oneBasedModifications: new Dictionary<int, List<Modification>> {
                    { 11, new List<Modification>{ DummyMod("M1") } }
                });

            var combined = SequenceVariation.CombineEquivalent(new[] { a1, a2 });
            Assert.That(combined, Has.Count.EqualTo(1));

            var merged = combined[0];
            Assert.Multiple(() =>
            {
                Assert.That(merged.Description, Does.StartWith("Combined(2):"));
                Assert.That(merged.OneBasedModifications, Has.Count.EqualTo(1));
                Assert.That(merged.OneBasedModifications.ContainsKey(11), Is.True);
            });
        }

        [Test]
        public void Equality_IgnoresDescriptionButRequiresCoreData()
        {
            var v1 = new SequenceVariation(5, 5, "A", "V", "d1");
            var v2 = new SequenceVariation(5, 5, "A", "V", "d2");
            var v3 = new SequenceVariation(5, 5, "A", "I", "d3");

            Assert.Multiple(() =>
            {
                Assert.That(v1.Equals(v2), Is.True);
                Assert.That(v1.Equals(v3), Is.False);
            });
        }

        [Test]
        public void ConvenienceCtor_SetsEndCoordinate()
        {
            var sv = new SequenceVariation(10, "ABC", "XYZ", "multi");
            Assert.Multiple(() =>
            {
                Assert.That(sv.OneBasedBeginPosition, Is.EqualTo(10));
                Assert.That(sv.OneBasedEndPosition, Is.EqualTo(12));
            });
        }

        [Test]
        public void SimpleString_PointAndSpanFormats()
        {
            var point = new SequenceVariation(4, 4, "A", "V", "p");
            var span = new SequenceVariation(10, 12, "ABC", "ADE", "s");

            Assert.Multiple(() =>
            {
                Assert.That(point.SimpleString(), Is.EqualTo("A4V"));
                Assert.That(span.SimpleString(), Is.EqualTo("ABC10-12ADE"));
            });
        }

        [Test]
        public void LegacyVariantDescription_ReturnsUnderlying()
        {
            string vcf = "1\t200\t.\tG\tC\t.\tPASS\tANN=C|missense_variant|LOW|G|G|transcript|TX|protein_coding|1/1|c.200G>C|p.G67A|200/900|67/300|67/100|0|\tGT:AD:DP\t0/1:3,6:9";
            var sv = new SequenceVariation(67, 67, "G", "A", "desc", vcf);
            Assert.That(sv.LegacyVariantDescription, Is.SameAs(sv.VariantCallFormatData));
        }

        [Test]
        public void StopGain_NotFrameshift()
        {
            var stop = new SequenceVariation(20, 22, "QWE", "QW*", "stop");
            Assert.Multiple(() =>
            {
                Assert.That(stop.IsStopGain, Is.True);
                Assert.That(stop.IsLikelyFrameshift, Is.False);
            });
        }

        [Test]
        public void Frameshift_NoInsertionDeletionOrStop()
        {
            var fs = new SequenceVariation(50, 52, "ABC", "AB", "fs");
            Assert.That(fs.IsLikelyFrameshift, Is.True);
        }
    }
}