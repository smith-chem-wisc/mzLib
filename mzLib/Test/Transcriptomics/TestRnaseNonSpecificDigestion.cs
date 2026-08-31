using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Non-specific digestion of nucleic acids: every oligo anchored at one terminus. Before this,
    /// Rnase.GetUnmodifiedOligos threw for anything but Full and None, which is what stopped
    /// non-specific search from working on RNA databases at all.
    /// </summary>
    [ExcludeFromCodeCoverage]
    public class TestRnaseNonSpecificDigestion
    {
        private static List<Modification> NoMods => new List<Modification>();

        [Test]
        public void SingleRnasesAreInTheDictionary()
        {
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("singleN"));
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("singleC"));
            Assert.That(RnaseDictionary.Dictionary["singleN"].CleavageSpecificity, Is.EqualTo(CleavageSpecificity.SingleN));
            Assert.That(RnaseDictionary.Dictionary["singleC"].CleavageSpecificity, Is.EqualTo(CleavageSpecificity.SingleC));
        }

        /// <summary>
        /// SingleN anchors the 5' end: one oligo per start position, each running as far toward the 3'
        /// end as maxLength allows. With no length cap that is one oligo per start position, each
        /// ending at the 3' terminus.
        /// </summary>
        [Test]
        public void SingleN_ProducesOneOligoPerStartPosition()
        {
            var rna = new RNA("GUACUGA");   // length 7
            var digestionParams = new RnaDigestionParams("singleN", minLength: 3);

            var oligos = rna.Digest(digestionParams, NoMods, NoMods).Cast<OligoWithSetMods>().ToList();

            // starts 1..5 give lengths 7..3; starts 6 and 7 would be too short
            Assert.That(oligos.Select(o => o.OneBasedStartResidue), Is.EqualTo(new[] { 1, 2, 3, 4, 5 }));
            Assert.That(oligos.Select(o => o.OneBasedEndResidue).Distinct().Single(), Is.EqualTo(7),
                "every SingleN oligo should run to the 3' end when no max length is set");
            Assert.That(oligos.Select(o => o.BaseSequence),
                Is.EqualTo(new[] { "GUACUGA", "UACUGA", "ACUGA", "CUGA", "UGA" }));
        }

        /// <summary>
        /// SingleC is the mirror image: anchored at the 3' end, one oligo per end position.
        /// </summary>
        [Test]
        public void SingleC_ProducesOneOligoPerEndPosition()
        {
            var rna = new RNA("GUACUGA");
            var digestionParams = new RnaDigestionParams("singleC", minLength: 3);

            var oligos = rna.Digest(digestionParams, NoMods, NoMods).Cast<OligoWithSetMods>().ToList();

            Assert.That(oligos.Select(o => o.OneBasedEndResidue), Is.EqualTo(new[] { 3, 4, 5, 6, 7 }));
            Assert.That(oligos.Select(o => o.OneBasedStartResidue).Distinct().Single(), Is.EqualTo(1),
                "every SingleC oligo should run back to the 5' end when no max length is set");
            Assert.That(oligos.Select(o => o.BaseSequence),
                Is.EqualTo(new[] { "GUA", "GUAC", "GUACU", "GUACUG", "GUACUGA" }));
        }

        [Test]
        public void SingleN_RespectsMaxLength()
        {
            var rna = new RNA("GUACUGA");
            var digestionParams = new RnaDigestionParams("singleN", minLength: 3, maxLength: 4);

            var oligos = rna.Digest(digestionParams, NoMods, NoMods).Cast<OligoWithSetMods>().ToList();

            Assert.That(oligos.Select(o => o.BaseSequence),
                Is.EqualTo(new[] { "GUAC", "UACU", "ACUG", "CUGA", "UGA" }));
            Assert.That(oligos.All(o => o.BaseSequence.Length <= 4));
        }

        [Test]
        public void SingleC_RespectsMaxLength()
        {
            var rna = new RNA("GUACUGA");
            var digestionParams = new RnaDigestionParams("singleC", minLength: 3, maxLength: 4);

            var oligos = rna.Digest(digestionParams, NoMods, NoMods).Cast<OligoWithSetMods>().ToList();

            Assert.That(oligos.Select(o => o.BaseSequence),
                Is.EqualTo(new[] { "GUA", "GUAC", "UACU", "ACUG", "CUGA" }));
            Assert.That(oligos.All(o => o.BaseSequence.Length <= 4));
        }

        /// <summary>
        /// The termini are the nucleic-acid-specific part of this. An oligo that does not reach an
        /// original terminus gets the digested default: 5'-OH, 3'-phosphate. Getting this wrong would
        /// shift every mass without any test failing on sequence alone.
        /// </summary>
        [Test]
        public void SingleN_TerminiFollowTheDigestionRules()
        {
            var rna = new RNA("GUACUGA");
            var digestionParams = new RnaDigestionParams("singleN", minLength: 3, maxLength: 4);

            var oligos = rna.Digest(digestionParams, NoMods, NoMods).Cast<OligoWithSetMods>().ToList();

            var startsAtFivePrime = oligos.First(o => o.OneBasedStartResidue == 1);
            Assert.That(startsAtFivePrime.FivePrimeTerminus, Is.EqualTo(rna.FivePrimeTerminus),
                "an oligo starting at residue 1 keeps the original 5' terminus");

            var interior = oligos.First(o => o.OneBasedStartResidue == 2);
            Assert.That(interior.FivePrimeTerminus.MonoisotopicMass,
                Is.EqualTo(Rnase.DefaultFivePrimeTerminus.MonoisotopicMass).Within(1e-9),
                "an oligo cut at the 5' side gets the default 5'-OH");

            var endsAtThreePrime = oligos.First(o => o.OneBasedEndResidue == rna.Length);
            Assert.That(endsAtThreePrime.ThreePrimeTerminus, Is.EqualTo(rna.ThreePrimeTerminus),
                "an oligo ending at the last residue keeps the original 3' terminus");

            var cutAtThreePrime = oligos.First(o => o.OneBasedEndResidue != rna.Length);
            Assert.That(cutAtThreePrime.ThreePrimeTerminus.MonoisotopicMass,
                Is.EqualTo(Rnase.DefaultThreePrimeTerminus.MonoisotopicMass).Within(1e-9),
                "an oligo cut at the 3' side gets the default 3'-phosphate");
        }

        /// <summary>
        /// Every oligo a full digestion produces is also reachable non-specifically -- the point of a
        /// non-specific search being that it does not assume where the RNase cut.
        /// </summary>
        [Test]
        public void SingleN_CoversWhatFullDigestionFinds()
        {
            var rna = new RNA("GUACUGAAGGUCCAUG");

            var full = rna.Digest(new RnaDigestionParams("RNase T1", minLength: 3), NoMods, NoMods)
                .Cast<OligoWithSetMods>()
                .Select(o => o.OneBasedStartResidue)
                .Distinct();

            var singleN = rna.Digest(new RnaDigestionParams("singleN", minLength: 3), NoMods, NoMods)
                .Cast<OligoWithSetMods>()
                .Select(o => o.OneBasedStartResidue)
                .ToHashSet();

            Assert.That(full.All(start => singleN.Contains(start)),
                "a start position found by full digestion was not reachable by SingleN");
        }

        [Test]
        public void SingleN_AndSingleC_HonourMinLength()
        {
            var rna = new RNA("GUACUGA");

            var singleN = rna.Digest(new RnaDigestionParams("singleN", minLength: 5), NoMods, NoMods)
                .Cast<OligoWithSetMods>().ToList();
            var singleC = rna.Digest(new RnaDigestionParams("singleC", minLength: 5), NoMods, NoMods)
                .Cast<OligoWithSetMods>().ToList();

            Assert.That(singleN, Is.Not.Empty);
            Assert.That(singleC, Is.Not.Empty);
            Assert.That(singleN.All(o => o.BaseSequence.Length >= 5));
            Assert.That(singleC.All(o => o.BaseSequence.Length >= 5));
        }

    }
}
