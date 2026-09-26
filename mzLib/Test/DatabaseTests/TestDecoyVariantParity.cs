using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    /// <summary>
    /// A target-decoy database has to hold one decoy per target, and a decoy has to be a reversal of the
    /// target it stands for. Both properties are stated here as one oracle: reversing a variant-applied
    /// protein must give the same sequence as applying the reversed variation to the reversed protein.
    /// </summary>
    [TestFixture]
    public static class TestDecoyVariantParity
    {
        private const string StartsWithMethionine = "MPEPTIDEKR";

        private static Protein Consensus(int begin, int end, string original, string variant) =>
            new Protein(StartsWithMethionine, "ACC", sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(begin, end, original, variant, "info")
            });

        private static string ReverseOf(Protein protein) =>
            DecoyProteinGenerator.GenerateDecoys(new List<Protein> { protein }, DecoyType.Reverse).Single().BaseSequence;

        /// <summary>
        /// Reversal must commute with variant application. Ordinary substitutions and indels already
        /// satisfied this; a variation that only inserts residues after the initiator methionine did not.
        /// </summary>
        [Test]
        [TestCase(5, 5, "T", "A", TestName = "SubstitutionInTheMiddle")]
        [TestCase(10, 10, "R", "A", TestName = "SubstitutionAtTheCTerminus")]
        [TestCase(5, 5, "T", "TGG", TestName = "InsertionInTheMiddle")]
        [TestCase(4, 6, "PTI", "P", TestName = "DeletionInTheMiddle")]
        [TestCase(3, 5, "EPT", "QQQ", TestName = "MultiResidueSubstitution")]
        [TestCase(1, 1, "M", "MKAV", TestName = "InsertionAfterTheInitiatorMethionine")]
        public static void ReversalCommutesWithVariantApplication(int begin, int end, string original, string variant)
        {
            var consensus = Consensus(begin, end, original, variant);

            var appliedTarget = consensus.GetVariantBioPolymers(4, 1)
                .Single(p => p.AppliedSequenceVariations.Any());
            string reverseOfApplied = ReverseOf((Protein)appliedTarget);

            var decoyConsensus = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { consensus }, DecoyType.Reverse).Single();
            var appliedDecoy = decoyConsensus.GetVariantBioPolymers(4, 1)
                .Single(p => p.AppliedSequenceVariations.Any());

            Assert.That(decoyConsensus.SequenceVariations.Single().AreValid(), Is.True,
                "a reversed variation that reads as invalid is discarded downstream");
            Assert.That(appliedDecoy.BaseSequence, Is.EqualTo(reverseOfApplied));
        }

        /// <summary>
        /// The parity failure this fixes. Reversing an insertion after the initiator methionine produced a
        /// variation beginning one past the end of the sequence, with an empty original sequence, which reads
        /// as invalid - and one invalid variation made every variation on that decoy be discarded, so a
        /// protein with two variants yielded four targets and a single decoy.
        /// </summary>
        [Test]
        public static void DecoyOfAProteinWithAnInsertionAfterMethionineKeepsEveryVariant()
        {
            var consensus = new Protein(StartsWithMethionine, "ACC", sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(1, 1, "M", "MKAV", "insertion after the initiator methionine"),
                new SequenceVariation(5, 5, "T", "A", "ordinary substitution")
            });

            int targets = consensus.GetVariantBioPolymers(4, 1).Count;
            int decoys = DecoyProteinGenerator.GenerateDecoys(new List<Protein> { consensus }, DecoyType.Reverse)
                .Single().GetVariantBioPolymers(4, 1).Count;

            Assert.That(targets, Is.EqualTo(4), "base, each variant alone, and both together");
            Assert.That(decoys, Is.EqualTo(targets));
        }

        /// <summary>
        /// A variation that cannot be expressed as one contiguous replacement on the reversed sequence must
        /// cost only itself. A start loss is such a case: it drops the initiator methionine and reverses the
        /// rest, changing both ends at once. Before, it took every other variation on the protein with it.
        /// </summary>
        [Test]
        public static void AnUnrepresentableVariationDoesNotDiscardTheOthers()
        {
            var consensus = new Protein(StartsWithMethionine, "ACC", sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(1, 1, "M", "A", "start loss"),
                new SequenceVariation(5, 5, "T", "A", "ordinary substitution")
            });

            var decoyVariations = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { consensus }, DecoyType.Reverse).Single().SequenceVariations;

            Assert.That(decoyVariations.Count, Is.EqualTo(1), "the start loss cannot be reversed, the substitution can");
            Assert.That(decoyVariations.Single().AreValid(), Is.True);
            Assert.That(decoyVariations.All(v => v.AreValid()), Is.True,
                "an invalid variation makes GetVariantBioPolymers discard all of them");
        }

        /// <summary>
        /// An accession has to identify one sequence, because it is the key everything downstream joins on -
        /// Spritz writes accession-to-sequence and accession-to-variant pivot tables straight from these.
        ///
        /// Two variations beginning at the same position are the case that broke it. Overlapping edits are
        /// spliced together rather than replacing one another, but the record of what was applied dropped any
        /// earlier variation whose span the later one covered, and a deletion's post-application span is the
        /// single residue it leaves behind. So the deletion vanished from the accession while staying in the
        /// sequence. Reversal maps variations sharing an end onto variations sharing a begin, which is why
        /// decoys hit this far more often than targets.
        /// </summary>
        [Test]
        public static void OverlappingVariationsAtOneBeginPositionGiveDistinctAccessions()
        {
            var consensus = new Protein("MABCDEFGHIJ", "ACC", sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(2, 6, "ABCDE", "X", "deletion collapsing to one residue"),
                new SequenceVariation(2, 4, "ABC", "YY", "substitution starting at the same position")
            });

            var entries = consensus.GetVariantBioPolymers(4, 1).ToList();

            var withBoth = entries.Where(p => p.AppliedSequenceVariations.Count == 2).ToList();
            Assert.That(withBoth, Is.Not.Empty,
                "premise: both variations are applied together, and both edits stay in the sequence");

            var accessions = entries.Select(p => p.Accession).ToList();
            Assert.That(accessions, Is.Unique);

            foreach (var group in entries.GroupBy(p => p.Accession))
            {
                Assert.That(group.Select(p => p.BaseSequence).Distinct().Count(), Is.EqualTo(1),
                    $"accession {group.Key} must not name two different sequences");
            }
        }

        /// <summary>
        /// The same property across a target-decoy pair, which is where it actually bit.
        /// </summary>
        [Test]
        public static void ADecoyDatabaseNamesEachSequenceOnce()
        {
            var consensus = new Protein("MABCDEFGHIJKLMNOP", "ACC", sequenceVariations: new List<SequenceVariation>
            {
                // sharing an end, so reversal makes them share a begin
                new SequenceVariation(4, 17, "DEFGHIJKLMNOP", "Q", "long deletion to the C terminus"),
                new SequenceVariation(12, 17, "LMNOP", "RR", "shorter edit ending at the same place")
            });

            var decoyConsensus = DecoyProteinGenerator.GenerateDecoys(
                new List<Protein> { consensus }, DecoyType.Reverse).Single();

            Assert.That(decoyConsensus.SequenceVariations.Select(v => v.OneBasedBeginPosition).Distinct().Count(),
                Is.EqualTo(1), "premise: reversal put both variations at the same begin position");

            var accessions = consensus.GetVariantBioPolymers(4, 1).Select(p => p.Accession)
                .Concat(decoyConsensus.GetVariantBioPolymers(4, 1).Select(p => p.Accession)).ToList();

            Assert.That(accessions, Is.Unique);
        }

        /// <summary>
        /// An applied stop gain truncates the sequence but keeps recording its position in the untruncated
        /// one, so reversing such a protein indexed past the end of the sequence and threw. Reachable from any
        /// caller that generates decoys from variant-applied proteins, which is a public entry point.
        /// </summary>
        [Test]
        public static void ReversingAnAppliedStopGainDoesNotThrow()
        {
            var consensus = Consensus(5, 5, "T", "*");
            var appliedTarget = (Protein)consensus.GetVariantBioPolymers(4, 1)
                .Single(p => p.AppliedSequenceVariations.Any());

            Assert.That(appliedTarget.BaseSequence.Length,
                Is.LessThan(appliedTarget.AppliedSequenceVariations.Single().OneBasedEndPosition),
                "premise: the recorded position is past the end of the truncated sequence");
            Assert.That(() => ReverseOf(appliedTarget), Throws.Nothing);
        }
    }
}
