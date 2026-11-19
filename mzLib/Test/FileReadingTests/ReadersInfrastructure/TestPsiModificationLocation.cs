using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="PsiModificationLocation"/> called directly, with the numbering mzIdentML and mzTab share: 0 is the
    /// N-terminus, 1 to length a residue, length + 1 the C-terminus. MzIdentMLResultFile resolves every modification
    /// through it, so TestMzIdentMLResultFile covers it in context; these cases pin the shared entry point itself.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestPsiModificationLocation
    {
        private static readonly UnimodModificationLookup Lookup = new();

        [TestCase(0, null)]
        [TestCase(1, null)]
        [TestCase(7, null)]
        [TestCase(8, null)]
        [TestCase(-1, "modification location -1 is not on the peptide")]
        [TestCase(9, "modification location 9 is not on the peptide")]
        public void ALocationIsOnThePeptideFromTheNTerminusToTheCTerminus(int location, string expectedReason)
        {
            Assert.That(PsiModificationLocation.ReasonNotOnPeptide("PEPTIDE", location), Is.EqualTo(expectedReason));
        }

        /// <summary>
        /// Every key is the location plus one. The residue a terminal modification sits on comes from the sequence.
        /// </summary>
        [Test]
        public void EachLocationIsKeyedOnePastIt()
        {
            var mods = new Dictionary<int, Modification>();

            Assert.Multiple(() =>
            {
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 0, "UNIMOD:1", 1, null, Lookup, out var nTerm), Is.True, nTerm);
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:21", 21, 79.966331, Lookup, out var residue), Is.True, residue);
                Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 8, "UNIMOD:2", 2, null, Lookup, out var cTerm), Is.True, cTerm);

                Assert.That(mods.Keys, Is.EquivalentTo(new[] { 1, 5, 9 }));
                Assert.That(mods[1].IdWithMotif, Is.EqualTo("Acetyl on X"));
                Assert.That(mods[5].IdWithMotif, Is.EqualTo("Phospho on T"));
                Assert.That(mods[9].IdWithMotif, Is.EqualTo("Amidated on X"));
                Assert.That(mods[9].LocationRestriction, Does.Contain("C-terminal"));
            });
        }

        /// <summary>
        /// Location length + 1 on a peptide ending in K is the C-terminal methyl, not the lysine one.
        /// </summary>
        [Test]
        public void TheLastLocationIsTheCTerminusEvenOnALysine()
        {
            var mods = new Dictionary<int, Modification>();

            Assert.That(PsiModificationLocation.TryAddUnimod(mods, "PEPTIDEK", 9, "UNIMOD:34", 34, null, Lookup, out _), Is.True);
            Assert.That(mods[10].IdWithMotif, Is.EqualTo("Methyl on X"));
            Assert.That(mods[10].LocationRestriction, Does.Contain("C-terminal"));
        }

        [TestCase(-1, "UNIMOD:35", 35, "modification location -1 is not on the peptide")]
        [TestCase(9, "UNIMOD:35", 35, "modification location 9 is not on the peptide")]
        [TestCase(2, "UNIMOD:999999", 999999, "UNIMOD:999999 at location 2 does not resolve to a Unimod modification")]
        public void AModificationThatCannotBePlacedIsNotAdded(int location, string accession, int unimodId, string expectedReason)
        {
            var mods = new Dictionary<int, Modification>();

            bool added = PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", location, accession, unimodId, null, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Is.EqualTo(expectedReason));
            Assert.That(mods, Is.Empty);
        }

        [Test]
        public void ASecondModificationAtOneLocationIsNotAddedAndTheFirstIsKept()
        {
            var mods = new Dictionary<int, Modification>();
            PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:21", 21, 79.966331, Lookup, out _);

            bool added = PsiModificationLocation.TryAddUnimod(mods, "PEPTIDE", 4, "UNIMOD:1", 1, null, Lookup, out var reason);

            Assert.That(added, Is.False);
            Assert.That(reason, Is.EqualTo("more than one modification at location 4"));
            Assert.That(mods[5].IdWithMotif, Is.EqualTo("Phospho on T"));
        }
    }
}
