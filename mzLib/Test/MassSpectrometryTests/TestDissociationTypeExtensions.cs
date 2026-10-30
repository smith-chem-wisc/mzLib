using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.MassSpectrometryTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestDissociationTypeExtensions
    {
        /// <summary>
        /// Exhaustive rather than hand-picked: IsCid has to be true for exactly CID and LowCID over every
        /// member of the enum, so adding a member cannot quietly widen or narrow it.
        /// </summary>
        [Test]
        public static void IsCidIsTrueForExactlyCidAndLowCid()
        {
            var expected = new HashSet<DissociationType> { DissociationType.CID, DissociationType.LowCID };

            var actual = Enum.GetValues<DissociationType>().Where(d => d.IsCid()).ToHashSet();

            Assert.That(actual, Is.EquivalentTo(expected));
        }

        /// <summary>
        /// The types most likely to be swept in by a looser reading of "is CID". HCD is beam-type
        /// collision-induced dissociation and ISCID is in-source, but neither shares the b1 behaviour
        /// that the CID check guards.
        /// </summary>
        [Test]
        public static void IsCidExcludesTheOtherCollisionInducedTypes()
        {
            Assert.Multiple(() =>
            {
                Assert.That(DissociationType.HCD.IsCid(), Is.False);
                Assert.That(DissociationType.ISCID.IsCid(), Is.False);
                Assert.That(DissociationType.EThcD.IsCid(), Is.False);
                Assert.That(DissociationType.AnyActivationType.IsCid(), Is.False);
                Assert.That(DissociationType.Unknown.IsCid(), Is.False);
            });
        }

        /// <summary>
        /// Pins the equivalence the call site in PeptideWithSetModifications relied on before it was
        /// replaced, so the refactor stays behaviour-preserving for every enum member.
        /// </summary>
        [Test]
        public static void IsCidMatchesTheInlineComparisonItReplaced()
        {
            foreach (var d in Enum.GetValues<DissociationType>())
            {
                bool inline = d == DissociationType.CID || d == DissociationType.LowCID;
                Assert.That(d.IsCid(), Is.EqualTo(inline), d.ToString());
            }
        }
    }
}
