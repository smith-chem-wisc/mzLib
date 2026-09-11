using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Pins the surface a derived group outside this repository binds to.
    ///
    /// MetaMorpheus's ProteinGroup overrides <see cref="BioPolymerGroup.GetTabSeparatedHeader"/> and
    /// reads the protected quantification gate. Removing either member is not a test failure here --
    /// it is a COMPILE failure in a different repository, which mzLib's own CI cannot see and which
    /// the integration job only catches if its clone of MetaMorpheus master happens to be current.
    ///
    /// <see cref="DerivedGroup"/> below stands in for that consumer. If the virtual header, the
    /// interface member, or the protected gate is removed, THIS FILE STOPS COMPILING -- so the break
    /// surfaces in mzLib's own build rather than in a downstream bump.
    ///
    /// Delete this file in the same PR that removes the members, once every implementation has
    /// switched to the schema and the writer.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class DerivedGroupRenderingTests
    {
        /// <summary>
        /// A group that renders itself the way an out-of-repository subclass does: overriding the
        /// header and consulting the protected gate to decide whether it has intensity columns.
        /// </summary>
        private sealed class DerivedGroup : BioPolymerGroup
        {
            public DerivedGroup(HashSet<IBioPolymer> bioPolymers,
                HashSet<IBioPolymerWithSetMods> all, HashSet<IBioPolymerWithSetMods> unique)
                : base(bioPolymers, all, unique)
            {
            }

            public override string GetTabSeparatedHeader() => "Derived\t" + base.GetTabSeparatedHeader();

            /// <summary>Reads the protected member, as ProteinGroup does at two sites.</summary>
            public bool SeesTheQuantificationGate => HasAssignedSampleIntensities;
        }

        private static DerivedGroup NewGroup()
        {
            var bioPolymer = new MockBioPolymer("ACGTACGT", "BP12345");
            var sequence = new MockBioPolymerWithSetMods("ACGT", "ACGT");
            return new DerivedGroup(
                new HashSet<IBioPolymer> { bioPolymer },
                new HashSet<IBioPolymerWithSetMods> { sequence },
                new HashSet<IBioPolymerWithSetMods> { sequence });
        }

        [Test]
        public void DerivedGroupCanOverrideTheHeaderAndReachTheBaseImplementation()
        {
            var group = NewGroup();

            Assert.That(group.GetTabSeparatedHeader(), Does.StartWith("Derived\t"));
            Assert.That(group.GetTabSeparatedHeader(), Does.Contain("BioPolymer Accession"),
                "the base implementation must still produce the group's own columns");
        }

        [Test]
        public void TheInterfaceStillCarriesTheHeaderMember()
        {
            IBioPolymerGroup group = NewGroup();

            // Through the interface, not the class -- this is the binding an implementation outside
            // this repository satisfies.
            Assert.That(group.GetTabSeparatedHeader(), Is.Not.Empty);
        }

        [Test]
        public void TheQuantificationGateIsFalseUntilBothHalvesAreAssigned()
        {
            var group = NewGroup();
            Assert.That(group.SeesTheQuantificationGate, Is.False, "neither half assigned");

            group.SamplesForQuantification = new List<ISampleInfo>();
            Assert.That(group.SeesTheQuantificationGate, Is.False, "an empty sample list is not quantification");

            group.IntensitiesBySample = new Dictionary<ISampleInfo, double>();
            Assert.That(group.SeesTheQuantificationGate, Is.False,
                "intensities without samples have nothing to label");
        }

        /// <summary>
        /// The point of keeping the members rather than reimplementing them: the header a group
        /// renders for itself and the header the schema builds for that same group are the same
        /// string, because one delegates to the other. Two implementations that could drift is the
        /// defect this whole refactor exists to remove.
        /// </summary>
        [Test]
        public void TheGroupsOwnHeaderIsTheSchemasHeaderForThatGroup()
        {
            var group = NewGroup();
            BioPolymerGroup asBase = group;

            Assert.That(asBase.ToString(),
                Is.EqualTo(TsvWriter.RowLine(BioPolymerGroupTsvSchema.For(new[] { asBase }), asBase)));
            Assert.That(group.GetTabSeparatedHeader(),
                Is.EqualTo("Derived\t" + TsvWriter.HeaderLine(BioPolymerGroupTsvSchema.For(new[] { asBase }))));
        }
    }
}
