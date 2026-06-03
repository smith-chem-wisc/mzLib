using System;
using System.Collections.Generic;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class ModificationLookupBaseCachingTests
{
    [Test]
    public void TryResolve_StringCachesPositiveResult()
    {
        var modification = CreateModification("SerPhos", "S");
        var lookup = new CountingLookup(new[] { modification });

        var first = lookup.TryResolve("SerPhos on S", 'S');
        var second = lookup.TryResolve("SerPhos on S", 'S');

        Assert.That(first, Is.Not.Null);
        Assert.That(second, Is.Not.Null);
        Assert.That(lookup.FilterEvaluations, Is.EqualTo(1));
    }

    [Test]
    public void TryResolve_CanonicalRespectsResidueKey()
    {
        var mods = new[]
        {
            CreateModification("SerMod", "S"),
            CreateModification("ThrMod", "T")
        };

        var lookup = new CountingLookup(mods);
        var serCanonical = CanonicalModification.AtResidue(0, 'S', "SerMod on S");
        var thrCanonical = CanonicalModification.AtResidue(1, 'T', "ThrMod on T");

        var resolvedSer = lookup.TryResolve(serCanonical);
        var resolvedThr = lookup.TryResolve(thrCanonical);

        Assert.That(resolvedSer, Is.Not.Null);
        Assert.That(resolvedThr, Is.Not.Null);
        Assert.That(resolvedSer!.Value.MzLibModification.IdWithMotif, Does.Contain("SerMod"));
        Assert.That(resolvedThr!.Value.MzLibModification.IdWithMotif, Does.Contain("ThrMod"));
    }

    private static Modification CreateModification(string originalId, string residue)
    {
        if (!ModificationMotif.TryGetMotif(residue, out var motif))
        {
            throw new ArgumentException($"Invalid residue '{residue}' for motif creation.", nameof(residue));
        }

        return new Modification(
            _originalId: originalId,
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.");
    }

    private sealed class CountingLookup : ModificationLookupBase
    {
        public CountingLookup(IEnumerable<Modification> candidateSet)
            : base(candidateSet, massTolerance: 0.001)
        {
        }

        public override string Name => "CountingLookup";

        public int FilterEvaluations { get; private set; }

        protected override IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod) => Array.Empty<Modification>();

        protected override IEnumerable<Modification> FilterByName(IEnumerable<Modification> source, string name, char? targetResidue)
        {
            FilterEvaluations++;
            return base.FilterByName(source, name, targetResidue);
        }
    }
}
