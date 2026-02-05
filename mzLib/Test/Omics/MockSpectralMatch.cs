using System.Collections.Generic;
using Omics;
using Omics.SpectralMatch;

namespace Test.Omics;

/// <summary>
/// Test implementation of ISpectralMatch for sequence coverage tests.
/// Also implements IHasSequenceCoverageFromFragments to support fragment coverage calculation.
/// </summary>
public class MockSpectralMatch : BaseSpectralMatch
{
    private readonly List<IBioPolymerWithSetMods> _identified;
    public MockSpectralMatch(
        string filePath,
        string fullSequence,
        string baseSequence,
        double score,
        int scanNumber,
        IEnumerable<IBioPolymerWithSetMods>? identified = null) : base(filePath, scanNumber, score, fullSequence, baseSequence)
    {
        _identified = identified != null ? [.. identified] : [];
    }

    public override IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
}