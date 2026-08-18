using System.Collections.Generic;
using Omics;
using Omics.SpectralMatch;

namespace Development.QuantificationDevelopment.TestHelpers;

/// <summary>
/// Test implementation of ISpectralMatch for sequence coverage tests.
/// Also implements IHasSequenceCoverageFromFragments to support fragment coverage calculation.
/// </summary>
public class MockSpectralMatch : BaseSpectralMatch
{
    private readonly List<IBioPolymerWithSetMods> _identified;

    public void AddIdentifiedBioPolymer(IBioPolymerWithSetMods bioPolymer) => _identified.Add(bioPolymer);

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

    /// <summary>
    /// Constructor matching the BaseSpectralMatch parameter order for compatibility with quantification tests.
    /// </summary>
    public MockSpectralMatch(
        string fullFilePath,
        int oneBasedScanNumber,
        double score,
        string fullSequence,
        string baseSequence,
        IEnumerable<IBioPolymerWithSetMods>? identifiedBioPolymers = null) : base(fullFilePath, oneBasedScanNumber, score, fullSequence, baseSequence)
    {
        _identified = identifiedBioPolymers != null ? [.. identifiedBioPolymers] : [];
    }

    public override IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
}