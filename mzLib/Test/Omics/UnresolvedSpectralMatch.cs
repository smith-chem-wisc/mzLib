using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Omics;
using Omics.SpectralMatch;

namespace Test.Omics;

/// <summary>
/// A match whose base sequence could not be resolved. Implements ISpectralMatch directly rather than
/// deriving from BaseSpectralMatch, which coerces a null base sequence to empty - this is the shape
/// MetaMorpheus's SpectralMatch has when its best matches disagree, down to the reference equality.
/// </summary>
[ExcludeFromCodeCoverage]
public class UnresolvedSpectralMatch : ISpectralMatch
{
    public UnresolvedSpectralMatch(string filePath, int scanNumber, double score = 1.0)
    {
        FullFilePath = filePath;
        OneBasedScanNumber = scanNumber;
        Score = score;
    }

    public string FullFilePath { get; }
    public int OneBasedScanNumber { get; }
    public double Score { get; }
    public string BaseSequence => null!;
    public string FullSequence => null!;
    public string Accession => null!;
    public bool IsDecoy => false;
    public double[]? Intensities => null;

    public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => [];

    public int CompareTo(ISpectralMatch? other) => other is null ? 1 : other.Score.CompareTo(Score);

    public bool Equals(ISpectralMatch? other) => ReferenceEquals(this, other);
}
