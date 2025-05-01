using Omics.Modifications;

namespace Readers;

public interface ISpectralMatch
{
    /// <summary>
    /// The scan number of the identification
    /// </summary>
    public int OneBasedScanNumber { get; }

    /// <summary>
    /// Primary Sequence
    /// </summary>
    public string BaseSequence { get; }

    /// <summary>
    /// Modified Sequence in MetaMorpheus format
    /// </summary>
    public string FullSequence { get; }

    /// <summary>
    /// The accession (unique identifier) of the identification
    /// </summary>
    public string Accession { get; }

    /// <summary>
    /// If the given Spectral Match is a decoy
    /// </summary>
    public bool IsDecoy { get; }

    /// <summary>
    /// The Mass Spec file name without the extension
    /// </summary>
    public string FileNameWithoutExtension { get; }

    /// <summary>
    /// Modifications on the spectral match
    /// </summary>
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
}