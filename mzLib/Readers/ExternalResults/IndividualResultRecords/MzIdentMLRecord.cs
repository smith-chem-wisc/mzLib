using MzIdentML;
using Omics.Modifications;

namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// One SpectrumIdentificationItem read from an mzIdentML file by <see cref="ResultFiles.MzIdentMLResultFile"/>.
/// The file does not choose which items matter, so every item is a record: filter on <see cref="Rank"/> and
/// <see cref="PassThreshold"/> for what the submitter reported. There is no common score across search
/// engines, so each engine's own terms are in <see cref="Scores"/>.
/// </summary>
public class MzIdentMLRecord : ISpectralMatch
{
    /// <summary>The item as the document gives it, with its references already resolved.</summary>
    public MzidSpectrumMatch Match { get; internal set; } = null!;

    /// <summary>
    /// The scan number from the nativeID in <see cref="SpectrumId"/>. "scan=N" (Thermo, Proteome Discoverer)
    /// gives N. CAREFUL: "index=N", which peak-list files use, is the zero-based position of the spectrum in
    /// the file, not an instrument scan number, and gives N + 1. -1 when the nativeID has neither.
    /// </summary>
    public int OneBasedScanNumber { get; internal set; }

    /// <summary>The nativeID of the spectrum, as written.</summary>
    public string SpectrumId { get; internal set; } = string.Empty;

    /// <summary>The spectrum title, when the result carries one (typical for peak-list input), otherwise null.</summary>
    public string? SpectrumTitle { get; internal set; }

    /// <summary>
    /// The spectra file the identification came from, as named in its SpectraData entry. For Proteome
    /// Discoverer's "scan=N file=K" this is the RAW file K, not the mzML the result references.
    /// </summary>
    public string SpectraFileLocation { get; internal set; } = string.Empty;

    public string FileNameWithoutExtension { get; internal set; } = string.Empty;

    public string BaseSequence { get; internal set; } = string.Empty;

    public string FullSequence { get; internal set; } = string.Empty;

    public Dictionary<int, Modification> AllModsOneIsNterminus { get; internal set; } = [];

    /// <summary>The accessions of every protein the peptide evidence names, in document order, joined with '|'.</summary>
    public string Accession { get; internal set; } = string.Empty;

    /// <summary>True when every peptide evidence the item references is a decoy.</summary>
    public bool IsDecoy { get; internal set; }

    public int ChargeState { get; internal set; }

    public double ExperimentalMassToCharge { get; internal set; }

    public double? CalculatedMassToCharge { get; internal set; }

    public int Rank { get; internal set; }

    public bool PassThreshold { get; internal set; }

    /// <summary>
    /// The item's numeric cvParams and userParams, keyed by name (for example "Mascot:score" or
    /// "MS-GF:SpecEValue"). A name the item repeats keeps its first value.
    /// </summary>
    public IReadOnlyDictionary<string, double> Scores { get; internal set; } = new Dictionary<string, double>();

    /// <summary>
    /// The PSM-level q-value (MS:1002354, or one of its engine-specific children), or null when the item
    /// reports none.
    /// </summary>
    public double? QValue { get; internal set; }
}
