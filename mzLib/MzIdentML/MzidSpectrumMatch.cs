using MzLibUtil;
using System.Collections.Generic;

namespace MzIdentML
{
    /// <summary>
    /// One SpectrumIdentificationItem with everything a reader needs from the rest of the document already
    /// resolved: its result, the spectra file the result points at, the peptide and its modifications, and
    /// every PeptideEvidence the item references. Returned by <see cref="MzidIdentifications.GetSpectrumMatches"/>,
    /// and the same shape for every mzIdentML schema version.
    /// </summary>
    /// <remarks>
    /// A reference the document does not resolve is reported as null (a SpectraData or Peptide) or left out
    /// (a PeptideEvidence), rather than failing the whole file.
    /// </remarks>
    public sealed record MzidSpectrumMatch
    {
        public string SpectrumIdentificationListId { get; init; }
        public string SpectrumIdentificationResultId { get; init; }
        public string SpectrumIdentificationItemId { get; init; }

        /// <summary>The result's spectrumID, the nativeID of the spectrum in <see cref="SpectraData"/>.</summary>
        public string SpectrumId { get; init; }

        /// <summary>The SpectraData the result references, or null if the reference does not resolve.</summary>
        public MzidSpectraData SpectraData { get; init; }

        /// <summary>
        /// The result's spectrum title (MS:1000796, or the obsolete MS:1001416), or null when the result carries
        /// none or it has no value.
        /// </summary>
        public string SpectrumTitle { get; init; }

        public int Rank { get; init; }
        public bool PassThreshold { get; init; }
        public int ChargeState { get; init; }
        public double ExperimentalMassToCharge { get; init; }

        /// <summary>Null when the item does not state it; the attribute is optional.</summary>
        public double? CalculatedMassToCharge { get; init; }

        /// <summary>The peptide's sequence, or null if the item's peptide reference does not resolve.</summary>
        public string PeptideSequence { get; init; }

        public IReadOnlyList<MzidModification> Modifications { get; init; } = [];

        /// <summary>Whether the peptide carries any SubstitutionModification (a residue swapped for another).</summary>
        public bool HasSubstitutionModifications { get; init; }

        /// <summary>The PeptideEvidence entries the item references, in document order.</summary>
        public IReadOnlyList<MzidPeptideEvidence> PeptideEvidence { get; init; } = [];

        public IReadOnlyList<CvParam> ResultCvParams { get; init; } = [];
        public IReadOnlyList<MzidUserParam> ResultUserParams { get; init; } = [];
        public IReadOnlyList<CvParam> ItemCvParams { get; init; } = [];
        public IReadOnlyList<MzidUserParam> ItemUserParams { get; init; } = [];
    }

    /// <summary>A SpectraData entry: one spectra file the identifications were made against.</summary>
    public sealed record MzidSpectraData(string Id, string Name, string Location, CvParam FileFormat, CvParam SpectrumIdFormat);

    /// <summary>
    /// A peptide modification. <see cref="Location"/> is 0 for the N-terminus, 1 to the sequence length for a
    /// residue, and the length plus one for the C-terminus; it is null when the document omits it.
    /// </summary>
    public sealed record MzidModification
    {
        public int? Location { get; init; }
        public IReadOnlyList<string> Residues { get; init; } = [];
        public double? MonoisotopicMassDelta { get; init; }
        public IReadOnlyList<CvParam> CvParams { get; init; } = [];
    }

    /// <summary>A PeptideEvidence entry, with its DBSequence's accession resolved (null if it does not resolve).</summary>
    public sealed record MzidPeptideEvidence(string Id, bool IsDecoy, string DBSequenceAccession, int? Start, int? End, string Pre, string Post);

    /// <summary>A userParam: a name/value pair outside any controlled vocabulary.</summary>
    public sealed record MzidUserParam(string Name, string Value, string Type);
}
