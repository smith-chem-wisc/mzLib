namespace Readers.ExternalResults.IndividualResultRecords;

/// <summary>
/// A SpectrumIdentificationItem that <see cref="ResultFiles.MzIdentMLResultFile"/> did not turn into a record,
/// and why: a crosslink, a modification that does not resolve to an mzLib modification, a substitution, or
/// a peptide reference the document does not resolve.
/// </summary>
public sealed record MzIdentMLSkippedMatch(string SpectrumIdentificationItemId, string SpectrumId, string Reason);
