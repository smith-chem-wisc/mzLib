using System.Collections.Generic;
using System.Linq;
using Omics;
using ISpectralMatch = Omics.ISpectralMatch; // disambiguate from Readers.ISpectralMatch

namespace Test.Quantification.TestHelpers;

/// <summary>
/// Reads a .psmtsv file and converts records to ISpectralMatch objects
/// suitable for the Quantification pipeline.
/// </summary>
public static class PsmTsvQuantAdapter
{
    /// <summary>
    /// Reads a psmtsv file and returns ISpectralMatch objects with QuantValues populated.
    /// Filters by q-value cutoff and optionally excludes decoys.
    /// </summary>
    /// <param name="psmtsvFilePath">Path to the AllPSMs.psmtsv file</param>
    /// <param name="qValueCutoff">Maximum q-value for filtering (default 0.01)</param>
    /// <param name="includeDecoys">Whether to include decoy matches (default false)</param>
    /// <returns>List of ISpectralMatch objects with QuantValues set</returns>
    public static List<ISpectralMatch> LoadSpectralMatches(
        string psmtsvFilePath,
        double qValueCutoff = 0.01,
        bool includeDecoys = false)
    {
        var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _);
        var results = new List<ISpectralMatch>();

        foreach (var record in psms)
        {
            if (record.QValue > qValueCutoff)
                continue;

            if (!includeDecoys && !record.DecoyContamTarget.Contains('T'))
                continue;

            var match = new BaseSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq,
                identifiedBioPolymers: null)
            {
                QuantValues = record.QuantValues
            };

            results.Add(match);
        }

        return results;
    }

    /// <summary>
    /// Extracts unique protein accessions from PSM records in a psmtsv file.
    /// </summary>
    /// <param name="psmtsvFilePath">Path to the AllPSMs.psmtsv file</param>
    /// <returns>HashSet of unique accession strings</returns>
    public static HashSet<string> GetUniqueAccessions(string psmtsvFilePath)
    {
        var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _);
        return new HashSet<string>(
            psms.Where(p => p.Accession != null)
                .Select(p => p.Accession));
    }
}
