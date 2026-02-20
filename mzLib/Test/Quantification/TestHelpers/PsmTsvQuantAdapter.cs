using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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
    /// Reads a psmtsv file and builds all objects needed to run the quantification pipeline:
    /// spectral matches (with identified biopolymers set), peptides, and protein groups.
    /// Only target PSMs passing the q-value cutoff with non-zero QuantValues are included.
    /// Multi-protein accessions (pipe-delimited) use the first non-decoy accession.
    /// </summary>
    /// <param name="psmtsvFilePath">Path to the AllPSMs.psmtsv file</param>
    /// <param name="qValueCutoff">Maximum q-value for filtering (default 0.01)</param>
    /// <returns>
    /// Tuple of:
    /// - spectralMatches: BaseSpectralMatch list with QuantValues and identified biopolymers set
    /// - peptides: all unique IBioPolymerWithSetMods created from the PSM records
    /// - proteinGroups: one IBioPolymerGroup per unique accession
    /// </returns>
    public static (
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups)
    BuildQuantificationInputs(string psmtsvFilePath, double qValueCutoff = 0.01)
    {
        var records = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
            .Where(r => r.QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && r.QuantValues != null
                        && r.QuantValues.Any(v => v > 0)
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        // Digestion parameters that allow recovering the full peptide even if it contains
        // internal K/R (e.g., missed cleavages present in original search).
        var digestionParams = new DigestionParams(maxMissedCleavages: 100, minPeptideLength: 3);
        var emptyMods = new List<Modification>();

        // Maps BaseSeq → (Protein, PeptideWithSetModifications).
        // We key by BaseSeq because all PSMs with the same base sequence should produce
        // the same peptide object regardless of TMT modification notation in FullSequence.
        var baseSeqToPeptide = new Dictionary<string, IBioPolymerWithSetMods>();
        var baseSeqToAccession = new Dictionary<string, string>();

        foreach (var record in records)
        {
            if (baseSeqToPeptide.ContainsKey(record.BaseSeq))
                continue;

            string accession = GetFirstNonDecoyAccession(record.Accession);
            baseSeqToAccession[record.BaseSeq] = accession;

            // Create a protein whose entire sequence IS the peptide's base sequence.
            // Digestion with high missed cleavages always yields the full sequence as one peptide.
            var protein = new Protein(record.BaseSeq, accession);
            var pep = protein.Digest(digestionParams, emptyMods, emptyMods)
                             .FirstOrDefault(p => p.BaseSequence == record.BaseSeq);

            if (pep != null)
                baseSeqToPeptide[record.BaseSeq] = pep;
        }

        // Build accession → peptides and accession → Protein maps
        var accessionToPeptides = new Dictionary<string, HashSet<IBioPolymerWithSetMods>>();
        var accessionToProtein  = new Dictionary<string, Protein>();

        foreach (var kvp in baseSeqToPeptide)
        {
            string accession = baseSeqToAccession[kvp.Key];
            if (!accessionToPeptides.ContainsKey(accession))
            {
                accessionToPeptides[accession] = new HashSet<IBioPolymerWithSetMods>();
                // Use the parent protein from the first peptide for this accession
                accessionToProtein[accession] = (Protein)((IBioPolymerWithSetMods)kvp.Value).Parent;
            }
            accessionToPeptides[accession].Add(kvp.Value);
        }

        // Create one BioPolymerGroup per unique accession
        var proteinGroups = new List<IBioPolymerGroup>();
        foreach (var kvp in accessionToPeptides)
        {
            string accession = kvp.Key;
            var protein = accessionToProtein[accession];
            var pepsForProtein = kvp.Value;

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                pepsForProtein,
                pepsForProtein);

            proteinGroups.Add(group);
        }

        // Build spectral matches with identified biopolymers set
        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in records)
        {
            if (!baseSeqToPeptide.TryGetValue(record.BaseSeq, out var peptide))
                continue;

            var match = new BaseSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq)
            {
                QuantValues = record.QuantValues
            };
            match.AddIdentifiedBioPolymer(peptide);

            spectralMatches.Add(match);
        }

        var allPeptides = baseSeqToPeptide.Values.ToList();
        return (spectralMatches, allPeptides, proteinGroups);
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

    private static string GetFirstNonDecoyAccession(string rawAccession)
    {
        // e.g. "O76070|P62937|P62988" → "O76070"
        // e.g. "DECOY_P16083|P16083" → "P16083"
        var parts = rawAccession.Split('|');
        return parts.FirstOrDefault(a => !a.StartsWith("DECOY_")) ?? parts[0];
    }
}
