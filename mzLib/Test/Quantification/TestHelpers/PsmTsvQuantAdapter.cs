using System.Collections.Generic;
using System.IO;
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

        // Maps BaseSeq -> (Protein, PeptideWithSetModifications).
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

        // Build accession -> peptides and accession -> Protein maps
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
    /// Reads a search results directory containing AllPSMs.psmtsv, AllPeptides.psmtsv, and
    /// AllQuantifiedProteinGroups.tsv. Each file is filtered independently:
    /// PSMs and peptides are filtered by PEP_QValue, protein groups by Protein QValue.
    /// Only PSMs whose accessions appear in passing protein groups are included.
    /// Only peptides (by base sequence) that appear in passing protein groups are included.
    /// </summary>
    /// <param name="resultsDirectory">Path to the Task1-SearchTask directory</param>
    /// <param name="qValueCutoff">Maximum q-value for filtering (default 0.01)</param>
    /// <returns>
    /// Tuple of:
    /// - spectralMatches: BaseSpectralMatch list with QuantValues and identified biopolymers set
    /// - peptides: all unique IBioPolymerWithSetMods from passing peptide records
    /// - proteinGroups: one IBioPolymerGroup per passing protein group
    /// </returns>
    public static (
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups)
    BuildQuantificationInputsFromDirectory(string resultsDirectory, double qValueCutoff = 0.01)
    {
        string psmtsvFilePath = Path.Combine(resultsDirectory, "AllPSMs.psmtsv");
        string peptidesFilePath = Path.Combine(resultsDirectory, "AllPeptides.psmtsv");
        string proteinFilePath = Path.Combine(resultsDirectory, "AllQuantifiedProteinGroups.tsv");

        // --- 1. Read protein groups and filter by Protein QValue ---
        var passingProteinAccessions = ReadPassingProteinAccessions(proteinFilePath, qValueCutoff);

        // --- 2. Read peptides (same format as PSMs), filter by PEP_QValue ---
        // The peptides file contains the best PSM per unique peptide sequence.
        // We use it to define the peptide universe, filtered by PEP_QValue.
        var peptideRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(peptidesFilePath, out _)
            .Where(r => !double.IsNaN(r.PEP_QValue)
                        && r.PEP_QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        // Only keep peptides whose accession is in a passing protein group
        var digestionParams = new DigestionParams(maxMissedCleavages: 100, minPeptideLength: 3);
        var emptyMods = new List<Modification>();
        var baseSeqToPeptide = new Dictionary<string, IBioPolymerWithSetMods>();
        var baseSeqToAccession = new Dictionary<string, string>();
        var passingPeptideBaseSeqs = new HashSet<string>();

        foreach (var record in peptideRecords)
        {
            if (baseSeqToPeptide.ContainsKey(record.BaseSeq))
                continue;

            string accession = GetFirstNonDecoyAccession(record.Accession);
            if (!passingProteinAccessions.Contains(accession))
                continue;

            var protein = new Protein(record.BaseSeq, accession);
            var pep = protein.Digest(digestionParams, emptyMods, emptyMods)
                             .FirstOrDefault(p => p.BaseSequence == record.BaseSeq);

            if (pep != null)
            {
                baseSeqToPeptide[record.BaseSeq] = pep;
                baseSeqToAccession[record.BaseSeq] = accession;
                passingPeptideBaseSeqs.Add(record.BaseSeq);
            }
        }

        // --- 3. Read PSMs, filter by PEP_QValue, restrict to passing peptides/proteins ---
        var psmRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
            .Where(r => !double.IsNaN(r.PEP_QValue)
                        && r.PEP_QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && r.QuantValues != null
                        && r.QuantValues.Any(v => v > 0)
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession)
                        && passingPeptideBaseSeqs.Contains(r.BaseSeq))
            .ToList();

        // --- 4. Build accession -> peptides and protein group objects ---
        var accessionToPeptides = new Dictionary<string, HashSet<IBioPolymerWithSetMods>>();
        var accessionToProtein = new Dictionary<string, Protein>();

        foreach (var kvp in baseSeqToPeptide)
        {
            string accession = baseSeqToAccession[kvp.Key];
            if (!accessionToPeptides.ContainsKey(accession))
            {
                accessionToPeptides[accession] = new HashSet<IBioPolymerWithSetMods>();
                accessionToProtein[accession] = (Protein)((IBioPolymerWithSetMods)kvp.Value).Parent;
            }
            accessionToPeptides[accession].Add(kvp.Value);
        }

        var proteinGroups = new List<IBioPolymerGroup>();
        foreach (var kvp in accessionToPeptides)
        {
            var protein = accessionToProtein[kvp.Key];
            var pepsForProtein = kvp.Value;
            proteinGroups.Add(new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                pepsForProtein,
                pepsForProtein));
        }

        // --- 5. Build spectral matches with identified biopolymers ---
        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in psmRecords)
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
    /// Combines quantification inputs from multiple search result directories into a single
    /// unified set of spectral matches, peptides, and protein groups. Shared peptides/proteins
    /// across directories are merged so that a single BioPolymerGroup represents each accession.
    /// </summary>
    public static (
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups)
    BuildQuantificationInputsFromMultipleDirectories(
        IEnumerable<string> resultsDirectories, double qValueCutoff = 0.01)
    {
        var digestionParams = new DigestionParams(maxMissedCleavages: 100, minPeptideLength: 3);
        var emptyMods = new List<Modification>();

        // Merged maps across all directories
        var baseSeqToPeptide = new Dictionary<string, IBioPolymerWithSetMods>();
        var baseSeqToAccession = new Dictionary<string, string>();

        // Collect all passing protein accessions and peptide base sequences from all directories
        var allPassingProteinAccessions = new HashSet<string>();
        var allPeptideRecords = new List<Readers.SpectrumMatchFromTsv>();
        var allPsmRecords = new List<Readers.SpectrumMatchFromTsv>();

        foreach (string resultsDirectory in resultsDirectories)
        {
            string psmtsvFilePath = Path.Combine(resultsDirectory, "AllPSMs.psmtsv");
            string peptidesFilePath = Path.Combine(resultsDirectory, "AllPeptides.psmtsv");
            string proteinFilePath = Path.Combine(resultsDirectory, "AllQuantifiedProteinGroups.tsv");

            // Read passing protein accessions from this directory
            var passingAccessions = ReadPassingProteinAccessions(proteinFilePath, qValueCutoff);
            allPassingProteinAccessions.UnionWith(passingAccessions);

            // Read peptide records
            var peptideRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(peptidesFilePath, out _)
                .Where(r => !double.IsNaN(r.PEP_QValue)
                            && r.PEP_QValue <= qValueCutoff
                            && r.DecoyContamTarget.Contains('T')
                            && !string.IsNullOrEmpty(r.BaseSeq)
                            && !string.IsNullOrEmpty(r.Accession))
                .ToList();
            allPeptideRecords.AddRange(peptideRecords);

            // Read PSM records
            var psmRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
                .Where(r => !double.IsNaN(r.PEP_QValue)
                            && r.PEP_QValue <= qValueCutoff
                            && r.DecoyContamTarget.Contains('T')
                            && r.QuantValues != null
                            && r.QuantValues.Any(v => v > 0)
                            && !string.IsNullOrEmpty(r.BaseSeq)
                            && !string.IsNullOrEmpty(r.Accession))
                .ToList();
            allPsmRecords.AddRange(psmRecords);
        }

        // Build peptide objects from merged peptide records, filtering to passing proteins
        var passingPeptideBaseSeqs = new HashSet<string>();
        foreach (var record in allPeptideRecords)
        {
            if (baseSeqToPeptide.ContainsKey(record.BaseSeq))
                continue;

            string accession = GetFirstNonDecoyAccession(record.Accession);
            if (!allPassingProteinAccessions.Contains(accession))
                continue;

            var protein = new Protein(record.BaseSeq, accession);
            var pep = protein.Digest(digestionParams, emptyMods, emptyMods)
                             .FirstOrDefault(p => p.BaseSequence == record.BaseSeq);

            if (pep != null)
            {
                baseSeqToPeptide[record.BaseSeq] = pep;
                baseSeqToAccession[record.BaseSeq] = accession;
                passingPeptideBaseSeqs.Add(record.BaseSeq);
            }
        }

        // Filter PSMs to passing peptides
        var filteredPsmRecords = allPsmRecords
            .Where(r => passingPeptideBaseSeqs.Contains(r.BaseSeq))
            .ToList();

        // Build protein groups
        var accessionToPeptides = new Dictionary<string, HashSet<IBioPolymerWithSetMods>>();
        var accessionToProtein = new Dictionary<string, Protein>();

        foreach (var kvp in baseSeqToPeptide)
        {
            string accession = baseSeqToAccession[kvp.Key];
            if (!accessionToPeptides.ContainsKey(accession))
            {
                accessionToPeptides[accession] = new HashSet<IBioPolymerWithSetMods>();
                accessionToProtein[accession] = (Protein)((IBioPolymerWithSetMods)kvp.Value).Parent;
            }
            accessionToPeptides[accession].Add(kvp.Value);
        }

        var proteinGroups = new List<IBioPolymerGroup>();
        foreach (var kvp in accessionToPeptides)
        {
            var protein = accessionToProtein[kvp.Key];
            var pepsForProtein = kvp.Value;
            proteinGroups.Add(new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                pepsForProtein,
                pepsForProtein));
        }

        // Build spectral matches
        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in filteredPsmRecords)
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

    /// <summary>
    /// Reads AllQuantifiedProteinGroups.tsv and returns the set of accessions passing the
    /// Protein QValue cutoff. Handles pipe-delimited accessions in the "Protein Accession" column.
    /// </summary>
    private static HashSet<string> ReadPassingProteinAccessions(string proteinGroupsFilePath, double qValueCutoff)
    {
        var passingAccessions = new HashSet<string>();
        using var reader = new StreamReader(proteinGroupsFilePath);

        string headerLine = reader.ReadLine();
        if (headerLine == null) return passingAccessions;

        var headers = headerLine.Split('\t');
        int accessionIdx = System.Array.IndexOf(headers, "Protein Accession");
        int qValueIdx = System.Array.IndexOf(headers, "Protein QValue");
        int dctIdx = System.Array.IndexOf(headers, "Protein Decoy/Contaminant/Target");

        if (accessionIdx < 0 || qValueIdx < 0)
            throw new System.InvalidOperationException(
                $"Protein groups file missing required columns. Found: {headerLine}");

        string line;
        while ((line = reader.ReadLine()) != null)
        {
            var fields = line.Split('\t');
            if (fields.Length <= System.Math.Max(accessionIdx, qValueIdx))
                continue;

            // Filter: target only
            if (dctIdx >= 0 && dctIdx < fields.Length && !fields[dctIdx].Contains('T'))
                continue;

            // Filter: Protein QValue
            if (!double.TryParse(fields[qValueIdx], System.Globalization.NumberStyles.Any,
                    System.Globalization.CultureInfo.InvariantCulture, out double proteinQValue))
                continue;

            if (proteinQValue > qValueCutoff)
                continue;

            // Accession may be pipe-delimited (e.g., "P02768|P02768-2")
            string rawAccession = fields[accessionIdx];
            string accession = GetFirstNonDecoyAccession(rawAccession);
            passingAccessions.Add(accession);
        }

        return passingAccessions;
    }

    private static string GetFirstNonDecoyAccession(string rawAccession)
    {
        // e.g. "O76070|P62937|P62988" -> "O76070"
        // e.g. "DECOY_P16083|P16083" -> "P16083"
        var parts = rawAccession.Split('|');
        return parts.FirstOrDefault(a => !a.StartsWith("DECOY_")) ?? parts[0];
    }
}
