using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Omics.SpectralMatch;
using ISpectralMatch = Omics.SpectralMatch.ISpectralMatch; // disambiguate from Readers.ISpectralMatch

namespace Development.QuantificationDevelopment.TestHelpers;

/// <summary>
/// Reads a .psmtsv file and converts records to ISpectralMatch objects
/// suitable for the Quantification pipeline.
/// </summary>
public static class PsmTsvQuantAdapter
{
    /// <summary>
    /// Reads a psmtsv file and returns ISpectralMatch objects with Intensities populated.
    /// Filters by q-value cutoff and optionally excludes decoys.
    /// </summary>
    /// <param name="psmtsvFilePath">Path to the AllPSMs.psmtsv file</param>
    /// <param name="qValueCutoff">Maximum q-value for filtering (default 0.01)</param>
    /// <param name="includeDecoys">Whether to include decoy matches (default false)</param>
    /// <returns>List of ISpectralMatch objects with Intensities set</returns>
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

            var match = new MockSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq,
                identifiedBioPolymers: null)
            {
                Intensities = record.Intensities
            };

            results.Add(match);
        }

        return results;
    }

    /// <summary>
    /// Reads a psmtsv file and builds all objects needed to run the quantification pipeline:
    /// spectral matches (with identified biopolymers set), peptides, and protein groups.
    /// Only target PSMs passing the q-value cutoff with non-zero Intensities are included.
    /// Multi-protein accessions (pipe-delimited) use the first non-decoy accession.
    /// </summary>
    public static (
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups)
    BuildQuantificationInputs(string psmtsvFilePath, double qValueCutoff = 0.01)
    {
        var records = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
            .Where(r => r.QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && r.Intensities != null
                        && r.Intensities.Any(v => v > 0)
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

        var digestionParams = new DigestionParams(maxMissedCleavages: 100, minPeptideLength: 3);
        var emptyMods = new List<Modification>();

        var baseSeqToPeptide = new Dictionary<string, IBioPolymerWithSetMods>();
        var baseSeqToAccession = new Dictionary<string, string>();

        foreach (var record in records)
        {
            if (baseSeqToPeptide.ContainsKey(record.BaseSeq))
                continue;

            string accession = GetFirstNonDecoyAccession(record.Accession);
            baseSeqToAccession[record.BaseSeq] = accession;

            var protein = new Protein(record.BaseSeq, accession);
            var pep = protein.Digest(digestionParams, emptyMods, emptyMods)
                             .FirstOrDefault(p => p.BaseSequence == record.BaseSeq);

            if (pep != null)
                baseSeqToPeptide[record.BaseSeq] = pep;
        }

        var accessionToPeptides = new Dictionary<string, HashSet<IBioPolymerWithSetMods>>();
        var accessionToProtein  = new Dictionary<string, Protein>();

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
            string accession = kvp.Key;
            var protein = accessionToProtein[accession];
            var pepsForProtein = kvp.Value;

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein },
                pepsForProtein,
                pepsForProtein);

            proteinGroups.Add(group);
        }

        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in records)
        {
            if (!baseSeqToPeptide.TryGetValue(record.BaseSeq, out var peptide))
                continue;

            var match = new MockSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq)
            {
                Intensities = record.Intensities
            };
            match.AddIdentifiedBioPolymer(peptide);

            spectralMatches.Add(match);
        }

        var allPeptides = baseSeqToPeptide.Values.ToList();
        return (spectralMatches, allPeptides, proteinGroups);
    }

    /// <summary>
    /// Reads a search results directory containing AllPSMs.psmtsv, AllPeptides.psmtsv, and
    /// AllQuantifiedProteinGroups.tsv. Each file is filtered independently.
    /// </summary>
    public static (
        List<ISpectralMatch> spectralMatches,
        List<IBioPolymerWithSetMods> peptides,
        List<IBioPolymerGroup> proteinGroups)
    BuildQuantificationInputsFromDirectory(string resultsDirectory, double qValueCutoff = 0.01)
    {
        string psmtsvFilePath = Path.Combine(resultsDirectory, "AllPSMs.psmtsv");
        string peptidesFilePath = Path.Combine(resultsDirectory, "AllPeptides.psmtsv");
        string proteinFilePath = Path.Combine(resultsDirectory, "AllQuantifiedProteinGroups.tsv");

        var passingProteinAccessions = ReadPassingProteinAccessions(proteinFilePath, qValueCutoff);

        var peptideRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(peptidesFilePath, out _)
            .Where(r => !double.IsNaN(r.PEP_QValue)
                        && r.PEP_QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession))
            .ToList();

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

        var psmRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
            .Where(r => !double.IsNaN(r.PEP_QValue)
                        && r.PEP_QValue <= qValueCutoff
                        && r.DecoyContamTarget.Contains('T')
                        && r.Intensities != null
                        && r.Intensities.Any(v => v > 0)
                        && !string.IsNullOrEmpty(r.BaseSeq)
                        && !string.IsNullOrEmpty(r.Accession)
                        && passingPeptideBaseSeqs.Contains(r.BaseSeq))
            .ToList();

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

        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in psmRecords)
        {
            if (!baseSeqToPeptide.TryGetValue(record.BaseSeq, out var peptide))
                continue;

            var match = new MockSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq)
            {
                Intensities = record.Intensities
            };
            match.AddIdentifiedBioPolymer(peptide);
            spectralMatches.Add(match);
        }

        var allPeptides = baseSeqToPeptide.Values.ToList();
        return (spectralMatches, allPeptides, proteinGroups);
    }

    /// <summary>
    /// Combines quantification inputs from multiple search result directories into a single
    /// unified set of spectral matches, peptides, and protein groups.
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

        var baseSeqToPeptide = new Dictionary<string, IBioPolymerWithSetMods>();
        var baseSeqToAccession = new Dictionary<string, string>();

        var allPassingProteinAccessions = new HashSet<string>();
        var allPeptideRecords = new List<Readers.SpectrumMatchFromTsv>();
        var allPsmRecords = new List<Readers.SpectrumMatchFromTsv>();

        foreach (string resultsDirectory in resultsDirectories)
        {
            string psmtsvFilePath = Path.Combine(resultsDirectory, "AllPSMs.psmtsv");
            string peptidesFilePath = Path.Combine(resultsDirectory, "AllPeptides.psmtsv");
            string proteinFilePath = Path.Combine(resultsDirectory, "AllQuantifiedProteinGroups.tsv");

            var passingAccessions = ReadPassingProteinAccessions(proteinFilePath, qValueCutoff);
            allPassingProteinAccessions.UnionWith(passingAccessions);

            var peptideRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(peptidesFilePath, out _)
                .Where(r => !double.IsNaN(r.PEP_QValue)
                            && r.PEP_QValue <= qValueCutoff
                            && r.DecoyContamTarget.Contains('T')
                            && !string.IsNullOrEmpty(r.BaseSeq)
                            && !string.IsNullOrEmpty(r.Accession))
                .ToList();
            allPeptideRecords.AddRange(peptideRecords);

            var psmRecords = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _)
                .Where(r => !double.IsNaN(r.PEP_QValue)
                            && r.PEP_QValue <= qValueCutoff
                            && r.DecoyContamTarget.Contains('T')
                            && r.Intensities != null
                            && r.Intensities.Any(v => v > 0)
                            && !string.IsNullOrEmpty(r.BaseSeq)
                            && !string.IsNullOrEmpty(r.Accession))
                .ToList();
            allPsmRecords.AddRange(psmRecords);
        }

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

        var filteredPsmRecords = allPsmRecords
            .Where(r => passingPeptideBaseSeqs.Contains(r.BaseSeq))
            .ToList();

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

        var spectralMatches = new List<ISpectralMatch>();
        foreach (var record in filteredPsmRecords)
        {
            if (!baseSeqToPeptide.TryGetValue(record.BaseSeq, out var peptide))
                continue;

            var match = new MockSpectralMatch(
                fullFilePath: record.FileNameWithoutExtension,
                oneBasedScanNumber: record.Ms2ScanNumber,
                score: record.Score,
                fullSequence: record.FullSequence,
                baseSequence: record.BaseSeq)
            {
                Intensities = record.Intensities
            };
            match.AddIdentifiedBioPolymer(peptide);
            spectralMatches.Add(match);
        }

        var allPeptides = baseSeqToPeptide.Values.ToList();
        return (spectralMatches, allPeptides, proteinGroups);
    }

    public static HashSet<string> GetUniqueAccessions(string psmtsvFilePath)
    {
        var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(psmtsvFilePath, out _);
        return new HashSet<string>(
            psms.Where(p => p.Accession != null)
                .Select(p => p.Accession));
    }

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

            if (dctIdx >= 0 && dctIdx < fields.Length && !fields[dctIdx].Contains('T'))
                continue;

            if (!double.TryParse(fields[qValueIdx], System.Globalization.NumberStyles.Any,
                    System.Globalization.CultureInfo.InvariantCulture, out double proteinQValue))
                continue;

            if (proteinQValue > qValueCutoff)
                continue;

            string rawAccession = fields[accessionIdx];
            string accession = GetFirstNonDecoyAccession(rawAccession);
            passingAccessions.Add(accession);
        }

        return passingAccessions;
    }

    private static string GetFirstNonDecoyAccession(string rawAccession)
    {
        var parts = rawAccession.Split('|');
        return parts.FirstOrDefault(a => !a.StartsWith("DECOY_")) ?? parts[0];
    }
}
