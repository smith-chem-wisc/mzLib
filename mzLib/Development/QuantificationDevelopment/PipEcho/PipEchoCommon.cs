using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;

namespace Development.QuantificationDevelopment.PipEcho;

/// <summary>
/// Which proteome a peptide/protein belongs to in the two-proteome and spike-in datasets from the
/// PIP-ECHO paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC12488043/).
/// </summary>
public enum Species { Human, Ecoli, Other }

/// <summary>
/// Shared infrastructure for the PIP-ECHO development evaluation tests. These tests reference large
/// local data (multi-GB psmtsv + ~1 GB calibrated mzML files each) that only exists on the developer's
/// machine, so the tests that use this are [Explicit] and never run in CI.
///
/// The psmtsv files are 1.6-2.4 GB, so identifications are STREAMED and filtered line-by-line rather
/// than loaded with SpectrumMatchTsvReader (which would materialize every row in memory).
/// </summary>
public static class PipEchoCommon
{
    // ── Default data locations (override any with the matching environment variable) ──────────────
    public const string SpikeInPsmtsvDefault =
        @"D:\PIP_ECHO_PRIDE\ShenIonstarDataset\IonStar_SearchResults-MetaMorpheus\Task2-SearchTask\AllPSMs.psmtsv";
    public const string SpikeInSpectraDirDefault =
        @"D:\PIP_ECHO_PRIDE\ShenIonstarDataset\CalibratedFiles-MetaMorpheus";
    public const string SpikeInExperimentalDesignDefault =
        @"D:\PIP_ECHO_PRIDE\ShenIonstarDataset\CalibratedFiles-MetaMorpheus\ExperimentalDesign.tsv";

    public const string EcoliPsmtsvDefault =
        @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_SearchResults-MetaMorpheus\MM_ConcatenatedHumanDb_Search\Task1-SearchTask\AllPSMs.psmtsv";
    public const string EcoliSpectraDirDefault =
        @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_CalibratedFiles-MetaMorpheus";

    public static string Env(string key, string fallback)
    {
        string? v = Environment.GetEnvironmentVariable(key);
        return string.IsNullOrWhiteSpace(v) ? fallback : v;
    }

    public static int EnvInt(string key, int fallback)
    {
        string? v = Environment.GetEnvironmentVariable(key);
        return int.TryParse(v, out int i) ? i : fallback;
    }

    public static double EnvDouble(string key, double fallback)
    {
        string? v = Environment.GetEnvironmentVariable(key);
        return double.TryParse(v, NumberStyles.Any, CultureInfo.InvariantCulture, out double d) ? d : fallback;
    }

    public static string OutputDir => Env("MZLIB_PIPECHO_OUTDIR", Path.GetTempPath());
    public static double QValueCutoff => EnvDouble("MZLIB_PIPECHO_QCUTOFF", 0.01);
    public static int MaxThreads => EnvInt("MZLIB_PIPECHO_MAXTHREADS", Math.Max(1, Environment.ProcessorCount - 1));

    private static readonly string[] SpectraExtensions = { ".mzML", ".mzml", ".raw", ".d" };

    /// <summary>
    /// Finds the spectra file on disk whose name (sans extension) matches a psmtsv "File Name".
    /// </summary>
    public static string? ResolveSpectraFile(string spectraDir, string fileNameNoExt)
    {
        foreach (var ext in SpectraExtensions)
        {
            string candidate = Path.Combine(spectraDir, fileNameNoExt + ext);
            if (File.Exists(candidate) || Directory.Exists(candidate)) // .d is a directory
                return candidate;
        }
        return null;
    }

    /// <summary>
    /// Classifies a MetaMorpheus "Organism Name" field (which is pipe-delimited for shared peptides).
    /// A peptide is E. coli or human only if EVERY organism token agrees; mixed or contaminant
    /// (e.g. Bos taurus) rows are Species.Other and excluded from the species-specific analyses.
    /// </summary>
    public static Species ClassifyOrganism(string? organismName)
    {
        if (string.IsNullOrWhiteSpace(organismName))
            return Species.Other;

        var tokens = organismName
            .Split('|', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Distinct()
            .ToList();
        if (tokens.Count == 0)
            return Species.Other;

        if (tokens.All(t => t.Contains("Escherichia coli", StringComparison.OrdinalIgnoreCase)))
            return Species.Ecoli;
        if (tokens.All(t => t.Equals("Homo sapiens", StringComparison.OrdinalIgnoreCase)))
            return Species.Human;
        return Species.Other;
    }

    /// <summary>
    /// The identifications and species lookups produced by streaming a psmtsv.
    /// </summary>
    public sealed class LoadResult
    {
        public List<Identification> Identifications { get; } = new();
        /// <summary>Full (modified) sequence → species. Conflicting rows collapse to Other.</summary>
        public Dictionary<string, Species> PeptideSpecies { get; } = new();
        /// <summary>Protein accession → species. Conflicting rows collapse to Other.</summary>
        public Dictionary<string, Species> AccessionSpecies { get; } = new();
        public int RowsRead, RowsKept, RowsSkippedNoFile, RowsSkippedQ, RowsSkippedAmbiguous, RowsSkippedContaminant;
    }

    /// <summary>
    /// Streams a MetaMorpheus AllPSMs.psmtsv and builds FlashLFQ Identifications for every PSM that
    /// (a) belongs to one of the requested spectra files, (b) passes the q-value cutoff, and (c) is an
    /// unambiguous, non-contaminant sequence. Targets and decoys are both kept (decoys flagged), which
    /// mirrors what MetaMorpheus hands FlashLFQ. Species lookups are populated from the Organism column.
    /// </summary>
    public static LoadResult StreamIdentifications(
        string psmtsvPath,
        IReadOnlyDictionary<string, SpectraFileInfo> fileInfoByName,
        double qValueCutoff)
    {
        var result = new LoadResult();
        var proteinGroupCache = new Dictionary<string, ProteinGroup>();

        using var reader = new StreamReader(psmtsvPath);
        string? header = reader.ReadLine();
        if (header == null)
            return result;

        var columns = header.Split('\t');
        int Idx(string name)
        {
            int i = Array.IndexOf(columns, name);
            if (i < 0) throw new InvalidOperationException($"Column '{name}' not found in {psmtsvPath}");
            return i;
        }

        int cFile = Idx("File Name");
        int cBase = Idx("Base Sequence");
        int cFull = Idx("Full Sequence");
        int cMass = Idx("Peptide Monoisotopic Mass");
        int cCharge = Idx("Precursor Charge");
        int cRt = Idx("Scan Retention Time");
        int cQ = Idx("QValue");
        int cAcc = Idx("Protein Accession");
        int cOrg = Idx("Organism Name");
        int cDct = Idx("Decoy/Contaminant/Target");
        int cScore = Idx("Score");
        int maxIdx = new[] { cFile, cBase, cFull, cMass, cCharge, cRt, cQ, cAcc, cOrg, cDct, cScore }.Max();

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            result.RowsRead++;
            var f = line.Split('\t');
            if (f.Length <= maxIdx)
                continue;

            if (!fileInfoByName.TryGetValue(f[cFile], out var fileInfo))
            {
                result.RowsSkippedNoFile++;
                continue;
            }

            string dct = f[cDct];
            bool isDecoy = dct.Contains('D');
            if (dct.Contains('C')) { result.RowsSkippedContaminant++; continue; } // drop contaminants

            if (!double.TryParse(f[cQ], NumberStyles.Any, CultureInfo.InvariantCulture, out double q))
                q = 1.0;
            if (q > qValueCutoff) { result.RowsSkippedQ++; continue; }

            string baseSeq = f[cBase];
            string fullSeq = f[cFull];
            if (string.IsNullOrEmpty(baseSeq) || string.IsNullOrEmpty(fullSeq))
                continue;
            if (fullSeq.Contains('|')) { result.RowsSkippedAmbiguous++; continue; } // ambiguous PSM

            if (!double.TryParse(f[cMass], NumberStyles.Any, CultureInfo.InvariantCulture, out double mass))
                continue;
            int charge = double.TryParse(f[cCharge], NumberStyles.Any, CultureInfo.InvariantCulture, out double z) ? (int)Math.Round(z) : 0;
            double rt = double.TryParse(f[cRt], NumberStyles.Any, CultureInfo.InvariantCulture, out double rtv) ? rtv : 0;
            double score = double.TryParse(f[cScore], NumberStyles.Any, CultureInfo.InvariantCulture, out double sv) ? sv : 0;

            Species species = ClassifyOrganism(f[cOrg]);
            Record(result.PeptideSpecies, fullSeq, species);

            var accessions = f[cAcc].Split('|', StringSplitOptions.RemoveEmptyEntries);
            var proteinGroups = new List<ProteinGroup>(accessions.Length);
            foreach (var acc in accessions)
            {
                if (!proteinGroupCache.TryGetValue(acc, out var pg))
                {
                    pg = new ProteinGroup(acc, "", "");
                    proteinGroupCache[acc] = pg;
                }
                proteinGroups.Add(pg);
                Record(result.AccessionSpecies, acc, species);
            }

            result.Identifications.Add(new Identification(
                fileInfo, baseSeq, fullSeq, mass, rt, charge, proteinGroups,
                psmScore: score, qValue: q, decoy: isDecoy));
            result.RowsKept++;
        }

        return result;

        static void Record(Dictionary<string, Species> map, string key, Species species)
        {
            if (map.TryGetValue(key, out var prev))
            {
                if (prev != species) map[key] = Species.Other; // conflicting assignment → ambiguous
            }
            else
            {
                map[key] = species;
            }
        }
    }

    public static double Median(IReadOnlyList<double> values)
    {
        var sorted = values.Where(v => !double.IsNaN(v) && !double.IsInfinity(v)).OrderBy(v => v).ToList();
        if (sorted.Count == 0) return double.NaN;
        int mid = sorted.Count / 2;
        return sorted.Count % 2 == 1 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2.0;
    }

    public static double Log2(double x) => Math.Log(x, 2);
}
