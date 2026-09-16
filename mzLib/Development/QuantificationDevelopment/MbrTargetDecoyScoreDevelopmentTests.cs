using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Development.QuantificationDevelopment;

/// <summary>
/// Development harness that runs FlashLFQ match-between-runs on LOCAL spectra files + a local
/// MetaMorpheus psmtsv, then compares the score distributions of target vs. decoy MBR transfers.
///
/// This is the "run it on your own data" companion to Test.FlashLFQ.MbrTargetDecoyScoreDistributionTest,
/// which is pinned to the tiny sliced test data. Point it at a real dataset to see whether the individual
/// MBR score components (RtScore, PpmScore, IntensityScore, ScanCountScore, IsotopicDistributionScore) and
/// the composite MbrScore actually separate targets from decoys, and to dump every scored peak to a CSV
/// for plotting.
///
/// These tests reference files that only exist on the developer's machine, so they are [Explicit] and are
/// NOT run in CI. Configure them either by editing the constants below or by setting environment variables:
///   MZLIB_MBR_PSMTSV      - full path to an AllPSMs.psmtsv (must include decoy PSMs)
///   MZLIB_MBR_SPECTRA_DIR - directory holding the spectra files named in the psmtsv (.raw / .mzML / .d)
///   MZLIB_MBR_OUTPUT_DIR  - optional; where the per-peak CSV is written (default: system temp)
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MbrTargetDecoyScoreDevelopmentTests
{
    // ── Edit these, or set the matching environment variables, to point at your data ──────────────
    private const string DefaultPsmtsvPath = @"";
    private const string DefaultSpectraDirectory = @"";

    // Q-value cutoff applied to every input PSM (targets and decoys alike), mirroring the filtered
    // identification list MetaMorpheus normally hands to FlashLFQ.
    private const double QValueCutoff = 0.01;

    // Spectra file extensions FlashLFQ can read, in the order we try to resolve a psmtsv file name.
    private static readonly string[] SpectraExtensions = { ".raw", ".mzML", ".mzml", ".d" };

    // ──────────────────────────────────────────────────────────────────────────────────────────────

    [Test]
    [Explicit("Runs FlashLFQ MBR on local files. Set MZLIB_MBR_PSMTSV and MZLIB_MBR_SPECTRA_DIR first.")]
    public void CompareTargetDecoyScoreDistributions()
    {
        string psmtsvPath = Environment.GetEnvironmentVariable("MZLIB_MBR_PSMTSV") ?? DefaultPsmtsvPath;
        string spectraDir = Environment.GetEnvironmentVariable("MZLIB_MBR_SPECTRA_DIR") ?? DefaultSpectraDirectory;
        string outputDir = Environment.GetEnvironmentVariable("MZLIB_MBR_OUTPUT_DIR") ?? Path.GetTempPath();

        Assert.That(string.IsNullOrWhiteSpace(psmtsvPath), Is.False,
            "No psmtsv configured. Set MZLIB_MBR_PSMTSV or DefaultPsmtsvPath.");
        Assert.That(File.Exists(psmtsvPath), Is.True, $"psmtsv not found: {psmtsvPath}");
        Assert.That(Directory.Exists(spectraDir), Is.True, $"Spectra directory not found: {spectraDir}");

        var sw = Stopwatch.StartNew();
        List<Identification> ids = BuildFlashLfqIdentifications(psmtsvPath, spectraDir, QValueCutoff);
        sw.Stop();

        int decoyIdCount = ids.Count(id => id.IsDecoy);
        TestContext.WriteLine($"Loaded {ids.Count} identifications ({ids.Count - decoyIdCount} target, {decoyIdCount} decoy) " +
                              $"across {ids.Select(i => i.FileInfo).Distinct().Count()} spectra files in {sw.Elapsed}.");
        Assert.That(ids.Count, Is.GreaterThan(0), "No identifications built. Check the psmtsv and q-value cutoff.");
        Assert.That(ids.Select(i => i.FileInfo).Distinct().Count(), Is.GreaterThan(1),
            "MBR needs at least two spectra files. Only one file's identifications were resolved.");

        // MBR-enabled run. requireMsmsIdInCondition:false so transfers aren't gated by condition layout,
        // matchBetweenRunsFdrThreshold high so we keep (and thus can inspect) low-scoring transfers too.
        var engine = new FlashLfqEngine(
            ids,
            matchBetweenRuns: true,
            requireMsmsIdInCondition: false,
            matchBetweenRunsFdrThreshold: 0.5,
            maxThreads: Math.Max(1, Environment.ProcessorCount - 1));

        sw.Restart();
        FlashLfqResults results = engine.Run();
        sw.Stop();
        TestContext.WriteLine($"FlashLFQ run finished in {sw.Elapsed}.");

        List<MbrChromatographicPeak> mbrPeaks = results.Peaks
            .SelectMany(kvp => kvp.Value)
            .OfType<MbrChromatographicPeak>()
            .ToList();

        Assert.That(mbrPeaks.Count, Is.GreaterThan(0),
            "No MBR transfers were produced. The files may not overlap in peptides, or MBR could not be scored.");

        AnalyzeAndReport(mbrPeaks, psmtsvPath, outputDir);
    }

    #region Analysis

    private static readonly (string Name, Func<MbrChromatographicPeak, double> Selector)[] ScoreComponents =
    {
        ("MbrScore",                  p => p.MbrScore),
        ("IntensityScore",            p => p.IntensityScore),
        ("RtScore",                   p => p.RtScore),
        ("PpmScore",                  p => p.PpmScore),
        ("ScanCountScore",            p => p.ScanCountScore),
        ("IsotopicDistributionScore", p => p.IsotopicDistributionScore),
    };

    /// <summary>
    /// Prints per-component target/decoy separation statistics and writes every scored peak to a CSV.
    /// Asserts only that the composite MbrScore ranks targets above decoys; sub-scores are diagnostic.
    /// </summary>
    private static void AnalyzeAndReport(List<MbrChromatographicPeak> mbrPeaks, string psmtsvPath, string outputDir)
    {
        // Four categories from (DecoyPeptide, RandomRt), matching FlashLfqEngine.CalculateFdrForMbrPeaks.
        var targets = mbrPeaks.Where(p => !p.DecoyPeptide && !p.RandomRt).ToList();
        var randomRtDecoys = mbrPeaks.Where(p => !p.DecoyPeptide && p.RandomRt).ToList();
        var decoyPeptides = mbrPeaks.Where(p => p.DecoyPeptide && !p.RandomRt).ToList();
        var doubleDecoys = mbrPeaks.Where(p => p.DecoyPeptide && p.RandomRt).ToList();
        var allDecoys = mbrPeaks.Where(p => p.DecoyPeptide || p.RandomRt).ToList();

        TestContext.WriteLine($"\nTotal MBR peaks: {mbrPeaks.Count}");
        TestContext.WriteLine($"  Targets:            {targets.Count}");
        TestContext.WriteLine($"  Random-RT decoys:   {randomRtDecoys.Count}");
        TestContext.WriteLine($"  Decoy peptides:     {decoyPeptides.Count}");
        TestContext.WriteLine($"  Double decoys:      {doubleDecoys.Count}");
        TestContext.WriteLine($"  All decoys:         {allDecoys.Count}");

        Assert.That(targets.Count, Is.GreaterThan(0), "No target MBR transfers were produced.");
        Assert.That(allDecoys.Count, Is.GreaterThan(0),
            "No decoy MBR transfers were produced, so target/decoy distributions cannot be compared. " +
            "Ensure the psmtsv includes decoy PSMs.");

        TestContext.WriteLine();
        TestContext.WriteLine($"{"Score component",-28}{"targetMed",12}{"decoyMed",12}{"rtDecoyMed",12}{"pepDecoyMed",12}{"AUC(t>d)",10}");

        double overallAuc = double.NaN;
        foreach (var (name, selector) in ScoreComponents)
        {
            List<double> targetScores = targets.Select(selector).ToList();
            List<double> decoyScores = allDecoys.Select(selector).ToList();

            double targetMedian = Median(targetScores);
            double decoyMedian = Median(decoyScores);
            double rtDecoyMedian = randomRtDecoys.Count > 0 ? Median(randomRtDecoys.Select(selector).ToList()) : double.NaN;
            double pepDecoyMedian = decoyPeptides.Count > 0 ? Median(decoyPeptides.Select(selector).ToList()) : double.NaN;
            double auc = RankSumAuc(targetScores, decoyScores);

            TestContext.WriteLine($"{name,-28}{targetMedian,12:F4}{decoyMedian,12:F4}{rtDecoyMedian,12:F4}{pepDecoyMedian,12:F4}{auc,10:F3}");

            if (name == "MbrScore")
                overallAuc = auc;
        }

        string csvPath = WritePeakCsv(mbrPeaks, psmtsvPath, outputDir);
        TestContext.WriteLine($"\nPer-peak scores written to: {csvPath}");

        Assert.That(overallAuc, Is.GreaterThan(0.5),
            "Targets did not rank above decoys on the composite MBR score (AUC <= 0.5).");
    }

    /// <summary>
    /// Writes one row per scored MBR peak with its classification and every score component, so the
    /// distributions can be plotted or analyzed externally.
    /// </summary>
    private static string WritePeakCsv(List<MbrChromatographicPeak> mbrPeaks, string psmtsvPath, string outputDir)
    {
        Directory.CreateDirectory(outputDir);
        string csvPath = Path.Combine(outputDir,
            $"MbrTargetDecoyScores_{Path.GetFileNameWithoutExtension(psmtsvPath)}_{DateTime.Now:yyyyMMdd_HHmmss}.csv");

        var sb = new StringBuilder();
        sb.AppendLine("File,ModifiedSequence,Class,DecoyPeptide,RandomRt,MbrPep,MbrScore," +
                      "IntensityScore,RtScore,PpmScore,ScanCountScore,IsotopicDistributionScore");

        foreach (var p in mbrPeaks.OrderByDescending(p => p.MbrScore))
        {
            string cls = (p.DecoyPeptide, p.RandomRt) switch
            {
                (false, false) => "Target",
                (true, false) => "DecoyPeptide",
                (false, true) => "RandomRtDecoy",
                (true, true) => "DoubleDecoy",
            };

            sb.Append(Csv(p.SpectraFileInfo.FilenameWithoutExtension)).Append(',');
            sb.Append(Csv(p.Identifications.First().ModifiedSequence)).Append(',');
            sb.Append(cls).Append(',');
            sb.Append(p.DecoyPeptide).Append(',');
            sb.Append(p.RandomRt).Append(',');
            sb.Append(p.MbrPep is double pep ? F(pep) : "").Append(',');
            sb.Append(F(p.MbrScore)).Append(',');
            sb.Append(F(p.IntensityScore)).Append(',');
            sb.Append(F(p.RtScore)).Append(',');
            sb.Append(F(p.PpmScore)).Append(',');
            sb.Append(F(p.ScanCountScore)).Append(',');
            sb.Append(F(p.IsotopicDistributionScore)).AppendLine();
        }

        File.WriteAllText(csvPath, sb.ToString());
        return csvPath;

        static string F(double d) => d.ToString("G17", CultureInfo.InvariantCulture);
        static string Csv(string s) => s.Contains(',') || s.Contains('"')
            ? "\"" + s.Replace("\"", "\"\"") + "\""
            : s;
    }

    #endregion

    #region Input building

    /// <summary>
    /// Reads a MetaMorpheus psmtsv and builds FlashLFQ Identifications (targets and decoys), mapping each
    /// PSM's file name to a spectra file in <paramref name="spectraDirectory"/>. PSMs whose spectra file
    /// cannot be found on disk are skipped (with a one-time warning per missing file).
    /// </summary>
    private static List<Identification> BuildFlashLfqIdentifications(string psmtsvPath, string spectraDirectory, double qValueCutoff)
    {
        var records = SpectrumMatchTsvReader.ReadPsmTsv(psmtsvPath, out var warnings);
        if (warnings != null && warnings.Count > 0)
            TestContext.WriteLine($"psmtsv reader reported {warnings.Count} warnings (first: {warnings[0]}).");

        var fileInfoByName = new Dictionary<string, SpectraFileInfo>(StringComparer.OrdinalIgnoreCase);
        var proteinGroupsByName = new Dictionary<string, ProteinGroup>();
        var missingFilesWarned = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

        var ids = new List<Identification>();
        int biologicalReplicate = 0;

        foreach (var record in records)
        {
            if (!double.IsNaN(record.QValue) && record.QValue > qValueCutoff)
                continue;
            if (string.IsNullOrEmpty(record.BaseSeq) || string.IsNullOrEmpty(record.FullSequence))
                continue;
            if (record.FullSequence.Contains('|')) // ambiguous identification; skip
                continue;

            string fileName = record.FileNameWithoutExtension;
            if (!fileInfoByName.TryGetValue(fileName, out var fileInfo))
            {
                string? resolvedPath = ResolveSpectraFile(spectraDirectory, fileName);
                if (resolvedPath == null)
                {
                    if (missingFilesWarned.Add(fileName))
                        TestContext.WriteLine($"  WARNING: no spectra file for '{fileName}' in {spectraDirectory}; its PSMs are skipped.");
                    continue;
                }

                // Give every file its own biological replicate so no file is treated as an MS/MS-required
                // condition partner. Condition is shared so MBR is attempted between all files.
                fileInfo = new SpectraFileInfo(resolvedPath, condition: "condition", biorep: biologicalReplicate++, techrep: 0, fraction: 0);
                fileInfoByName[fileName] = fileInfo;
            }

            var proteinGroups = new List<ProteinGroup>();
            foreach (var accession in (record.Accession ?? "").Split('|', StringSplitOptions.RemoveEmptyEntries))
            {
                if (!proteinGroupsByName.TryGetValue(accession, out var pg))
                {
                    pg = new ProteinGroup(accession, "", "");
                    proteinGroupsByName[accession] = pg;
                }
                proteinGroups.Add(pg);
            }

            ids.Add(new Identification(
                fileInfo,
                record.BaseSeq,
                record.FullSequence,
                record.MonoisotopicMass,
                record.RetentionTime,
                record.PrecursorCharge,
                proteinGroups,
                psmScore: record.Score,
                qValue: double.IsNaN(record.QValue) ? 0 : record.QValue,
                decoy: record.IsDecoy));
        }

        return ids;
    }

    /// <summary>
    /// Finds a spectra file on disk matching a psmtsv file name, trying each supported extension.
    /// </summary>
    private static string? ResolveSpectraFile(string spectraDirectory, string fileNameWithoutExtension)
    {
        foreach (var ext in SpectraExtensions)
        {
            string candidate = Path.Combine(spectraDirectory, fileNameWithoutExtension + ext);
            if (File.Exists(candidate) || Directory.Exists(candidate)) // .d is a directory
                return candidate;
        }
        return null;
    }

    #endregion

    #region Statistics

    private static double Median(List<double> values)
    {
        if (values.Count == 0) return double.NaN;
        var sorted = values.Where(v => !double.IsNaN(v)).OrderBy(v => v).ToList();
        if (sorted.Count == 0) return double.NaN;
        int mid = sorted.Count / 2;
        return sorted.Count % 2 == 1 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2.0;
    }

    /// <summary>
    /// Rank-based AUC: the probability that a randomly chosen target score exceeds a randomly chosen
    /// decoy score (ties count as half). 0.5 means indistinguishable, 1.0 means perfect separation with
    /// targets on top. Equivalent to the normalized Mann-Whitney U statistic.
    /// </summary>
    private static double RankSumAuc(IReadOnlyList<double> targetScores, IReadOnlyList<double> decoyScores)
    {
        if (targetScores.Count == 0 || decoyScores.Count == 0)
            return double.NaN;

        var pooled = targetScores.Select(s => (score: s, isTarget: true))
            .Concat(decoyScores.Select(s => (score: s, isTarget: false)))
            .OrderBy(x => x.score)
            .ToList();

        double[] ranks = new double[pooled.Count];
        int i = 0;
        while (i < pooled.Count)
        {
            int j = i;
            while (j + 1 < pooled.Count && pooled[j + 1].score == pooled[i].score)
                j++;
            double averageRank = ((i + 1) + (j + 1)) / 2.0; // ranks are 1-based
            for (int k = i; k <= j; k++)
                ranks[k] = averageRank;
            i = j + 1;
        }

        double targetRankSum = 0;
        for (int k = 0; k < pooled.Count; k++)
            if (pooled[k].isTarget)
                targetRankSum += ranks[k];

        int nTarget = targetScores.Count;
        int nDecoy = decoyScores.Count;
        double u = targetRankSum - (double)nTarget * (nTarget + 1) / 2.0;
        return u / ((double)nTarget * nDecoy);
    }

    #endregion
}
