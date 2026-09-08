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

namespace Development.QuantificationDevelopment.PipEcho;

/// <summary>
/// Development harness for the two-proteome MATCH-BETWEEN-RUNS experiment from the PIP-ECHO paper
/// (https://pmc.ncbi.nlm.nih.gov/articles/PMC12488043/), using the in-house E. coli dataset (PXD057758):
/// 10 pure-human runs and 10 human+E.coli (10:1) mixed runs.
///
/// The logic: MBR transfers peptide identities from donor runs (where a peptide was MS2-identified) to
/// acceptor runs. A PURE human sample contains no E. coli, so any E. coli peptide that MBR quantifies in
/// a pure-human run is, by construction, a FALSE transfer that we can directly observe. The observed
/// foreign-transfer FDP = (E. coli transfers into pure runs) / (all species-assignable transfers into
/// pure runs) is a lower bound on the true MBR error rate, and should track the nominal MBR FDR threshold.
///
/// Run this before and after an MBR change: a better method accepts more transfers overall while keeping
/// the observed foreign FDP at or below the nominal threshold.
///
/// Local-data, [Explicit], not run in CI. Configure via constants in PipEchoCommon or environment
/// variables (MZLIB_PIPECHO_ECOLI_PSMTSV, MZLIB_PIPECHO_ECOLI_SPECTRA, MZLIB_PIPECHO_MBR_FDR,
/// MZLIB_PIPECHO_MAX_PERGROUP, MZLIB_PIPECHO_OUTDIR).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TwoProteomeMbrEvaluationTest
{
    private const string MixedMarker = "Human_Ecoli_10to1"; // donor runs containing both proteomes

    [Test]
    [Explicit("Runs FlashLFQ MBR on the 20-run E. coli two-proteome set from local disk. Set MZLIB_PIPECHO_ECOLI_* first.")]
    public void EvaluateForeignTransferRate()
    {
        string psmtsv = PipEchoCommon.Env("MZLIB_PIPECHO_ECOLI_PSMTSV", PipEchoCommon.EcoliPsmtsvDefault);
        string spectraDir = PipEchoCommon.Env("MZLIB_PIPECHO_ECOLI_SPECTRA", PipEchoCommon.EcoliSpectraDirDefault);
        double mbrFdr = PipEchoCommon.EnvDouble("MZLIB_PIPECHO_MBR_FDR", 0.01);

        Assert.That(File.Exists(psmtsv), Is.True, $"psmtsv not found: {psmtsv}");
        Assert.That(Directory.Exists(spectraDir), Is.True, $"spectra dir not found: {spectraDir}");

        int maxPerGroup = PipEchoCommon.EnvInt("MZLIB_PIPECHO_MAX_PERGROUP", int.MaxValue);
        var fileInfoByName = BuildFileInfos(spectraDir, maxPerGroup, out var pureFiles, out var mixedFiles);

        TestContext.WriteLine($"Files: {pureFiles.Count} pure-human, {mixedFiles.Count} mixed (human+E.coli).");
        Assert.That(pureFiles.Count, Is.GreaterThan(0), "No pure-human acceptor files found.");
        Assert.That(mixedFiles.Count, Is.GreaterThan(0), "No mixed donor files found.");

        var sw = Stopwatch.StartNew();
        var load = PipEchoCommon.StreamIdentifications(psmtsv, fileInfoByName, PipEchoCommon.QValueCutoff);
        sw.Stop();
        int ecoliIds = load.PeptideSpecies.Count(kv => kv.Value == Species.Ecoli);
        TestContext.WriteLine($"Loaded {load.Identifications.Count} identifications " +
            $"(kept {load.RowsKept}/{load.RowsRead}; skipped file={load.RowsSkippedNoFile}, q={load.RowsSkippedQ}, " +
            $"ambiguous={load.RowsSkippedAmbiguous}, contaminant={load.RowsSkippedContaminant}) in {sw.Elapsed}. " +
            $"Distinct peptides: {load.PeptideSpecies.Count} ({ecoliIds} E. coli).");
        Assert.That(load.Identifications.Count, Is.GreaterThan(0), "No identifications loaded.");
        Assert.That(ecoliIds, Is.GreaterThan(0), "No E. coli peptides identified — foreign transfers can't be measured.");

        TestContext.WriteLine($"Running FlashLFQ MBR (fdr={mbrFdr}, threads={PipEchoCommon.MaxThreads})...");
        var engine = new FlashLfqEngine(load.Identifications, matchBetweenRuns: true,
            requireMsmsIdInCondition: false, matchBetweenRunsFdrThreshold: mbrFdr, maxThreads: PipEchoCommon.MaxThreads);
        sw.Restart();
        FlashLfqResults results = engine.Run();
        sw.Stop();
        TestContext.WriteLine($"FlashLFQ MBR finished in {sw.Elapsed}.");

        // We evaluate two transfer populations landing in the pure-human acceptor runs:
        //
        //   (1) ALL MBR transfers produced, pre-FDR — read straight off results.Peaks. This is robust
        //       (non-zero whenever MBR ran at all) and shows the raw foreign contamination.
        //   (2) Transfers ACCEPTED at the nominal MBR FDR — via the peptide-level detection type, which
        //       applies the internal MbrQValue < threshold filter. This is the headline error rate, but
        //       on a small file subset it can legitimately be empty if nothing clears the threshold.
        //
        // A foreign (E. coli) transfer into a pure-human run is, by construction, a false transfer.

        // ── Population (1): pre-FDR, from results.Peaks ───────────────────────────────────────────────
        int rawForeign = 0, rawNative = 0, rawOther = 0;
        var perFileRaw = new Dictionary<string, (int foreign, int native)>();
        var foreignRows = new List<(string file, string peptide, double? pep)>();

        foreach (var pure in pureFiles)
        {
            if (!results.Peaks.TryGetValue(pure, out var peaks))
                continue;
            foreach (var peak in peaks.OfType<MbrChromatographicPeak>().Where(p => !p.RandomRt && !p.DecoyPeptide))
            {
                string modSeq = peak.Identifications.First().ModifiedSequence;
                Species sp = load.PeptideSpecies.GetValueOrDefault(modSeq, Species.Other);
                if (sp == Species.Ecoli)
                {
                    rawForeign++;
                    perFileRaw[pure.FilenameWithoutExtension] = Add(perFileRaw.GetValueOrDefault(pure.FilenameWithoutExtension), 1, 0);
                    foreignRows.Add((pure.FilenameWithoutExtension, modSeq, peak.MbrPep));
                }
                else if (sp == Species.Human)
                {
                    rawNative++;
                    perFileRaw[pure.FilenameWithoutExtension] = Add(perFileRaw.GetValueOrDefault(pure.FilenameWithoutExtension), 0, 1);
                }
                else rawOther++;
            }
        }
        int rawAssignable = rawForeign + rawNative;
        double rawFdp = rawAssignable > 0 ? (double)rawForeign / rawAssignable : double.NaN;

        // ── Population (2): accepted at the nominal MBR FDR, via peptide detection type ────────────────
        int accForeign = 0, accNative = 0, accOther = 0;
        foreach (var (modSeq, peptide) in results.PeptideModifiedSequences)
        {
            Species sp = load.PeptideSpecies.GetValueOrDefault(modSeq, Species.Other);
            foreach (var pure in pureFiles)
            {
                if (peptide.GetDetectionType(pure) != DetectionType.MBR) continue;
                if (sp == Species.Ecoli) accForeign++;
                else if (sp == Species.Human) accNative++;
                else accOther++;
            }
        }
        int accAssignable = accForeign + accNative;
        double accFdp = accAssignable > 0 ? (double)accForeign / accAssignable : double.NaN;

        TestContext.WriteLine($"\n=== Foreign-transfer error into PURE-human runs ===");
        TestContext.WriteLine($"  All MBR transfers produced (pre-FDR):");
        TestContext.WriteLine($"    native(human)={rawNative}  foreign(E.coli)={rawForeign}  unassignable={rawOther}  " +
            $"observed foreign FDP={rawFdp:P3}");
        TestContext.WriteLine($"  Transfers accepted at nominal MBR FDR = {mbrFdr:P2}:");
        TestContext.WriteLine($"    native(human)={accNative}  foreign(E.coli)={accForeign}  unassignable={accOther}  " +
            $"observed foreign FDP={accFdp:P3}");

        TestContext.WriteLine($"\n  Per pure-human run (pre-FDR):");
        foreach (var pure in pureFiles)
        {
            var c = perFileRaw.GetValueOrDefault(pure.FilenameWithoutExtension);
            int tot = c.foreign + c.native;
            TestContext.WriteLine($"    {pure.FilenameWithoutExtension,-60} transfers={tot,6}  foreign={c.foreign,4}  fdp={(tot > 0 ? (double)c.foreign / tot : 0):P2}");
        }

        WriteForeignCsv(foreignRows, psmtsv);

        // Robust across subset sizes: MBR must have produced species-assignable transfers into pure runs.
        Assert.That(rawAssignable, Is.GreaterThan(0),
            "MBR produced no species-assignable transfers into pure-human runs.");
        // Loose ceiling only catches a catastrophic FDR-control regression, not normal drift. Prefer the
        // accepted-at-FDR FDP when the subset is large enough for transfers to clear the threshold.
        double fdpForSanity = accAssignable > 0 ? accFdp : rawFdp;
        Assert.That(fdpForSanity, Is.LessThan(0.5),
            $"Observed foreign FDP ({fdpForSanity:P2}) is implausibly high — MBR FDR control looks broken.");

        static (int foreign, int native) Add((int foreign, int native) c, int f, int n) => (c.foreign + f, c.native + n);
    }

    /// <summary>
    /// Enumerates the calibrated spectra files and splits them into pure-human acceptors and mixed
    /// (human+E.coli) donors by file name. All share one condition so MBR is attempted between every
    /// pair; each file gets its own biological replicate.
    /// </summary>
    private static Dictionary<string, SpectraFileInfo> BuildFileInfos(
        string spectraDir, int maxPerGroup,
        out List<SpectraFileInfo> pureFiles, out List<SpectraFileInfo> mixedFiles)
    {
        var fileInfoByName = new Dictionary<string, SpectraFileInfo>();
        pureFiles = new List<SpectraFileInfo>();
        mixedFiles = new List<SpectraFileInfo>();
        int biorep = 0;

        var files = Directory.EnumerateFiles(spectraDir, "*.mzML")
            .Concat(Directory.EnumerateFiles(spectraDir, "*.raw"))
            .Where(p => p.Contains("Human", StringComparison.OrdinalIgnoreCase))
            .OrderBy(p => p)
            .ToList();

        foreach (var path in files)
        {
            string nameNoExt = Path.GetFileNameWithoutExtension(path);
            bool mixed = nameNoExt.Contains(MixedMarker, StringComparison.OrdinalIgnoreCase);
            var targetList = mixed ? mixedFiles : pureFiles;
            if (targetList.Count >= maxPerGroup)
                continue;

            // Both groups share condition "proteome"; each file is its own biorep so no file is treated
            // as another's technical/biological replicate.
            var info = new SpectraFileInfo(path, condition: "proteome", biorep: biorep++, techrep: 0, fraction: 0);
            fileInfoByName[nameNoExt] = info;
            targetList.Add(info);
        }

        return fileInfoByName;
    }

    private static void WriteForeignCsv(List<(string file, string peptide, double? pep)> foreignRows, string psmtsvPath)
    {
        Directory.CreateDirectory(PipEchoCommon.OutputDir);
        string path = Path.Combine(PipEchoCommon.OutputDir,
            $"TwoProteome_ForeignTransfers_{DateTime.Now:yyyyMMdd_HHmmss}.csv");
        var sb = new StringBuilder();
        sb.AppendLine("PureRun,TransferredEcoliPeptide,MbrPep");
        foreach (var (file, peptide, pep) in foreignRows.OrderBy(r => r.file).ThenBy(r => r.peptide))
            sb.AppendLine($"{file},{peptide},{(pep is double p ? p.ToString("G17", CultureInfo.InvariantCulture) : "")}");
        File.WriteAllText(path, sb.ToString());
        TestContext.WriteLine($"\nForeign (E. coli → pure-human) transfers written to: {path}");
    }
}
