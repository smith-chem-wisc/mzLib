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
/// Development harness for the two MBR error types that the foreign-transfer test in
/// <see cref="TwoProteomeMbrEvaluationTest"/> cannot see, using the in-house E. coli two-proteome
/// data (PXD057758) searched with censoring + an E. coli entrapment DB
/// (folder MM_ConcatenatedHumanDb_Search_CensoredFiles). Both metrics come from ONE MBR run:
///
///  (1) NATIVE-PEAK (mis-)localization — measurable only with the censored files.
///      500 confidently MS2-identified human PSMs were removed from each of the 10 pure-human runs
///      (ground truth: CensoredPsms.psmtsv, with the true retention time). With the MS2 ID gone, MBR
///      must transfer the peptide back. We then ask: did MBR land on the RIGHT peak? For every censored
///      peptide we compare the MBR peak's apex RT to the true RT. |apexRT - trueRT| within tolerance =
///      correct recovery; a large error = a native mis-localization error (the peptide belongs in the
///      run, but MBR quantified the wrong feature). The two-proteome foreign test assumes all native
///      transfers are correct, so it is blind to exactly this.
///
///  (2) PROPAGATED FALSE-ID transfers — measurable with the concatenated entrapment DB.
///      Every human protein in the search DB is paired with an entrapment sequence, so a peptide whose
///      base sequence is in the entrapment list is a false identification by construction. Entrapment
///      peptides are LABELLED "Homo sapiens" in the psmtsv, so they are recognised only by sequence.
///      Entrapment peptides that MBR transfers are propagated false-ID transfers; a better method should
///      not increase them. We report the entrapment fraction at the ID level and after MBR. (E. coli
///      transferred into pure-human runs is reported too, but as a DISTINCT cross-proteome-contamination
///      metric — E. coli is real in the mixed donors, so it is not the entrapment.)
///
/// Run before/after an MBR change: a better method recovers more censored peptides ON the correct peak
/// (higher recovery, lower RT error, lower mis-localization), while not increasing propagated entrapment.
///
/// Local-data, [Explicit], not run in CI. Configure via constants in PipEchoCommon or environment
/// variables (MZLIB_PIPECHO_CENS_PSMTSV, MZLIB_PIPECHO_CENS_SPECTRA, MZLIB_PIPECHO_CENS_TRUTH,
/// MZLIB_PIPECHO_MBR_FDR, MZLIB_PIPECHO_RT_TOL, MZLIB_PIPECHO_MAX_PERGROUP, MZLIB_PIPECHO_OUTDIR).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class CensoredEntrapmentEvaluationTest
{
    private const string MixedMarker = "Ecoli_10to1"; // mixed (human+E.coli) donor runs; others are pure human

    [Test]
    [Explicit("Runs FlashLFQ MBR on the censored+entrapment E. coli set from local disk. Set MZLIB_PIPECHO_CENS_* first.")]
    public void EvaluateNativeRecoveryAndPropagatedTransfers()
    {
        string psmtsv = PipEchoCommon.Env("MZLIB_PIPECHO_CENS_PSMTSV", PipEchoCommon.CensoredEntrapmentPsmtsvDefault);
        string spectraDir = PipEchoCommon.Env("MZLIB_PIPECHO_CENS_SPECTRA", PipEchoCommon.CensoredSpectraDirDefault);
        string truthPath = PipEchoCommon.Env("MZLIB_PIPECHO_CENS_TRUTH", PipEchoCommon.CensoredGroundTruthDefault);
        string entrapmentPath = PipEchoCommon.Env("MZLIB_PIPECHO_ENTRAPMENT", PipEchoCommon.EntrapmentPeptidesDefault);
        double mbrFdr = PipEchoCommon.EnvDouble("MZLIB_PIPECHO_MBR_FDR", 0.01);
        double rtTol = PipEchoCommon.EnvDouble("MZLIB_PIPECHO_RT_TOL", 0.5); // minutes

        Assert.That(File.Exists(psmtsv), Is.True, $"psmtsv not found: {psmtsv}");
        Assert.That(Directory.Exists(spectraDir), Is.True, $"spectra dir not found: {spectraDir}");
        Assert.That(File.Exists(truthPath), Is.True, $"censored ground truth not found: {truthPath}");
        Assert.That(File.Exists(entrapmentPath), Is.True, $"entrapment sequence list not found: {entrapmentPath}");

        var entrapment = PipEchoCommon.LoadEntrapmentSequences(entrapmentPath);
        TestContext.WriteLine($"Loaded {entrapment.Count} entrapment peptide base sequences from {Path.GetFileName(entrapmentPath)}.");

        int maxPerGroup = PipEchoCommon.EnvInt("MZLIB_PIPECHO_MAX_PERGROUP", int.MaxValue);
        var fileInfoByName = BuildFileInfos(spectraDir, maxPerGroup, out var pureFiles, out var mixedFiles);
        var pureSet = new HashSet<SpectraFileInfo>(pureFiles);
        TestContext.WriteLine($"Files: {pureFiles.Count} pure-human, {mixedFiles.Count} mixed (human+E.coli).");
        Assert.That(pureFiles.Count, Is.GreaterThan(0), "No pure-human runs found.");
        Assert.That(mixedFiles.Count, Is.GreaterThan(0), "No mixed donor runs found.");

        var sw = Stopwatch.StartNew();
        var load = PipEchoCommon.StreamIdentifications(psmtsv, fileInfoByName, PipEchoCommon.QValueCutoff);
        sw.Stop();
        int ecoliPeptides = load.PeptideSpecies.Count(kv => kv.Value == Species.Ecoli);
        TestContext.WriteLine($"Loaded {load.Identifications.Count} identifications " +
            $"(kept {load.RowsKept}/{load.RowsRead}; skipped file={load.RowsSkippedNoFile}, q={load.RowsSkippedQ}, " +
            $"ambiguous={load.RowsSkippedAmbiguous}, contaminant={load.RowsSkippedContaminant}) in {sw.Elapsed}. " +
            $"Distinct peptides: {load.PeptideSpecies.Count} ({ecoliPeptides} E. coli).");
        Assert.That(load.Identifications.Count, Is.GreaterThan(0), "No identifications loaded.");

        TestContext.WriteLine($"Running FlashLFQ MBR (fdr={mbrFdr}, threads={PipEchoCommon.MaxThreads})...");
        var engine = new FlashLfqEngine(load.Identifications, matchBetweenRuns: true,
            requireMsmsIdInCondition: false, matchBetweenRunsFdrThreshold: mbrFdr, maxThreads: PipEchoCommon.MaxThreads);
        sw.Restart();
        FlashLfqResults results = engine.Run();
        sw.Stop();
        TestContext.WriteLine($"FlashLFQ MBR finished in {sw.Elapsed}.");

        EvaluateNativeRecovery(results, pureSet, truthPath, mbrFdr, rtTol, entrapment, psmtsv);
        EvaluateEntrapmentPropagation(load, results, pureFiles, entrapment, psmtsv);

        // Sanity: MBR must have recovered at least some censored peptides on-peak.
        // (assertions live inside the helpers via the returned summary)
    }

    // ── (1) Native peak recovery / mis-localization ───────────────────────────────────────────────
    private static void EvaluateNativeRecovery(
        FlashLfqResults results, HashSet<SpectraFileInfo> pureFiles, string truthPath,
        double mbrFdr, double rtTol, HashSet<string> entrapment, string psmtsvPath)
    {
        // For each pure file, index peaks by peptide: whether an MS2 peak exists (censoring incomplete)
        // and the best non-decoy MBR peak (highest MbrScore).
        var byFileSeq = new Dictionary<(string file, string seq), (bool hasMsms, MbrChromatographicPeak? mbr)>();
        foreach (var file in pureFiles)
        {
            if (!results.Peaks.TryGetValue(file, out var peaks)) continue;
            foreach (var peak in peaks)
            {
                if (peak?.Identifications == null || peak.Identifications.Count == 0) continue;
                string seq = peak.Identifications.First().ModifiedSequence;
                var key = (file.FilenameWithoutExtension, seq);
                var cur = byFileSeq.GetValueOrDefault(key);
                if (peak is MbrChromatographicPeak m && !m.RandomRt && !m.DecoyPeptide)
                {
                    if (cur.mbr == null || m.MbrScore > cur.mbr.MbrScore) cur.mbr = m;
                }
                else if (peak.DetectionType == DetectionType.MSMS)
                {
                    cur.hasMsms = true;
                }
                byFileSeq[key] = cur;
            }
        }

        var pureNames = new HashSet<string>(pureFiles.Select(f => f.FilenameWithoutExtension));
        var allTruth = PipEchoCommon.LoadCensoredGroundTruth(truthPath)
            .Where(t => pureNames.Contains(t.FileNoExt)) // only files we actually loaded
            .ToList();
        // Exclude entrapment sequences from the ground truth: "recovering" a false peptide on its
        // (false) MS2 retention time is not a meaningful localization test.
        int nEntrapTruth = allTruth.Count(t => entrapment.Contains(t.BaseSequence));
        var truth = allTruth.Where(t => !entrapment.Contains(t.BaseSequence)).ToList();

        int nStillMsms = 0, nNoPeak = 0, nMbrNoApex = 0, nRecovered = 0, nAcceptedFdr = 0;
        var apexErrors = new List<double>();       // signed apexRT - trueRT for recovered peaks
        var predErrors = new List<double>();       // signed predictedRT - trueRT
        int nMislocalized = 0;                     // recovered but |apexErr| > rtTol
        var rows = new List<(string file, string seq, double trueRt, double apexRt, double predRt, double mbrScore, bool acceptedFdr)>();

        foreach (var t in truth)
        {
            var key = (t.FileNoExt, t.FullSequence);
            if (!byFileSeq.TryGetValue(key, out var hit)) { nNoPeak++; continue; }
            if (hit.hasMsms) { nStillMsms++; continue; }        // not a clean MBR test — peptide still MS2-IDed
            if (hit.mbr == null) { nNoPeak++; continue; }
            if (hit.mbr.Apex == null) { nMbrNoApex++; continue; }

            double apexRt = hit.mbr.ApexRetentionTime;
            double predRt = hit.mbr.PredictedRetentionTime;
            double apexErr = apexRt - t.RetentionTime;
            nRecovered++;
            apexErrors.Add(apexErr);
            predErrors.Add(predRt - t.RetentionTime);
            if (Math.Abs(apexErr) > rtTol) nMislocalized++;

            bool acceptedFdr = results.PeptideModifiedSequences.TryGetValue(t.FullSequence, out var pep)
                               && pep.GetDetectionType(FirstFileByName(pureFiles, t.FileNoExt)) == DetectionType.MBR;
            if (acceptedFdr) nAcceptedFdr++;
            rows.Add((t.FileNoExt, t.FullSequence, t.RetentionTime, apexRt, predRt, hit.mbr.MbrScore, acceptedFdr));
        }

        int nEligible = truth.Count - nStillMsms;
        double recoveryRate = nEligible > 0 ? (double)nRecovered / nEligible : double.NaN;
        double mislocRate = nRecovered > 0 ? (double)nMislocalized / nRecovered : double.NaN;
        var absApex = apexErrors.Select(Math.Abs).ToList();
        var absPred = predErrors.Select(Math.Abs).ToList();

        TestContext.WriteLine("\n=== (1) NATIVE peak recovery of censored human peptides in pure-human runs ===");
        TestContext.WriteLine($"  censored ground-truth peptides (in loaded files): {truth.Count} " +
            $"(excluded {nEntrapTruth} entrapment-sequence ground-truth peptides)");
        TestContext.WriteLine($"  excluded (still MS2-identified, not a clean MBR test): {nStillMsms}  -> eligible {nEligible}");
        TestContext.WriteLine($"  recovered on an MBR peak: {nRecovered}  (recovery {recoveryRate:P1})");
        TestContext.WriteLine($"    of which accepted at MBR FDR {mbrFdr:P1}: {nAcceptedFdr}");
        TestContext.WriteLine($"  not recovered (no peak): {nNoPeak}   MBR peak but no apex: {nMbrNoApex}");
        TestContext.WriteLine($"  |apexRT - trueRT| (min): median={PipEchoCommon.Median(absApex):F4}  " +
            $"p90={Percentile(absApex, 0.90):F4}  mean={(absApex.Count > 0 ? absApex.Average() : double.NaN):F4}");
        TestContext.WriteLine($"  |predictedRT - trueRT| (min): median={PipEchoCommon.Median(absPred):F4}  " +
            $"p90={Percentile(absPred, 0.90):F4}");
        TestContext.WriteLine($"  mis-localized (|apexErr| > {rtTol} min): {nMislocalized}  ({mislocRate:P2} of recovered)");
        foreach (double tol in new[] { 0.1, 0.25, 0.5, 1.0 })
        {
            int within = absApex.Count(e => e <= tol);
            TestContext.WriteLine($"    within {tol,4} min: {within,6}  ({(nRecovered > 0 ? (double)within / nRecovered : 0):P1} of recovered)");
        }

        WriteRecoveryCsv(rows, psmtsvPath);

        Assert.That(nRecovered, Is.GreaterThan(0), "MBR recovered no censored peptides on-peak.");
        Assert.That(recoveryRate, Is.GreaterThan(0.1), $"Suspiciously low censored-peptide recovery ({recoveryRate:P1}).");
    }

    // ── (2) Propagated false-ID (entrapment) transfers ────────────────────────────────────────────
    private static void EvaluateEntrapmentPropagation(
        PipEchoCommon.LoadResult load, FlashLfqResults results,
        List<SpectraFileInfo> pureFiles, HashSet<string> entrapment, string psmtsvPath)
    {
        // Entrapment peptides (base sequence in the entrapment list) are false by construction and false
        // in EVERY run. ID-level: entrapment among target PSMs. MBR-level: entrapment among the transfers
        // MBR produced anywhere — these are the PROPAGATED false-ID transfers.
        int idEntrap = 0, idReal = 0;
        foreach (var id in load.Identifications)
        {
            if (id.IsDecoy) continue;
            if (entrapment.Contains(id.BaseSequence)) idEntrap++; else idReal++;
        }
        int idTotal = idEntrap + idReal;
        double idFrac = idTotal > 0 ? (double)idEntrap / idTotal : double.NaN;

        int mbrEntrap = 0, mbrReal = 0;
        var entrapRows = new List<(string file, string seq, double mbrScore, double? pep)>();
        foreach (var (file, peaks) in results.Peaks)
        {
            foreach (var peak in peaks.OfType<MbrChromatographicPeak>().Where(p => !p.RandomRt && !p.DecoyPeptide))
            {
                var id0 = peak.Identifications.First();
                if (entrapment.Contains(id0.BaseSequence))
                {
                    mbrEntrap++;
                    entrapRows.Add((file.FilenameWithoutExtension, id0.ModifiedSequence, peak.MbrScore, peak.MbrPep));
                }
                else mbrReal++;
            }
        }
        int mbrTotal = mbrEntrap + mbrReal;
        double mbrFrac = mbrTotal > 0 ? (double)mbrEntrap / mbrTotal : double.NaN;

        TestContext.WriteLine("\n=== (2) ENTRAPMENT propagation (false peptide SEQUENCES, labelled human) ===");
        TestContext.WriteLine($"  ID level (target PSMs, q<{PipEchoCommon.QValueCutoff}):");
        TestContext.WriteLine($"    real={idReal}  entrapment(false)={idEntrap}  entrapment fraction={idFrac:P3}");
        TestContext.WriteLine($"  MBR level (all transfers produced, pre-FDR) — propagated false-ID:");
        TestContext.WriteLine($"    real={mbrReal}  entrapment(false)={mbrEntrap}  entrapment fraction={mbrFrac:P3}");

        WriteEntrapmentCsv(entrapRows, psmtsvPath);

        // Secondary, DISTINCT from entrapment: E. coli transferred into pure-human runs is cross-proteome
        // contamination (E. coli is real in the mixed donors, wrong in a pure-human acceptor).
        int mbrEcoli = 0, mbrHuman = 0;
        foreach (var pure in pureFiles)
        {
            if (!results.Peaks.TryGetValue(pure, out var peaks)) continue;
            foreach (var peak in peaks.OfType<MbrChromatographicPeak>().Where(p => !p.RandomRt && !p.DecoyPeptide))
            {
                switch (load.PeptideSpecies.GetValueOrDefault(peak.Identifications.First().ModifiedSequence, Species.Other))
                {
                    case Species.Ecoli: mbrEcoli++; break;
                    case Species.Human: mbrHuman++; break;
                }
            }
        }
        int ecoliAssignable = mbrEcoli + mbrHuman;
        double ecoliFdp = ecoliAssignable > 0 ? (double)mbrEcoli / ecoliAssignable : double.NaN;
        TestContext.WriteLine($"  [secondary] E. coli cross-proteome transfers into pure-human runs: " +
            $"human={mbrHuman}  E.coli={mbrEcoli}  foreign FDP={ecoliFdp:P3}");

        Assert.That(mbrTotal, Is.GreaterThan(0), "MBR produced no transfers.");
    }

    private static void WriteEntrapmentCsv(
        List<(string file, string seq, double mbrScore, double? pep)> rows, string psmtsvPath)
    {
        Directory.CreateDirectory(PipEchoCommon.OutputDir);
        string path = Path.Combine(PipEchoCommon.OutputDir,
            $"EntrapmentTransfers_{DateTime.Now:yyyyMMdd_HHmmss}.csv");
        var sb = new StringBuilder();
        sb.AppendLine("File,EntrapmentPeptide,MbrScore,MbrPep");
        foreach (var r in rows.OrderBy(r => r.file).ThenBy(r => r.seq))
            sb.AppendLine(string.Join(",", r.file, r.seq,
                r.mbrScore.ToString("G17", CultureInfo.InvariantCulture),
                r.pep is double p ? p.ToString("G17", CultureInfo.InvariantCulture) : ""));
        File.WriteAllText(path, sb.ToString());
        TestContext.WriteLine($"\nEntrapment (false) MBR transfers written to: {path}");
    }

    private static SpectraFileInfo FirstFileByName(HashSet<SpectraFileInfo> files, string nameNoExt)
        => files.First(f => f.FilenameWithoutExtension == nameNoExt);

    private static double Percentile(List<double> values, double p)
    {
        if (values.Count == 0) return double.NaN;
        var sorted = values.OrderBy(v => v).ToList();
        int idx = (int)Math.Ceiling(p * sorted.Count) - 1;
        return sorted[Math.Clamp(idx, 0, sorted.Count - 1)];
    }

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
            if (targetList.Count >= maxPerGroup) continue;

            var info = new SpectraFileInfo(path, condition: "proteome", biorep: biorep++, techrep: 0, fraction: 0);

            // The censored search renamed the (censored) PURE runs with a "-censored" suffix in its
            // psmtsv "File Name", while the mixed donor runs kept their names. Censoring only removes MS2
            // identifications, so the UNCENSORED mzML here carries the same MS1 signal — we point the
            // censored identifications at it by registering the pure files under their "-censored" name.
            // The SpectraFileInfo itself keeps the uncensored FilenameWithoutExtension, which matches the
            // CensoredPsms ground-truth "File Name" used for recovery lookup.
            string psmtsvKey = mixed ? nameNoExt : nameNoExt + "-censored";
            fileInfoByName[psmtsvKey] = info;
            targetList.Add(info);
        }
        return fileInfoByName;
    }

    private static void WriteRecoveryCsv(
        List<(string file, string seq, double trueRt, double apexRt, double predRt, double mbrScore, bool acceptedFdr)> rows,
        string psmtsvPath)
    {
        Directory.CreateDirectory(PipEchoCommon.OutputDir);
        string path = Path.Combine(PipEchoCommon.OutputDir,
            $"CensoredRecovery_{DateTime.Now:yyyyMMdd_HHmmss}.csv");
        var sb = new StringBuilder();
        sb.AppendLine("File,Peptide,TrueRT,ApexRT,PredictedRT,ApexRtError,PredRtError,MbrScore,AcceptedAtFdr");
        foreach (var r in rows.OrderByDescending(r => Math.Abs(r.apexRt - r.trueRt)))
        {
            sb.AppendLine(string.Join(",",
                r.file, r.seq,
                r.trueRt.ToString("G17", CultureInfo.InvariantCulture),
                r.apexRt.ToString("G17", CultureInfo.InvariantCulture),
                r.predRt.ToString("G17", CultureInfo.InvariantCulture),
                (r.apexRt - r.trueRt).ToString("G17", CultureInfo.InvariantCulture),
                (r.predRt - r.trueRt).ToString("G17", CultureInfo.InvariantCulture),
                r.mbrScore.ToString("G17", CultureInfo.InvariantCulture),
                r.acceptedFdr));
        }
        File.WriteAllText(path, sb.ToString());
        TestContext.WriteLine($"\nPer-peptide censored recovery written to: {path}");
    }
}
