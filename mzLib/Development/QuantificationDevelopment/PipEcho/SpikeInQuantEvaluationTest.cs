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
/// Development harness for the Shen "IonStar" E. coli-in-human SPIKE-IN experiment (PXD003881), as used
/// in the PIP-ECHO paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC12488043/).
///
/// Human is held constant across all runs; E. coli is spiked at five increasing levels, conditions A–E:
/// A = 1×, B = 1.5×, C = 2×, D = 2.5×, E = 3× (4 biological replicates each, 20 runs total). After
/// quantification, E. coli proteins should recover the known fold-change between conditions
/// (log2(levelHi/levelLo)), while human proteins should stay flat (log2 FC ≈ 0).
///
/// The test reports, per condition comparison: the observed median log2 FC for E. coli vs. human, the
/// error against the expected E. coli FC, and how many E. coli proteins are recovered before the false
/// discovery proportion (human proteins ranked as "changed") exceeds 5%. Run it before and after a
/// quant/MBR change to see whether accuracy improves.
///
/// Local-data, [Explicit], not run in CI. Configure via constants in PipEchoCommon or environment
/// variables (MZLIB_PIPECHO_SPIKEIN_PSMTSV, MZLIB_PIPECHO_SPIKEIN_SPECTRA,
/// MZLIB_PIPECHO_SPIKEIN_DESIGN, MZLIB_PIPECHO_CONDITIONS, MZLIB_PIPECHO_MAX_BIOREPS,
/// MZLIB_PIPECHO_MBR, MZLIB_PIPECHO_NORMALIZE, MZLIB_PIPECHO_OUTDIR).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SpikeInQuantEvaluationTest
{
    // E. coli relative spike level per condition. Human is constant (expected FC = 1, log2 = 0).
    private static readonly Dictionary<string, double> ConditionLevel = new()
    {
        ["A"] = 1.0, ["B"] = 1.5, ["C"] = 2.0, ["D"] = 2.5, ["E"] = 3.0,
    };

    [Test]
    [Explicit("Runs FlashLFQ on the 20-run Shen spike-in from local disk. Set MZLIB_PIPECHO_SPIKEIN_* first.")]
    public void EvaluateSpikeInQuantificationAccuracy()
    {
        string psmtsv = PipEchoCommon.Env("MZLIB_PIPECHO_SPIKEIN_PSMTSV", PipEchoCommon.SpikeInPsmtsvDefault);
        string spectraDir = PipEchoCommon.Env("MZLIB_PIPECHO_SPIKEIN_SPECTRA", PipEchoCommon.SpikeInSpectraDirDefault);
        string designPath = PipEchoCommon.Env("MZLIB_PIPECHO_SPIKEIN_DESIGN", PipEchoCommon.SpikeInExperimentalDesignDefault);

        Assert.That(File.Exists(psmtsv), Is.True, $"psmtsv not found: {psmtsv}");
        Assert.That(Directory.Exists(spectraDir), Is.True, $"spectra dir not found: {spectraDir}");
        Assert.That(File.Exists(designPath), Is.True, $"experimental design not found: {designPath}");

        // Optional subsetting for quick iteration.
        var conditionFilter = PipEchoCommon.Env("MZLIB_PIPECHO_CONDITIONS", "")
            .Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(s => s.ToUpperInvariant()).ToHashSet();
        int maxBioreps = PipEchoCommon.EnvInt("MZLIB_PIPECHO_MAX_BIOREPS", int.MaxValue);

        var fileInfoByName = BuildFileInfos(designPath, spectraDir, conditionFilter, maxBioreps,
            out var conditionToFiles);
        Assert.That(fileInfoByName.Count, Is.GreaterThan(1), "Need at least two spectra files.");
        Assert.That(conditionToFiles.Count, Is.GreaterThan(1), "Need at least two conditions to compute fold-changes.");
        TestContext.WriteLine($"Files: {fileInfoByName.Count} across conditions " +
            string.Join(", ", conditionToFiles.OrderBy(k => k.Key).Select(k => $"{k.Key}({k.Value.Count})")));

        var sw = Stopwatch.StartNew();
        var load = PipEchoCommon.StreamIdentifications(psmtsv, fileInfoByName, PipEchoCommon.QValueCutoff);
        sw.Stop();
        TestContext.WriteLine($"Loaded {load.Identifications.Count} identifications " +
            $"(kept {load.RowsKept}/{load.RowsRead}; skipped file={load.RowsSkippedNoFile}, q={load.RowsSkippedQ}, " +
            $"ambiguous={load.RowsSkippedAmbiguous}, contaminant={load.RowsSkippedContaminant}) in {sw.Elapsed}.");
        Assert.That(load.Identifications.Count, Is.GreaterThan(0), "No identifications loaded.");

        bool mbr = PipEchoCommon.Env("MZLIB_PIPECHO_MBR", "true").Equals("true", StringComparison.OrdinalIgnoreCase);
        bool normalize = PipEchoCommon.Env("MZLIB_PIPECHO_NORMALIZE", "true").Equals("true", StringComparison.OrdinalIgnoreCase);
        TestContext.WriteLine($"Running FlashLFQ (MBR={mbr}, normalize={normalize}, threads={PipEchoCommon.MaxThreads})...");

        var engine = new FlashLfqEngine(load.Identifications, matchBetweenRuns: mbr, normalize: normalize,
            maxThreads: PipEchoCommon.MaxThreads);
        sw.Restart();
        FlashLfqResults results = engine.Run();
        sw.Stop();
        TestContext.WriteLine($"FlashLFQ finished in {sw.Elapsed}. Quantified {results.ProteinGroups.Count} protein groups.");

        // Mean protein intensity per condition (averaging non-zero replicate intensities).
        var conditionMeans = new Dictionary<string, Dictionary<string, double>>(); // condition -> accession -> mean intensity
        foreach (var (condition, files) in conditionToFiles)
        {
            var perProtein = new Dictionary<string, double>();
            foreach (var (accession, pg) in results.ProteinGroups)
            {
                var vals = files.Select(pg.GetIntensity).Where(v => v > 0).ToList();
                if (vals.Count > 0)
                    perProtein[accession] = vals.Average();
            }
            conditionMeans[condition] = perProtein;
        }

        // Compare every non-reference condition against A (the 1× reference), plus report each pair vs A.
        string reference = ConditionLevel.ContainsKey("A") && conditionMeans.ContainsKey("A")
            ? "A"
            : conditionToFiles.Keys.OrderBy(k => ConditionLevel.GetValueOrDefault(k, double.MaxValue)).First();

        var comparisons = conditionMeans.Keys
            .Where(c => c != reference)
            .OrderBy(c => ConditionLevel.GetValueOrDefault(c, 0))
            .ToList();

        TestContext.WriteLine($"\nReference condition: {reference} (level {ConditionLevel.GetValueOrDefault(reference, double.NaN)}×)\n");
        TestContext.WriteLine($"{"Comparison",-14}{"expEcoliFC",12}{"medEcoliFC",12}{"medHumanFC",12}{"ecoliMAE",10}{"nEcoli",8}{"nHuman",8}{"Ecoli@5%FDP",13}");

        ComparisonResult? headline = null;
        foreach (var cond in comparisons)
        {
            var cr = EvaluateComparison(cond, reference, conditionMeans, load.AccessionSpecies);
            TestContext.WriteLine(
                $"{cond + " vs " + reference,-14}{cr.ExpectedEcoliFc,12:F3}{cr.MedianEcoliFc,12:F3}{cr.MedianHumanFc,12:F3}" +
                $"{cr.EcoliMae,10:F3}{cr.NEcoli,8}{cr.NHuman,8}{cr.EcoliDiscoveredAt5Fdp,13}");
            if (cond == "E" || headline == null)
                headline = cr; // prefer E vs A (largest, 3×) as the headline
        }

        Assert.That(headline, Is.Not.Null);
        WritePerProteinCsv(headline!, psmtsv);

        // Sanity (direction of effect), tolerant so a quant regression is visible without machine-pinning:
        Assert.That(headline!.NEcoli, Is.GreaterThan(0), "No E. coli proteins quantified in the headline comparison.");
        Assert.That(headline.MedianEcoliFc, Is.GreaterThan(headline.MedianHumanFc),
            "E. coli proteins did not rank above human proteins in fold-change — quantification is not recovering the spike.");
        Assert.That(Math.Abs(headline.MedianHumanFc), Is.LessThan(0.5),
            $"Human median log2 FC ({headline.MedianHumanFc:F3}) is far from 0 — normalization/quant is skewed.");
    }

    private sealed record ComparisonResult(
        string Condition, string Reference, double ExpectedEcoliFc,
        double MedianEcoliFc, double MedianHumanFc, double EcoliMae,
        int NEcoli, int NHuman, int EcoliDiscoveredAt5Fdp,
        List<(string accession, Species species, double log2Fc)> PerProtein);

    private static ComparisonResult EvaluateComparison(
        string condition, string reference,
        Dictionary<string, Dictionary<string, double>> conditionMeans,
        Dictionary<string, Species> accessionSpecies)
    {
        var hi = conditionMeans[condition];
        var lo = conditionMeans[reference];
        double expectedEcoliFc = PipEchoCommon.Log2(
            ConditionLevel.GetValueOrDefault(condition, double.NaN) / ConditionLevel.GetValueOrDefault(reference, double.NaN));

        var perProtein = new List<(string accession, Species species, double log2Fc)>();
        foreach (var (accession, hiVal) in hi)
        {
            if (!lo.TryGetValue(accession, out double loVal) || loVal <= 0 || hiVal <= 0)
                continue;
            Species sp = accessionSpecies.GetValueOrDefault(accession, Species.Other);
            if (sp == Species.Other)
                continue;
            perProtein.Add((accession, sp, PipEchoCommon.Log2(hiVal / loVal)));
        }

        var ecoli = perProtein.Where(p => p.species == Species.Ecoli).Select(p => p.log2Fc).ToList();
        var human = perProtein.Where(p => p.species == Species.Human).Select(p => p.log2Fc).ToList();

        double ecoliMae = ecoli.Count > 0 ? ecoli.Average(fc => Math.Abs(fc - expectedEcoliFc)) : double.NaN;
        int discovered = EcoliDiscoveredAtFdp(perProtein, 0.05);

        return new ComparisonResult(condition, reference, expectedEcoliFc,
            PipEchoCommon.Median(ecoli), PipEchoCommon.Median(human), ecoliMae,
            ecoli.Count, human.Count, discovered, perProtein);
    }

    /// <summary>
    /// Ranks proteins by |observed log2 FC| descending (most apparently-changed first) and walks down
    /// counting E. coli as true positives and human as false. Returns the number of E. coli proteins
    /// recovered in the longest prefix whose false-discovery proportion (human / total) stays ≤ fdp.
    /// </summary>
    private static int EcoliDiscoveredAtFdp(List<(string accession, Species species, double log2Fc)> perProtein, double fdp)
    {
        var ordered = perProtein
            .Where(p => p.species is Species.Ecoli or Species.Human)
            .OrderByDescending(p => Math.Abs(p.log2Fc))
            .ToList();

        int ecoli = 0, human = 0, best = 0;
        for (int k = 0; k < ordered.Count; k++)
        {
            if (ordered[k].species == Species.Ecoli) ecoli++; else human++;
            if ((double)human / (ecoli + human) <= fdp)
                best = ecoli;
        }
        return best;
    }

    /// <summary>
    /// Parses the MetaMorpheus ExperimentalDesign.tsv and builds a SpectraFileInfo per run, resolving
    /// each file name to a spectra file on disk. Applies optional condition/biorep subsetting.
    ///
    /// Biological replicates are RE-INDEXED to contiguous 0-based integers within each condition (and
    /// fraction/techrep set to 0, single-fraction data). FlashLFQ's intensity normalization assumes
    /// 0-based, contiguous, gap-free biorep/fraction numbering — passing the design file's 1-based
    /// values (or a non-contiguous subset) makes NormalizeFractions throw on an empty replicate group.
    /// </summary>
    private static Dictionary<string, SpectraFileInfo> BuildFileInfos(
        string designPath, string spectraDir, HashSet<string> conditionFilter, int maxBioreps,
        out Dictionary<string, List<SpectraFileInfo>> conditionToFiles)
    {
        var fileInfoByName = new Dictionary<string, SpectraFileInfo>();
        conditionToFiles = new Dictionary<string, List<SpectraFileInfo>>();

        // Parse rows, keeping the original biorep so we can order replicates deterministically.
        var rows = new List<(string file, string condition, int origBiorep)>();
        foreach (var line in File.ReadAllLines(designPath).Skip(1)) // header: FileName Condition Biorep Fraction Techrep
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            var f = line.Split('\t');
            if (f.Length < 5) continue;

            string condition = f[1].Trim().ToUpperInvariant();
            if (conditionFilter.Count > 0 && !conditionFilter.Contains(condition))
                continue;
            rows.Add((f[0], condition, int.TryParse(f[2], out int b) ? b : 0));
        }

        foreach (var conditionGroup in rows.GroupBy(r => r.condition).OrderBy(g => g.Key))
        {
            string condition = conditionGroup.Key;
            int biorepIndex = 0; // contiguous, 0-based, per condition
            foreach (var row in conditionGroup.OrderBy(r => r.origBiorep).ThenBy(r => r.file))
            {
                if (biorepIndex >= maxBioreps)
                    break;

                string? path = PipEchoCommon.ResolveSpectraFile(spectraDir, row.file);
                if (path == null)
                {
                    TestContext.WriteLine($"  WARNING: no spectra file for '{row.file}' in {spectraDir}; skipped.");
                    continue;
                }

                var info = new SpectraFileInfo(path, condition, biorep: biorepIndex, techrep: 0, fraction: 0);
                fileInfoByName[row.file] = info;
                if (!conditionToFiles.TryGetValue(condition, out var list))
                    conditionToFiles[condition] = list = new List<SpectraFileInfo>();
                list.Add(info);
                biorepIndex++;
            }
        }

        return fileInfoByName;
    }

    private static void WritePerProteinCsv(ComparisonResult cr, string psmtsvPath)
    {
        Directory.CreateDirectory(PipEchoCommon.OutputDir);
        string path = Path.Combine(PipEchoCommon.OutputDir,
            $"SpikeIn_FoldChanges_{cr.Condition}vs{cr.Reference}_{DateTime.Now:yyyyMMdd_HHmmss}.csv");
        var sb = new StringBuilder();
        sb.AppendLine("Accession,Species,ObservedLog2FC,ExpectedLog2FC");
        double expectedHuman = 0;
        foreach (var (accession, species, log2Fc) in cr.PerProtein.OrderByDescending(p => Math.Abs(p.log2Fc)))
        {
            double expected = species == Species.Ecoli ? cr.ExpectedEcoliFc : expectedHuman;
            sb.AppendLine($"{accession},{species},{log2Fc.ToString("G17", CultureInfo.InvariantCulture)},{expected.ToString("G17", CultureInfo.InvariantCulture)}");
        }
        File.WriteAllText(path, sb.ToString());
        TestContext.WriteLine($"\nPer-protein fold changes ({cr.Condition} vs {cr.Reference}) written to: {path}");
    }
}
