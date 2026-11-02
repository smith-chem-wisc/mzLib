using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text.Json;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Development.QuantificationDevelopment.PipEcho
{
    /// <summary>
    /// PIP-ECHO censored-holdout FDP evaluation, driven by a per-dataset config so the single-cell (primary) and
    /// in-house E. coli cases are two [Test]s over one pipeline.
    ///
    /// The FDP computation mirrors the reference analysis at
    /// https://github.com/Alexander-Sol/PIP-ECHOanalysis (Python/FDPAnalysis.ipynb: get_fdp_flash +
    /// get_native_peak_error_rate), "*_dd_5" configuration (donor PEP q-value &lt;= 0.01, MBR PIP q-value cutoff 0.05).
    /// An earlier revision cross-checked a line-for-line Python port and confirmed the C# here is field-for-field
    /// identical, so that port was removed; this C# is now the single source of truth.
    ///
    /// Pipeline (see <see cref="Run"/>):
    ///   0. Resolve one spectra file per PSM File Name. If a "-censored" holdout file is absent (in-house case), the
    ///      uncensored calibrated file is substituted under the censored name via a hardlink (MBR quantifies from MS1,
    ///      which censoring - a removal of MS2 identifications - leaves unchanged). Single-cell has the real censored
    ///      spectra, so nothing is substituted there.
    ///   1. Read AllPSMs.psmtsv (censored search) -> MzLibExtensions.MakeIdentifications.
    ///   2. Read AllPeptides.psmtsv; target peptides with PEP_QValue &lt;= 0.01 -> peptideSequencesToQuantify.
    ///   3. Run FlashLfqEngine with MBR on, donorQValueThreshold 0.01, matchBetweenRunsFdrThreshold (MBR q cutoff) 0.05.
    ///   4. Write QuantifiedPeaks.tsv (== "c_peaks") and QuantifiedPeptides.tsv.
    ///   5. Compute FDP from the written peaks (foreign/entrapment/human counts from THIS run's own human-file peaks;
    ///      native-peak-error true retention times from the reference uncensored run). Reported at the dataset's main
    ///      RT tolerance and a 10 s secondary, each with the corrected AND the reproduced-bug native rate (the
    ///      notebook's discarded, non-inplace sort_values keeps an arbitrary first peak); results written to fdp_csharp.json.
    ///
    /// Reads/writes multi-GB local data, hence [Explicit]. The full dataset is archived on the network share
    /// \\bison.chem.wisc.edu\share\Projects\LFQEvalData; lay it out under D:\PIP_ECHO_PRIDE as the config paths expect.
    /// </summary>
    [TestFixture]
    public static class CensoredFdpEvaluationTest
    {
        // ----- shared FlashLFQ / analysis parameters (identical for every dataset) -----
        private const double DonorPepQCutoff = 0.01;      // AllPeptides PEP_QValue cutoff -> peptideSequencesToQuantify
        private const double DonorQValueThreshold = 0.01; // FlashLFQ donor gate on Identification QValue
        private const double MbrQValueCutoff = 0.05;      // FlashLFQ MBR FDR threshold AND the PIP Q-Value analysis gate
        private const double RtDeltaAltMin = 10.0 / 60.0; // 10 s secondary RT tolerance (both datasets)
        private const string HumanOrganismDefault = "Homo sapiens";

        // scaling_factor from the notebook (Cell 2): 1 + entrapment/human peptide-sequence-count ratio.
        private const double ScalingFactor = 1.0 + 3105275.0 / 3608159.0;

        // Entrapment peptide sequences (scrambled human) appended to the concatenated DBs (shared across datasets).
        private const string EntrapmentSeqPath =
            @"D:\PIP_ECHO_PRIDE\Proteomes\EntrapmentProteinPeptideSequences_50percent.txt";

        private static readonly string[] FdpFieldOrder =
        {
            "human", "foreign", "arabida", "total", "msms_count", "scaled_arabida", "scaled_nper",
            "sensitivity", "error_rate", "fper", "nper", "fder", "false_discovery_proportion"
        };

        /// <summary>Everything that differs between the PIP-ECHO datasets.</summary>
        private sealed record DatasetConfig(
            string Name,
            string AllPsmsPath,
            string AllPeptidesPath,
            string CensoredPsmsPath,
            string ReferencePeaksPath,
            string[] SpectraDirs,
            string StagingDir,
            string HumanFilePattern,
            string ForeignOrganism,
            double RtDeltaMainMin,
            string OutputDir,
            string HumanOrganism = HumanOrganismDefault);

        // Single-cell Human / Yeast (8 human "_1x02nguL_" + 2 "HeYe" mixed). The human "-calib-censored.mzML" holdout
        // spectra ARE present, so nothing is substituted. foreign = S. cerevisiae; notebook kelly_rt = 40 * 0.01 = 0.4.
        private static readonly DatasetConfig SingleCell = new(
            Name: "singlecell-yeast",
            AllPsmsPath: @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell_SearchResults-MetaMorpheus\Task1-SearchTask-Censored\AllPSMs.psmtsv",
            AllPeptidesPath: @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell_SearchResults-MetaMorpheus\Task1-SearchTask-Censored\AllPeptides.psmtsv",
            CensoredPsmsPath: @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell-CensoredFiles-MetaMorpheus\CensoredPsms.psmtsv",
            ReferencePeaksPath: @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell_QuantResults-FlashLFQ_PIP-ECHO\FlashLFQ_7772_DonorPepQ_1\QuantifiedPeaks.tsv",
            SpectraDirs: new[]
            {
                @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell-CensoredFiles-MetaMorpheus",
                @"D:\PIP_ECHO_PRIDE\SingleCellDataset\CalibratedFiles-MetaMorpheus",
            },
            StagingDir: @"D:\PIP_ECHO_PRIDE\SingleCellDataset\SingleCell-CensoredFiles-MetaMorpheus",
            HumanFilePattern: "_1x02nguL_",
            ForeignOrganism: "Saccharomyces cerevisiae",
            RtDeltaMainMin: 0.4,
            OutputDir: @"F:\pipecho-eval\fdp-singlecell");

        // In-house E. coli / Human (10 mixed "Human_Ecoli..._C18" + 10 human "Human_C18"). The human "-calib-censored"
        // holdout spectra are ABSENT, so the resolver substitutes the uncensored "-calib.mzML". The "Human_C18" pattern
        // excludes the mixed files (which read "Human_Ecoli..._C18"). foreign = E. coli; notebook inhouse_rt = 0.6.
        private static readonly DatasetConfig Inhouse = new(
            Name: "inhouse-ecoli",
            AllPsmsPath: @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_SearchResults-MetaMorpheus\MM_ConcatenatedHumanDb_Search_CensoredFiles\Task1-SearchTask\AllPSMs.psmtsv",
            AllPeptidesPath: @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_SearchResults-MetaMorpheus\MM_ConcatenatedHumanDb_Search_CensoredFiles\Task1-SearchTask\AllPeptides.psmtsv",
            CensoredPsmsPath: @"D:\PIP_ECHO_PRIDE\CensoredPsms\Inhouse_MetaMorpheus_CensoredPsms.psmtsv",
            ReferencePeaksPath: @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_QuantResults-FlashLFQ_PIP-ECHO\Human_FlashLFQ_7772_DonorPepQ_1\QuantifiedPeaks.tsv",
            SpectraDirs: new[]
            {
                @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_CensoredFiles-MetaMorpheus",
                @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_CalibratedFiles-MetaMorpheus",
            },
            StagingDir: @"D:\PIP_ECHO_PRIDE\EcoliDataset\Ecoli_CensoredFiles-MetaMorpheus",
            HumanFilePattern: "Human_C18",
            ForeignOrganism: "Escherichia coli",
            RtDeltaMainMin: 0.6,
            OutputDir: @"F:\pipecho-eval\fdp");

        [Test, Explicit("Requires the local single-cell PIP-ECHO dataset on D:.")]
        public static void EvaluateCensoredFdp_SingleCell() => Run(SingleCell);

        [Test, Explicit("Requires the local in-house E. coli PIP-ECHO dataset on D:.")]
        public static void EvaluateCensoredFdp_Inhouse() => Run(Inhouse);

        private static void Run(DatasetConfig cfg)
        {
            Directory.CreateDirectory(cfg.OutputDir);
            TestContext.WriteLine($"=== dataset: {cfg.Name} ===");

            // 1. SpectraFileInfo for every File Name in AllPSMs, all one condition so MBR matches across every file.
            var (spectraFiles, missing) = BuildSpectraFileInfos(cfg);
            if (missing.Count > 0)
                Assert.Fail("Missing source spectra files (FlashLFQ cannot run MBR without them):" +
                            Environment.NewLine + "  " + string.Join(Environment.NewLine + "  ", missing) +
                            Environment.NewLine + "Searched: " + string.Join("; ", cfg.SpectraDirs));
            TestContext.WriteLine($"Resolved {spectraFiles.Count} spectra files.");

            // 2. Identifications from AllPSMs (censored search) via MakeIdentifications.
            var quantFile = new SpectrumMatchFromTsvFile(cfg.AllPsmsPath);
            List<Identification> ids = quantFile.MakeIdentifications(spectraFiles);
            TestContext.WriteLine($"Loaded {ids.Count} identifications " +
                                  $"({ids.Select(i => i.ModifiedSequence).Distinct().Count()} distinct peptides).");

            // 3. Donor peptides: AllPeptides target rows with PEP_QValue <= 0.01.
            List<string> donorSeqs = LoadDonorSequences(cfg.AllPeptidesPath, DonorPepQCutoff);
            TestContext.WriteLine($"Donor peptides (PEP_QValue <= {DonorPepQCutoff}): {donorSeqs.Count}.");

            // 4. Run FlashLFQ with MBR.
            var engine = new FlashLfqEngine(
                allIdentifications: ids,
                peptideSequencesToQuantify: donorSeqs,
                matchBetweenRuns: true,
                matchBetweenRunsFdrThreshold: MbrQValueCutoff,
                donorQValueThreshold: DonorQValueThreshold,
                requireMsmsIdInCondition: false,
                maxThreads: Math.Max(1, Environment.ProcessorCount - 1),
                silent: false);

            FlashLfqResults results = engine.Run();

            string cPeaksPath = Path.Combine(cfg.OutputDir, "QuantifiedPeaks.tsv");
            string peptidesPath = Path.Combine(cfg.OutputDir, "QuantifiedPeptides.tsv");
            results.WriteResults(cPeaksPath, peptidesPath, proteinOutputPath: null, bayesianProteinQuantOutput: null, silent: false);
            TestContext.WriteLine($"Wrote peaks -> {cPeaksPath}");

            // 5. FDP at both RT tolerances, each with the corrected and the reproduced-buggy native rate.
            HashSet<string> entrapment = LoadEntrapmentSet(EntrapmentSeqPath);
            var (mainCorrected, mainBuggy) = ComputeFdp(cfg, cPeaksPath, entrapment, cfg.RtDeltaMainMin);
            var (altCorrected, altBuggy) = ComputeFdp(cfg, cPeaksPath, entrapment, RtDeltaAltMin);

            var csResults = new Dictionary<string, Dictionary<string, double>>
            {
                ["main"] = mainCorrected,
                ["alt"] = altCorrected,
                ["main_buggy"] = mainBuggy,
                ["alt_buggy"] = altBuggy,
            };

            double mainSec = cfg.RtDeltaMainMin * 60;
            PrintFdp($"main, corrected (RT {cfg.RtDeltaMainMin} min / {mainSec:0} s)", mainCorrected);
            PrintFdp($"main, original-bug (RT {cfg.RtDeltaMainMin} min / {mainSec:0} s)", mainBuggy);
            PrintFdp("alt, corrected (RT 10 s)", altCorrected);
            PrintFdp("alt, original-bug (RT 10 s)", altBuggy);

            File.WriteAllText(Path.Combine(cfg.OutputDir, "fdp_csharp.json"),
                JsonSerializer.Serialize(csResults, new JsonSerializerOptions { WriteIndented = true }));

            // 6. Sanity checks: every variant counted peaks and produced finite FDP terms.
            foreach (var (name, r) in csResults)
            {
                Assert.That(r["total"], Is.GreaterThan(0), $"[{name}] no MBR peaks were counted");
                foreach (var k in FdpFieldOrder)
                    Assert.That(double.IsFinite(r[k]), $"[{name}] field '{k}' is not finite ({r[k]})");
            }
            TestContext.WriteLine("FDP computed for all variants; wrote fdp_csharp.json.");
        }

        // ---------- setup helpers ----------

        private static (List<SpectraFileInfo> infos, List<string> missing) BuildSpectraFileInfos(DatasetConfig cfg)
        {
            var names = new List<string>();
            var seen = new HashSet<string>(StringComparer.Ordinal);
            using (var r = new StreamReader(cfg.AllPsmsPath))
            {
                r.ReadLine(); // header
                string? line;
                while ((line = r.ReadLine()) != null)
                {
                    int tab = line.IndexOf('\t');
                    if (tab <= 0) continue;
                    string fn = line.Substring(0, tab);
                    if (seen.Add(fn)) names.Add(fn);
                }
            }
            names.Sort(StringComparer.Ordinal);

            var infos = new List<SpectraFileInfo>();
            var missing = new List<string>();
            int biorep = 0;
            foreach (var name in names)
            {
                string? path = EnsureStagedSpectra(cfg, name);
                if (path == null) { missing.Add(name + ".mzML"); continue; }
                infos.Add(new SpectraFileInfo(path, condition: "proteome", biorep: biorep++, techrep: 0, fraction: 0));
            }
            return (infos, missing);
        }

        private static readonly string[] SpectraExtensions = { ".mzML", ".mzml", ".raw", ".d" };

        /// <summary>
        /// Returns a spectra file named exactly <paramref name="name"/> (matching the PSM File Name, so identifications
        /// and the native-error "-censored" join both line up). If the file exists as-named in any spectra dir it is
        /// used directly; otherwise, for a "-censored" holdout whose real spectrum is absent, the uncensored file is
        /// substituted by hardlinking it under the censored name in StagingDir. Returns null if no source exists.
        /// </summary>
        private static string? EnsureStagedSpectra(DatasetConfig cfg, string name)
        {
            // 1. exact match on disk
            string? direct = Resolve(cfg.SpectraDirs, name);
            if (direct != null) return direct;

            // 2. absent "-censored" holdout -> substitute the uncensored file under the censored name
            if (!name.Contains("-censored", StringComparison.Ordinal)) return null;
            string uncensored = name.Replace("-censored", "");
            string? source = Resolve(cfg.SpectraDirs, uncensored);
            if (source == null) return null;

            Directory.CreateDirectory(cfg.StagingDir);
            string staged = Path.Combine(cfg.StagingDir, name + Path.GetExtension(source));
            if (!File.Exists(staged))
            {
                if (!CreateHardLinkW(staged, source, IntPtr.Zero))
                    File.Copy(source, staged); // fallback if a hardlink can't be made (e.g. across volumes)
            }
            return staged;
        }

        private static string? Resolve(string[] dirs, string nameNoExt)
        {
            foreach (var dir in dirs)
                foreach (var ext in SpectraExtensions)
                {
                    string cand = Path.Combine(dir, nameNoExt + ext);
                    if (File.Exists(cand) || Directory.Exists(cand)) return cand; // .d is a directory
                }
            return null;
        }

        [DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Unicode)]
        [return: MarshalAs(UnmanagedType.Bool)]
        private static extern bool CreateHardLinkW(string lpFileName, string lpExistingFileName, IntPtr lpSecurityAttributes);

        private static List<string> LoadDonorSequences(string allPeptidesPath, double pepQCutoff)
        {
            List<PsmFromTsv> peptides = SpectrumMatchTsvReader.ReadPsmTsv(
                allPeptidesPath, out _, new SpectrumMatchParsingParameters { ParseMatchedFragmentIons = false });
            return peptides
                .Where(p => !IsDecoy(p.DecoyContamTarget) && p.PEP_QValue <= pepQCutoff)
                .Select(p => p.FullSequence)
                .Where(s => !string.IsNullOrEmpty(s) && !s.Contains('|'))
                .Distinct()
                .ToList();
        }

        private static bool IsDecoy(string? decoyContamTarget) =>
            decoyContamTarget != null && decoyContamTarget.Contains('D');

        private static HashSet<string> LoadEntrapmentSet(string path) =>
            new HashSet<string>(File.ReadAllLines(path).Select(l => l.Trim()).Where(l => l.Length > 0), StringComparer.Ordinal);

        // ---------- FDP (mirrors get_fdp_flash, foreign/entrapment/human counted from c_peaks) ----------

        private static (Dictionary<string, double> corrected, Dictionary<string, double> buggy) ComputeFdp(
            DatasetConfig cfg, string cPeaksPath, HashSet<string> entrapment, double rtDeltaMin)
        {
            var (col, rows) = ReadTsv(cPeaksPath);
            int cFile = col["File Name"], cFull = col["Full Sequence"], cBase = col["Base Sequence"],
                cOrg = col["Organism"], cDet = col["Peak Detection Type"], cQ = col["PIP Q-Value"],
                cDecoy = col["Decoy Peptide"], cRand = col["Random RT"];

            var human = rows.Where(r => Get(r, cFile).Contains(cfg.HumanFilePattern, StringComparison.Ordinal)).ToList();

            // Foreign peptides also MS/MS-identified in human files; these MBR "transfers" are legitimate, not errors.
            var msmsForeignSeqs = new HashSet<string>(
                human.Where(r => Get(r, cDet) == "MSMS" && Get(r, cOrg).Contains(cfg.ForeignOrganism, StringComparison.Ordinal))
                     .Select(r => Get(r, cFull)), StringComparer.Ordinal);

            int msmsCount = human.Count(r => Get(r, cDet) == "MSMS" && !ParseBool(Get(r, cDecoy)) &&
                (Get(r, cOrg).Contains(cfg.ForeignOrganism, StringComparison.Ordinal) ||
                 Get(r, cOrg).Contains(cfg.HumanOrganism, StringComparison.Ordinal)));

            var mbr = human.Where(r => Get(r, cDet) == "MBR"
                                    && TryDouble(Get(r, cQ), out double q) && q < MbrQValueCutoff
                                    && !ParseBool(Get(r, cDecoy))
                                    && !ParseBool(Get(r, cRand))).ToList();

            var humanMbr = mbr.Where(r => Get(r, cOrg).Contains(cfg.HumanOrganism, StringComparison.Ordinal)).ToList();
            int arabidaCount = humanMbr.Count(r => entrapment.Contains(Get(r, cFull)));
            int humanCount = humanMbr.Count(r => !entrapment.Contains(Get(r, cBase)));

            var foreignMbr = mbr.Where(r => Get(r, cOrg).Contains(cfg.ForeignOrganism, StringComparison.Ordinal)
                                         && !Get(r, cOrg).Contains(cfg.HumanOrganism, StringComparison.Ordinal)).ToList();
            int foreignCount = foreignMbr.Count - foreignMbr.Count(r => msmsForeignSeqs.Contains(Get(r, cFull)));

            var (sensitivity, correctedErr, buggyErr) =
                ComputeNativeErrorRate(cPeaksPath, cfg.ReferencePeaksPath, cfg.CensoredPsmsPath, rtDeltaMin);

            double scaledArabida = arabidaCount * ScalingFactor;
            double total = Math.Max(1, humanCount + foreignCount + arabidaCount);
            double fper = 100.0 * foreignCount / total;
            double fder = 100.0 * scaledArabida / total;

            // foreign (fper) and false-detection (fder) are identical either way; only the native term (nper) differs.
            Dictionary<string, double> Build(double errorRate)
            {
                double scaledNper = humanCount * errorRate;
                double nper = 100.0 * scaledNper / total;
                return new Dictionary<string, double>(StringComparer.Ordinal)
                {
                    ["human"] = humanCount,
                    ["foreign"] = foreignCount,
                    ["arabida"] = arabidaCount,
                    ["total"] = total,
                    ["msms_count"] = msmsCount,
                    ["scaled_arabida"] = scaledArabida,
                    ["scaled_nper"] = scaledNper,
                    ["sensitivity"] = sensitivity,
                    ["error_rate"] = errorRate,
                    ["fper"] = fper,
                    ["nper"] = nper,
                    ["fder"] = fder,
                    ["false_discovery_proportion"] = fper + nper + fder,
                };
            }

            return (Build(correctedErr), Build(buggyErr));
        }

        // mirrors get_native_peak_error_rate: o_peaks = reference (true RT), c_peaks = this run (MBR recovery).
        // Returns the corrected error rate (per peptide, prefer the RT-agreeing peak) AND the "buggy" error rate that
        // reproduces the notebook's non-inplace sort_values no-op: keep the FIRST peak per peptide in file order.
        private static (double sensitivity, double corrected, double buggy) ComputeNativeErrorRate(
            string cPeaksPath, string oPeaksPath, string censoredPsmPath, double rtDeltaMin)
        {
            // Reference true RT: (file, full) -> apex list (only peaks with a measured apex).
            var (ocol, orows) = ReadTsv(oPeaksPath);
            int oFile = ocol["File Name"], oFull = ocol["Full Sequence"], oApex = ocol["Peak RT Apex"];
            var trueRt = new Dictionary<(string, string), List<double>>();
            foreach (var r in orows)
            {
                string apex = Get(r, oApex);
                if (apex == "-" || !TryDouble(apex, out double rt)) continue;
                var key = (Get(r, oFile), Get(r, oFull));
                if (!trueRt.TryGetValue(key, out var list)) trueRt[key] = list = new List<double>();
                list.Add(rt);
            }

            // Censored PSMs that have a reference peak; join-key file name gets the "-censored" suffix.
            var (ccol, crows) = ReadTsv(censoredPsmPath);
            int ccFile = ccol["File Name"], ccFull = ccol["Full Sequence"];
            var censoredPeaks = new Dictionary<(string, string), List<double>>();
            foreach (var r in crows)
            {
                var refKey = (Get(r, ccFile), Get(r, ccFull));
                if (!trueRt.TryGetValue(refKey, out var trues)) continue;
                var joinKey = (Get(r, ccFile) + "-censored", Get(r, ccFull));
                if (!censoredPeaks.TryGetValue(joinKey, out var list)) censoredPeaks[joinKey] = list = new List<double>();
                list.AddRange(trues);
            }
            int censoredDistinct = censoredPeaks.Count;

            // Our MBR recoveries passing the gate: (file, full) -> list of (apex, apexParsed).
            var (col, rows) = ReadTsv(cPeaksPath);
            int cFile = col["File Name"], cFull = col["Full Sequence"], cDet = col["Peak Detection Type"],
                cQ = col["PIP Q-Value"], cDecoy = col["Decoy Peptide"], cRand = col["Random RT"], cApex = col["Peak RT Apex"];
            var newPeaks = new Dictionary<(string, string), List<(double apex, bool ok)>>();
            foreach (var r in rows)
            {
                if (Get(r, cDet) != "MBR") continue;
                if (ParseBool(Get(r, cRand)) || ParseBool(Get(r, cDecoy))) continue;
                if (!TryDouble(Get(r, cQ), out double q) || !(q < MbrQValueCutoff)) continue;
                bool ok = TryDouble(Get(r, cApex), out double apex);
                var key = (Get(r, cFile), Get(r, cFull));
                if (!newPeaks.TryGetValue(key, out var list)) newPeaks[key] = list = new List<(double, bool)>();
                list.Add((apex, ok));
            }

            // Inner join on (file, seq). Two verdicts per peptide:
            //   corrected: prefer agreement (any RT-agreeing peak counts the peptide as a good transfer);
            //   buggy: keep the FIRST peak in file order (what the notebook's discarded sort_values leaves behind).
            int good = 0, bad = 0, joined = 0;   // corrected
            int goodB = 0, badB = 0;             // buggy
            foreach (var kv in censoredPeaks)
            {
                if (!newPeaks.TryGetValue(kv.Key, out var mbrs)) continue;
                joined++;

                bool anyAgree = false, anyScorable = false;
                foreach (var trueApex in kv.Value)
                    foreach (var (apex, ok) in mbrs)
                    {
                        if (!ok) continue;
                        anyScorable = true;
                        if (Math.Abs(trueApex - apex) < rtDeltaMin) anyAgree = true;
                    }
                if (anyAgree) good++;
                else if (anyScorable) bad++;
                // else: unscorable (Agreement == -1) -> counts toward sensitivity only

                // buggy: first true apex vs first MBR peak, both in file order
                var (firstApex, firstOk) = mbrs[0];
                if (firstOk)
                {
                    if (Math.Abs(kv.Value[0] - firstApex) < rtDeltaMin) goodB++;
                    else badB++;
                }
            }

            double sensitivity = censoredDistinct == 0 ? 0 : (double)joined / censoredDistinct;
            double corrected = (good + bad) == 0 ? 0 : (double)bad / (good + bad);
            double buggy = (goodB + badB) == 0 ? 0 : (double)badB / (goodB + badB);
            return (sensitivity, corrected, buggy);
        }

        // ---------- tiny TSV + parsing helpers ----------

        private static (Dictionary<string, int> col, List<string[]> rows) ReadTsv(string path)
        {
            var col = new Dictionary<string, int>(StringComparer.Ordinal);
            var rows = new List<string[]>();
            using var r = new StreamReader(path);
            string? header = r.ReadLine();
            if (header == null) return (col, rows);
            var h = header.Split('\t');
            for (int i = 0; i < h.Length; i++) col[h[i]] = i;
            string? line;
            while ((line = r.ReadLine()) != null)
            {
                if (line.Length == 0) continue;
                rows.Add(line.Split('\t'));
            }
            return (col, rows);
        }

        private static string Get(string[] row, int idx) => idx >= 0 && idx < row.Length ? row[idx] : "";
        private static bool ParseBool(string s) => string.Equals(s, "True", StringComparison.OrdinalIgnoreCase);
        private static bool TryDouble(string s, out double v) =>
            double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v);

        private static void PrintFdp(string label, Dictionary<string, double> f)
        {
            TestContext.WriteLine($"--- {label} ---");
            foreach (var k in FdpFieldOrder)
                TestContext.WriteLine($"  {k,-28} {f[k].ToString("0.######", CultureInfo.InvariantCulture)}");
        }
    }
}
