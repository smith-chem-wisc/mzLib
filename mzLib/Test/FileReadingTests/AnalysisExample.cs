using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Plotly.NET.CSharp;
using Proteomics.AminoAcidPolymer;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Serialization;
using Test.FileReadingTests;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.FileReadingTests
{
    [TestFixture]
    internal class AnalysisExample
    {
        private enum SimulationRunMode
        {
            QuickDev,
            FullFidelity,
        }

        private sealed record SimulationRunProfile(
            SimulationRunMode Mode,
            int MaxRecords,
            double RtHalfWidth,
            double IntensityThresholdFraction,
            double MinIntensityThreshold);

        private static SimulationRunProfile GetSimulationRunProfile()
        {
            string? mode = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_MODE");
            if (string.Equals(mode, "full", StringComparison.OrdinalIgnoreCase)
                || string.Equals(mode, "fullfidelity", StringComparison.OrdinalIgnoreCase))
            {
                return new SimulationRunProfile(
                    Mode: SimulationRunMode.FullFidelity,
                    MaxRecords: 200,
                    RtHalfWidth: 0.40,
                    IntensityThresholdFraction: 1e-5,
                    MinIntensityThreshold: 0.1);
            }

            return new SimulationRunProfile(
                Mode: SimulationRunMode.QuickDev,
                MaxRecords: 25,
                RtHalfWidth: 0.25,
                IntensityThresholdFraction: 1e-4,
                MinIntensityThreshold: 1.0);
        }

        private static bool GetDeduplicateProteoforms()
        {
            if (IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_NO_DEDUP"))
                return false;

            return true;
        }

        private static bool GetGlobalAbundanceRefitEnabled()
        {
            if (IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_NO_GLOBAL_ABUNDANCE_REFIT"))
                return false;

            return true;
        }

        private static int GetGlobalAbundanceRefitMaxModels()
        {
            const int defaultMaxModels = 200;
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_GLOBAL_REFIT_MAX_MODELS");
            if (string.IsNullOrWhiteSpace(raw))
                return defaultMaxModels;

            return int.TryParse(raw, out int parsed) && parsed > 0
                ? parsed
                : defaultMaxModels;
        }

        private static bool IsTrueEnvironmentVariable(string variableName)
        {
            var value = Environment.GetEnvironmentVariable(variableName);
            return string.Equals(value, "1", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "true", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "yes", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "on", StringComparison.OrdinalIgnoreCase);
        }


        [Test]
        [Explicit("Interactive Plotly demo for a synthetic MS1 spectrum")]
        public void PlotSyntheticSpectrum()
        {
            var simulation = BuildSimulation();
            var scan = simulation.Scans[simulation.Scans.Length / 2];

            Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Synthetic MS1 Spectrum")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        [Explicit("Interactive Plotly demo for a synthetic charge XIC")]
        public void PlotSyntheticChargeXic()
        {
            const double mass = 10000.0;
            var simulation = BuildSimulation();
            var extractor = new GroundTruthExtractor(simulation.Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            var truth = extractor.Extract(mass, rtCenter: 20.0, rtHalfWidth: 1.0, minCharge: 6, maxCharge: 11);

            int chargeOffset = 8 - truth.MinCharge;
            Chart.Line<double, double, string>(truth.ScanTimes, truth.ChargeXics[chargeOffset])
                .WithTraceInfo("Charge 8 XIC")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        private static SimulationResult BuildSimulation()
        {
            var model = new ProteoformModel(
                MonoisotopicMass: 10000.0,
                Abundance: 1.5e6,
                RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

            double[] scanTimes = Enumerable.Range(0, 31).Select(i => 18.5 + i * 0.1).ToArray();
            return new Simulator().Simulate(new[] { model }, minCharge: 6, maxCharge: 11, sigmaMz: 0.012, scanTimes: scanTimes);
        }

        [Test]
        public static void ControlXIC()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder > 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];

            for (int i = 0; i < ms2Scans.Count; i++)
            {
                rtArray[i] = ms2Scans[i].RetentionTime;
                intensityArray[i] = 0;
                int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
                if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

                int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
                if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

                intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            }

            Chart.Line<double, double, string>(rtArray, intensityArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity of Diagnostic Ions"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        public static void RawDataLoadingTimer()
        {
            var path = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            long size = new FileInfo(path).Length;
            var sizeMb = size / (1024.0 * 1024.0);
            var sw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            sw.Stop();
            Console.WriteLine($"Loaded {sizeMb:F2} MB .raw file in {sw.Elapsed.ToString(@"hh\:mm\:ss\.fff")}");

            sw = Stopwatch.StartNew();
            var ms1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            sw.Stop();
            Console.WriteLine($"Collected {ms1Scans.Length} MS1 scans in {sw.Elapsed.ToString(@"hh\:mm\:ss\.fff")}");
        }

        [Test]
        public static void MzmlDataLoadingTimer()
        {
            var path = @"D:\JurkatTopdown\02-17-20_jurkat_td_rep2_fract2-calib-averaged.mzML";
            long size = new FileInfo(path).Length;
            var sizeMb = size / (1024.0 * 1024.0);
            var sw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            sw.Stop();
            Console.WriteLine($"Loaded {sizeMb:F2} MB .mzML file in {sw.Elapsed}");
        }

        [Test]
        public static void GetIsoEnv()
        {
            string histoneSeq = "MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK";

            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(histoneSeq).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
            double[] mz = dist.Masses.Select(v => v.ToMz(20)).ToArray();
            double[] intensities = dist.Intensities.Select(v => v * 100).ToArray();
            double rt = 1;

            ChemicalFormula methyl = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("CH2"), cf } );
            ChemicalFormula acetyl = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("C2H2O"), cf });
            ChemicalFormula phospho = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("PO4H3"), cf });




            double[] plotMz = new double[mz.Length * 3];
            double[] plotIntensity = new double[mz.Length * 3];

            for (int i = 0; i < mz.Length; i++)
            {
                int first = (3 * i);
                int second = (3 * i + 1);
                int third = (3 * i + 2);

                plotMz[first] = mz[i] - 0.00001;
                plotIntensity[first] = 0;

                plotMz[second] = mz[i];
                plotIntensity[second] = intensities[i];

                plotMz[third] = mz[i] + 0.00001;
                plotIntensity[third] = 0;
            }
            

            // add the scan
            MsDataScan scan = new MsDataScan(massSpectrum: new MzSpectrum(plotMz, plotIntensity, false), oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1));

            Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Theoretical Envelope")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();





        }

        [Test]
        public static void Ms1Example()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder == 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];


            var scan = ms2Scans[30 + ms2Scans.Count / 2];

            //for (int i = 0; i < ms2Scans.Count; i++)
            //{
            //    rtArray[i] = ms2Scans[i].RetentionTime;
            //    intensityArray[i] = 0;
            //    int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

            //    int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

            //    intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            //}

            Chart.Point<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        public static void ConvertJurkatTopdownRawToMzml()
        {
            string rawPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            string outPath = Path.ChangeExtension(rawPath, ".mzML");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, outPath, false);
        }

        [Test]
        public static void Ms2LogExample()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder == 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];


            var scan = ms2Scans[30 + ms2Scans.Count / 2];

            //for (int i = 0; i < ms2Scans.Count; i++)
            //{
            //    rtArray[i] = ms2Scans[i].RetentionTime;
            //    intensityArray[i] = 0;
            //    int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

            //    int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

            //    intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            //}

            Chart.Bar<double, double, string>(scan.MassSpectrum.XArray.Select(x => Math.Log(x)), scan.MassSpectrum.YArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        [Explicit("Fits TopDownSimulator models from rep2 and writes centroided simulated mzML")]
        public static void SimulateJurkatRep2AndWriteMzml()
        {
            var totalSw = Stopwatch.StartNew();
            var profile = GetSimulationRunProfile();
            bool deduplicate = GetDeduplicateProteoforms();

            Console.WriteLine($"Simulation mode: {profile.Mode} (set MZLIB_TOPDOWN_SIM_MODE=full for full-fidelity)");
            Console.WriteLine($"Deduplicate proteoforms: {deduplicate} (set MZLIB_TOPDOWN_SIM_NO_DEDUP=1 to disable)");

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string mmPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\AllProteoforms.psmtsv");

            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;
            string mzmlOutPath = Path.Combine(outDir, stem + ".simulated.mzML");

            var stageSw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            Console.WriteLine($"Loaded source file in {stageSw.Elapsed}");

            stageSw.Restart();
            var ms1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Console.WriteLine($"Collected {ms1Scans.Length} MS1 scans in {stageSw.Elapsed}");

            stageSw.Restart();
            var extractor = new GroundTruthExtractor(ms1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            Console.WriteLine($"Built ground-truth extractor in {stageSw.Elapsed}");

            stageSw.Restart();
            var loader = new MmResultLoader();
            var records = loader.Load(mmPath)
                .Where(r => string.Equals(r.FileNameWithoutExtension, stem, StringComparison.OrdinalIgnoreCase))
                .OrderByDescending(r => r.Score)
                .Take(profile.MaxRecords)
                .ToArray();

            var simulationRecords = deduplicate ? DeduplicateByProteoform(records) : records;
            Console.WriteLine($"Loaded and filtered {records.Length} MM records in {stageSw.Elapsed}");
            if (simulationRecords.Length != records.Length)
                Console.WriteLine($"Deduplicated records: {simulationRecords.Length}");

            Assert.That(simulationRecords, Is.Not.Empty, "No MetaMorpheus proteoform rows matched the rep2 raw filename.");

            var fitter = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(fallbackSigmaMz: 0.012));
            var fitted = new List<FittedProteoform>(simulationRecords.Length);

            int globalMinCharge = int.MaxValue;
            int globalMaxCharge = int.MinValue;

            int counter = 0;
            foreach (var record in simulationRecords)
            {
                int minCharge = Math.Max(2, record.PrecursorCharge - 2);
                int maxCharge = Math.Min(80, record.PrecursorCharge + 2);
                if (minCharge > maxCharge)
                    continue;

                Console.WriteLine($"Starting fit {counter}/{simulationRecords.Length}: {record.Identifier}, charge range {minCharge}-{maxCharge}");

                var truth = extractor.Extract(record.MonoisotopicMass, record.RetentionTime, rtHalfWidth: profile.RtHalfWidth, minCharge: minCharge, maxCharge: maxCharge);
                FittedProteoform fit;
                try
                {
                    fit = fitter.Fit(truth, record.Identifier);
                }
                catch (InvalidOperationException)
                {
                    continue;
                }

                if (double.IsNaN(fit.Model.Abundance) || fit.Model.Abundance <= 0)
                    continue;

                fitted.Add(fit);
                globalMinCharge = Math.Min(globalMinCharge, minCharge);
                globalMaxCharge = Math.Max(globalMaxCharge, maxCharge);
                counter++;
                Console.WriteLine($"Completed fit {counter}/{simulationRecords.Length}: {record.Identifier}, charge range {minCharge}-{maxCharge}, fitted abundance {fit.Model.Abundance:F2}, sigmaMz {fit.SigmaMz:F6}");
            }

            Assert.That(fitted, Is.Not.Empty, "No fit records were produced from rep2 MM results.");

            var sigmaCandidates = fitted
                .Select(f => f.SigmaMz)
                .Where(s => !double.IsNaN(s) && !double.IsInfinity(s) && s > 0)
                .OrderBy(s => s)
                .ToArray();
            double sigmaMz = sigmaCandidates.Length == 0 ? 0.012 : sigmaCandidates[sigmaCandidates.Length / 2];

            if (globalMinCharge == int.MaxValue)
                globalMinCharge = 2;
            if (globalMaxCharge == int.MinValue)
                globalMaxCharge = 80;

            var scanTimes = ms1Scans.Select(s => s.RetentionTime).ToArray();
            var models = fitted.Select(f => f.Model).ToArray();

            stageSw.Restart();
            var export = new Simulator().WriteMzml(
                models,
                globalMinCharge,
                globalMaxCharge,
                sigmaMz,
                scanTimes,
                mzmlOutPath,
                new ScanReductionOptions
                {
                    RelativeIntensityThreshold = profile.IntensityThresholdFraction,
                    MinimumIntensity = profile.MinIntensityThreshold,
                });
            Console.WriteLine($"Simulated and wrote mzML in {stageSw.Elapsed}");

            Console.WriteLine($"Simulated models: {models.Length}");
            Console.WriteLine($"Fitted sigmaMz (median): {sigmaMz:F6}");
            Console.WriteLine($"Charge range: {globalMinCharge}-{globalMaxCharge}");
            Console.WriteLine($"Wrote simulated mzML: {export.MzmlPath}");
            Console.WriteLine($"Ground truth: {export.GroundTruthPath}");
            Console.WriteLine($"Simulated scans: {export.ScanCount}, peak count: {export.PeakCount}");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        [Test]
        [Explicit("Writes centroided simulated mzML for the 31-35 min slice and the full q<=0.01 run of rep2 fract7")]
        public static void ExportRep2SliceAndQValueSimulations()
        {
            const double rtStart = 31.0;
            const double rtEnd = 35.0;
            const double qValueThreshold = 0.01;
            const double rtHalfWidth = 0.25;
            bool deduplicate = GetDeduplicateProteoforms();
            bool useGlobalAbundanceRefit = GetGlobalAbundanceRefitEnabled();
            int globalRefitMaxModels = GetGlobalAbundanceRefitMaxModels();

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string resultPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv");
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;

            string simSliceMzmlPath = Path.Combine(outDir, stem + ".rt31-35.simulated.q001.mzML");
            string simFullMzmlPath = Path.Combine(outDir, stem + ".full.simulated.q001.mzML");

            var reduction = new ScanReductionOptions
            {
                RelativeIntensityThreshold = 1e-4,
                MinimumIntensity = 1.0,
            };

            var totalSw = Stopwatch.StartNew();
            Console.WriteLine($"Global abundance refit: {useGlobalAbundanceRefit} (set MZLIB_TOPDOWN_SIM_NO_GLOBAL_ABUNDANCE_REFIT=1 to disable)");
            Console.WriteLine($"Global abundance refit max models: {globalRefitMaxModels} (override with MZLIB_TOPDOWN_SIM_GLOBAL_REFIT_MAX_MODELS)");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var allMs1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Assert.That(allMs1Scans, Is.Not.Empty, "No MS1 scans were found in the raw file.");

            var rtSliceMs1Scans = allMs1Scans
                .Where(s => s.RetentionTime >= rtStart && s.RetentionTime <= rtEnd)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Assert.That(rtSliceMs1Scans, Is.Not.Empty, "No MS1 scans were found in the 31-35 min window.");

            // Real data is not re-exported here; ConvertJurkatTopdownRawToMzml handles raw -> mzML.
            // Its centroid status is logged because mzLib cannot read a profile-mode mzML back.
            Console.WriteLine($"Real slice scans: {rtSliceMs1Scans.Length} " +
                              $"(centroided: {rtSliceMs1Scans.Count(s => s.IsCentroid)}/{rtSliceMs1Scans.Length})");

            var qFilteredRecords = LoadQualifiedMmRecords(resultPath, stem, qValueThreshold, rtStart: null, rtEnd: null);
            if (deduplicate)
            {
                qFilteredRecords = DeduplicateByProteoform(qFilteredRecords);
                Console.WriteLine($"Deduplicated q<=0.01 record count: {qFilteredRecords.Length}");
            }
            var qFilteredSliceRecords = qFilteredRecords
                .Where(r => r.RetentionTime >= rtStart && r.RetentionTime <= rtEnd)
                .ToArray();

            Assert.That(qFilteredSliceRecords, Is.Not.Empty, "No proteoforms passed q<=0.01 inside 31-35 min.");
            Assert.That(qFilteredRecords, Is.Not.Empty, "No proteoforms passed q<=0.01 for full run.");

            var extractor = new GroundTruthExtractor(allMs1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

            // The slice records are a subset of the full set, so everything is extracted and fitted
            // once and the slice is then selected out of it. The global refit is a joint fit, so it
            // is still applied separately to each output set — the slice refit is conditioned on the
            // slice, exactly as before.
            var stageSw = Stopwatch.StartNew();
            var allFits = FitProteoforms(qFilteredRecords, extractor, rtHalfWidth);
            Console.WriteLine($"Fitted {allFits.Fits.Length}/{qFilteredRecords.Length} records in {stageSw.Elapsed}");
            Assert.That(allFits.Fits, Is.Not.Empty, "No simulated models were fitted for full q<=0.01 run.");

            var sliceSelection = Enumerable.Range(0, allFits.Fits.Length)
                .Where(i => allFits.Records[i].RetentionTime >= rtStart && allFits.Records[i].RetentionTime <= rtEnd)
                .ToArray();
            Assert.That(sliceSelection, Is.Not.Empty, "No fitted proteoforms fell inside 31-35 min.");

            var sliceModels = ApplyGlobalRefit(
                sliceSelection.Select(i => allFits.Fits[i]).ToArray(),
                sliceSelection.Select(i => allFits.Truths[i]).ToArray(),
                allFits.MinCharge, allFits.MaxCharge, allFits.SigmaMz,
                useGlobalAbundanceRefit, globalRefitMaxModels, "31-35 min slice");

            var sliceScanTimes = rtSliceMs1Scans.Select(s => s.RetentionTime).ToArray();
            var sliceExport = new Simulator().WriteMzml(
                sliceModels,
                allFits.MinCharge,
                allFits.MaxCharge,
                allFits.SigmaMz,
                sliceScanTimes,
                simSliceMzmlPath,
                reduction);

            Console.WriteLine($"Simulated slice mzML written: {sliceExport.MzmlPath}");
            Console.WriteLine($"Ground truth: {sliceExport.GroundTruthPath}");
            Console.WriteLine($"Simulated slice models: {sliceModels.Length}, scans: {sliceExport.ScanCount}, peaks: {sliceExport.PeakCount}");

            var fullModels = ApplyGlobalRefit(
                allFits.Fits, allFits.Truths,
                allFits.MinCharge, allFits.MaxCharge, allFits.SigmaMz,
                useGlobalAbundanceRefit, globalRefitMaxModels, "full q<=0.01 run");

            var fullScanTimes = allMs1Scans.Select(s => s.RetentionTime).ToArray();
            var fullExport = new Simulator().WriteMzml(
                fullModels,
                allFits.MinCharge,
                allFits.MaxCharge,
                allFits.SigmaMz,
                fullScanTimes,
                simFullMzmlPath,
                reduction);

            Console.WriteLine($"Simulated full mzML written: {fullExport.MzmlPath}");
            Console.WriteLine($"Ground truth: {fullExport.GroundTruthPath}");
            Console.WriteLine($"Simulated full models: {fullModels.Length}, scans: {fullExport.ScanCount}, peaks: {fullExport.PeakCount}");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        private static MmResultRecord[] LoadQualifiedMmRecords(
            string psmTsvPath,
            string expectedFileStem,
            double qValueThreshold,
            double? rtStart,
            double? rtEnd)
        {
            var file = new PsmFromTsvFile(psmTsvPath, new SpectrumMatchParsingParameters
            {
                ParseMatchedFragmentIons = false,
            });
            file.LoadResults();

            var records = file.Results
                .Where(p => p is not null)
                .Where(p => string.Equals(p.FileNameWithoutExtension, expectedFileStem, StringComparison.OrdinalIgnoreCase))
                .Where(p => p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
                .Where(p => !double.IsNaN(p.QValue) && p.QValue <= qValueThreshold)
                .Where(p => !rtStart.HasValue || p.RetentionTime >= rtStart.Value)
                .Where(p => !rtEnd.HasValue || p.RetentionTime <= rtEnd.Value)
                .Select(p => new MmResultRecord(
                    FileNameWithoutExtension: p.FileNameWithoutExtension,
                    PrecursorScanNumber: p.PrecursorScanNum,
                    Ms2ScanNumber: p.Ms2ScanNumber,
                    PrecursorCharge: p.PrecursorCharge,
                    MonoisotopicMass: p.MonoisotopicMass,
                    RetentionTime: p.RetentionTime,
                    Score: p.Score,
                    FullSequence: p.FullSequence,
                    Accession: p.Accession,
                    Identifier: BuildIdentifier(p)))
                .OrderByDescending(r => r.Score)
                .ThenBy(r => r.RetentionTime)
                .ToArray();

            Console.WriteLine($"Loaded qualified records from {psmTsvPath}: {records.Length} (q<={qValueThreshold})");
            return records;
        }

        private sealed record FitBatch(
            FittedProteoform[] Fits,
            ProteoformGroundTruth[] Truths,
            MmResultRecord[] Records,
            int MinCharge,
            int MaxCharge,
            double SigmaMz);

        /// <summary>
        /// Extracts and fits every record. No global refit is applied here — that is a joint fit
        /// over a chosen set of models, so it belongs to whoever assembles that set.
        /// </summary>
        /// <remarks>
        /// The record loop runs in parallel. Each iteration builds its own IsotopeEnvelopeKernel
        /// and its own fitters, and GroundTruthExtractor is read-only after construction, so
        /// nothing is shared across threads. Results are written into preallocated slots and then
        /// compacted in record order, so the surviving set, the model ordering and the median sigma
        /// are all independent of how the work happened to be scheduled.
        /// </remarks>
        private static FitBatch FitProteoforms(
            IReadOnlyList<MmResultRecord> records,
            GroundTruthExtractor extractor,
            double rtHalfWidth)
        {
            var fits = new FittedProteoform[records.Count];
            var truths = new ProteoformGroundTruth[records.Count];
            var chargeRanges = new (int Min, int Max)[records.Count];
            int completed = 0;

            Parallel.For(0, records.Count, i =>
            {
                var record = records[i];
                int minCharge = Math.Max(2, record.PrecursorCharge - 2);
                int maxCharge = Math.Min(80, record.PrecursorCharge + 2);
                if (minCharge > maxCharge)
                    return;

                var truth = extractor.Extract(record.MonoisotopicMass, record.RetentionTime, rtHalfWidth, minCharge, maxCharge);

                FittedProteoform fit;
                try
                {
                    fit = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(fallbackSigmaMz: 0.012))
                        .Fit(truth, record.Identifier);
                }
                catch (InvalidOperationException)
                {
                    return;
                }

                if (double.IsNaN(fit.Model.Abundance) || fit.Model.Abundance <= 0)
                    return;

                fits[i] = fit;
                truths[i] = truth;
                chargeRanges[i] = (minCharge, maxCharge);

                int done = Interlocked.Increment(ref completed);
                if (done % 100 == 0)
                    Console.WriteLine($"Fitted {done}/{records.Count} records");
            });

            var keptFits = new List<FittedProteoform>(records.Count);
            var keptTruths = new List<ProteoformGroundTruth>(records.Count);
            var keptRecords = new List<MmResultRecord>(records.Count);

            int minChargeGlobal = int.MaxValue;
            int maxChargeGlobal = int.MinValue;

            for (int i = 0; i < records.Count; i++)
            {
                if (fits[i] is null)
                    continue;

                keptFits.Add(fits[i]);
                keptTruths.Add(truths[i]);
                keptRecords.Add(records[i]);
                minChargeGlobal = Math.Min(minChargeGlobal, chargeRanges[i].Min);
                maxChargeGlobal = Math.Max(maxChargeGlobal, chargeRanges[i].Max);
            }

            if (minChargeGlobal == int.MaxValue)
                minChargeGlobal = 2;
            if (maxChargeGlobal == int.MinValue)
                maxChargeGlobal = 80;

            var sigmaCandidates = keptFits
                .Select(f => f.SigmaMz)
                .Where(s => !double.IsNaN(s) && !double.IsInfinity(s) && s > 0)
                .OrderBy(s => s)
                .ToArray();
            double sigmaMz = sigmaCandidates.Length == 0 ? 0.012 : sigmaCandidates[sigmaCandidates.Length / 2];

            return new FitBatch(
                keptFits.ToArray(),
                keptTruths.ToArray(),
                keptRecords.ToArray(),
                minChargeGlobal,
                maxChargeGlobal,
                sigmaMz);
        }

        /// <summary>
        /// Runs the joint non-negative abundance refit over one set of fitted proteoforms and
        /// returns their models. The refit is conditioned on exactly the set passed in.
        /// </summary>
        private static ProteoformModel[] ApplyGlobalRefit(
            FittedProteoform[] fits,
            ProteoformGroundTruth[] truths,
            int minCharge,
            int maxCharge,
            double sigmaMz,
            bool useGlobalAbundanceRefit,
            int globalRefitMaxModels,
            string label)
        {
            if (!useGlobalAbundanceRefit || fits.Length <= 1)
                return fits.Select(f => f.Model).ToArray();

            if (fits.Length > globalRefitMaxModels)
            {
                Console.WriteLine($"Skipping global abundance refit for {label}: {fits.Length} models exceeds max {globalRefitMaxModels}.");
                return fits.Select(f => f.Model).ToArray();
            }

            var sw = Stopwatch.StartNew();
            var refitter = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(
                MaxIterations: 8,
                ConvergenceTolerance: 1e-3,
                MinimumAbundance: 0,
                Verbose: true));

            var refitResult = refitter.Refit(fits, truths, minCharge, maxCharge, sigmaMz);
            Console.WriteLine($"Global abundance refit ({label}) over {fits.Length} models took {sw.Elapsed}");
            Console.WriteLine($"  iterations: {refitResult.IterationsCompleted}, converged: {refitResult.Converged}");
            Console.WriteLine($"  residual fraction: {refitResult.InitialResidualFraction:G6} -> {refitResult.FinalResidualFraction:G6}");

            return refitResult.FittedProteoforms.Select(f => f.Model).ToArray();
        }

        private static MmResultRecord[] DeduplicateByProteoform(IReadOnlyList<MmResultRecord> records)
        {
            return records
                .GroupBy(r => (
                    Accession: string.IsNullOrWhiteSpace(r.Accession) ? string.Empty : r.Accession,
                    Sequence: r.FullSequence,
                    Charge: r.PrecursorCharge))
                .Select(g => g
                    .OrderByDescending(r => r.Score)
                    .ThenBy(r => r.RetentionTime)
                    .First())
                .OrderByDescending(r => r.Score)
                .ThenBy(r => r.RetentionTime)
                .ToArray();
        }

        private static string BuildIdentifier(PsmFromTsv psm)
        {
            if (!string.IsNullOrWhiteSpace(psm.Accession))
                return $"{psm.Accession}:{psm.FullSequence}:{psm.Ms2ScanNumber}";

            return $"{psm.FileNameWithoutExtension}:{psm.Ms2ScanNumber}";
        }


        private static string ResolveLocalPath(string preferredPath)
        {
            if (File.Exists(preferredPath) || Directory.Exists(preferredPath))
                return preferredPath;

            if (preferredPath.Length >= 3 && preferredPath[1] == ':' && (preferredPath[2] == '\\' || preferredPath[2] == '/'))
            {
                char drive = char.ToLowerInvariant(preferredPath[0]);
                string remainder = preferredPath.Substring(3).Replace('\\', '/');
                string wslPath = $"/mnt/{drive}/{remainder}";
                if (File.Exists(wslPath) || Directory.Exists(wslPath))
                    return wslPath;
            }

            return preferredPath;
        }

    }
}
