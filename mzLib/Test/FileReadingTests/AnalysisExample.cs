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
            int PointsPerSigma,
            double MzPaddingInSigmas,
            double ImspThresholdFraction,
            double MinImspThreshold);

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
                    PointsPerSigma: 3,
                    MzPaddingInSigmas: 6.0,
                    ImspThresholdFraction: 1e-5,
                    MinImspThreshold: 0.1);
            }

            return new SimulationRunProfile(
                Mode: SimulationRunMode.QuickDev,
                MaxRecords: 25,
                RtHalfWidth: 0.25,
                PointsPerSigma: 1,
                MzPaddingInSigmas: 4.0,
                ImspThresholdFraction: 1e-4,
                MinImspThreshold: 1.0);
        }

        private static bool GetSimulateCentroidedOutput()
        {
            return IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_CENTROID");
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
            var index = PeakIndexingEngine.InitializeIndexingEngine(simulation.Scans)!;
            var extractor = new GroundTruthExtractor(index, simulation.Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
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

            var outDir = Path.GetDirectoryName(path);
            var outPath = Path.Combine(outDir!, Path.GetFileNameWithoutExtension(path) + "_test.imsp");

            int peakCount = WriteImspFile(ms1Scans, outPath, intensityThreshold: 10000);
            sw.Stop();
            Console.WriteLine($"Wrote {peakCount} peaks in {sw.Elapsed.ToString(@"hh\:mm\:ss\.fff")}");


            if(File.Exists(outPath))
                File.Delete(outPath);
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

        /// <summary>
        /// Reads the Jurkat top-down raw file, bins every MS1 peak by m/z (100 bins per Dalton,
        /// matching the FlashLFQ PeakIndexingEngine strategy), and writes the result as a
        /// portable binary .imsp file plus a human-readable format specification (.md).
        ///
        /// File layout overview – see IMSP_Format.md for the full spec:
        ///   Section 1 – File header         (24 bytes)
        ///   Section 2 – Scan table          (ScanCount × 16 bytes)
        ///   Section 3 – Bin directory       (NonEmptyBinCount × 12 bytes)
        ///   Section 4 – Peak array          (TotalPeakCount × 12 bytes)
        /// </summary>
        [Test]
        public static void IndexJurkatTopdownMs1Peaks()
        {
            string rawPath  = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            string outDir   = Path.GetDirectoryName(rawPath);
            string stem     = Path.GetFileNameWithoutExtension(rawPath);
            string imspPath = Path.Combine(outDir, stem + ".imsp");
            string mdPath   = Path.Combine(outDir, "IMSP_Format.md");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var ms1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            int peakCount = WriteImspFile(ms1Scans, imspPath, intensityThreshold: 10000);
            Console.WriteLine($"Total MS1 peaks written: {peakCount}");

            File.WriteAllText(mdPath, ImspFormatMarkdown);
        }

        [Test]
        [Explicit("Fits TopDownSimulator models from rep2 and writes simulated IMSP")]
        public static void SimulateJurkatRep2AndWriteImsp()
        {
            var totalSw = Stopwatch.StartNew();
            var profile = GetSimulationRunProfile();
            bool centroidOutput = GetSimulateCentroidedOutput();
            bool deduplicate = GetDeduplicateProteoforms();

            Console.WriteLine($"Simulation mode: {profile.Mode} (set MZLIB_TOPDOWN_SIM_MODE=full for full-fidelity)");
            Console.WriteLine($"Centroid output: {centroidOutput} (set MZLIB_TOPDOWN_SIM_CENTROID=1 to enable)");
            Console.WriteLine($"Deduplicate proteoforms: {deduplicate} (set MZLIB_TOPDOWN_SIM_NO_DEDUP=1 to disable)");

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string mmPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\AllProteoforms.psmtsv");

            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;
            string imspOutPath = Path.Combine(outDir, stem + ".simulated.imsp");

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
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(ms1Scans)!;
            var extractor = new GroundTruthExtractor(indexingEngine, ms1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            Console.WriteLine($"Built indexing engine in {stageSw.Elapsed}");

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
            var simulation = new Simulator().Simulate(
                models,
                globalMinCharge,
                globalMaxCharge,
                sigmaMz,
                scanTimes,
                pointsPerSigma: profile.PointsPerSigma,
                mzPaddingInSigmas: profile.MzPaddingInSigmas);
            Console.WriteLine($"Simulated {simulation.Scans.Length} scans in {stageSw.Elapsed}");

            var scansForExport = centroidOutput
                ? CentroidizeSimulatedScans(simulation.Scans, relativeIntensityThreshold: 1e-4)
                : simulation.Scans;

            stageSw.Restart();
            double maxSimIntensity = 0;
            foreach (var scan in scansForExport)
            {
                if (!scan.MassSpectrum.YArray.Any())
                    continue;

                double localMax = scan.MassSpectrum.YArray.Max();
                if (localMax > maxSimIntensity)
                    maxSimIntensity = localMax;
            }

            double intensityThreshold = Math.Max(profile.MinImspThreshold, maxSimIntensity * profile.ImspThresholdFraction);
            int peakCount = WriteImspFile(scansForExport, imspOutPath, intensityThreshold: intensityThreshold);
            Console.WriteLine($"Wrote IMSP in {stageSw.Elapsed} (threshold={intensityThreshold:F2})");

            Console.WriteLine($"Simulated models: {models.Length}");
            Console.WriteLine($"Fitted sigmaMz (median): {sigmaMz:F6}");
            Console.WriteLine($"Charge range: {globalMinCharge}-{globalMaxCharge}");
            Console.WriteLine($"Wrote simulated IMSP: {imspOutPath}");
            Console.WriteLine($"Simulated peak count: {peakCount}");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        [Test]
        [Explicit("Writes real 31-35 min IMSP slice, simulated 31-35 slice, and full q<=0.01 simulation for rep2 fract7")]
        public static void ExportRep2SliceAndQValueSimulations()
        {
            const double rtStart = 31.0;
            const double rtEnd = 35.0;
            const double qValueThreshold = 0.01;
            const double rtHalfWidth = 0.25;
            bool centroidOutput = GetSimulateCentroidedOutput();
            bool deduplicate = GetDeduplicateProteoforms();
            bool useGlobalAbundanceRefit = GetGlobalAbundanceRefitEnabled();
            int globalRefitMaxModels = GetGlobalAbundanceRefitMaxModels();

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string resultPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv");
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;

            string realSliceImspPath = Path.Combine(outDir, stem + ".rt31-35.real.imsp");
            string simSliceImspPath = Path.Combine(outDir, stem + ".rt31-35.simulated.q001.imsp");
            string simFullImspPath = Path.Combine(outDir, stem + ".full.simulated.q001.imsp");

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

            int realSlicePeakCount = WriteImspFile(rtSliceMs1Scans, realSliceImspPath, intensityThreshold: 10000);
            Console.WriteLine($"Real slice IMSP written: {realSliceImspPath}");
            Console.WriteLine($"Real slice scans: {rtSliceMs1Scans.Length}, peaks: {realSlicePeakCount}");

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

            var index = PeakIndexingEngine.InitializeIndexingEngine(allMs1Scans)!;
            var extractor = new GroundTruthExtractor(index, allMs1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

            var sliceFit = FitProteoforms(qFilteredSliceRecords, extractor, rtHalfWidth, useGlobalAbundanceRefit, globalRefitMaxModels);
            Assert.That(sliceFit.Models, Is.Not.Empty, "No simulated models were fitted for 31-35 min slice.");

            var sliceScanTimes = rtSliceMs1Scans.Select(s => s.RetentionTime).ToArray();
            var sliceSimulation = new Simulator().Simulate(
                sliceFit.Models,
                sliceFit.MinCharge,
                sliceFit.MaxCharge,
                sliceFit.SigmaMz,
                sliceScanTimes,
                pointsPerSigma: 1,
                mzPaddingInSigmas: 4.0);

            var sliceScansForExport = centroidOutput
                ? CentroidizeSimulatedScans(sliceSimulation.Scans, relativeIntensityThreshold: 1e-4)
                : sliceSimulation.Scans;

            int simSlicePeakCount = WriteImspFile(
                sliceScansForExport,
                simSliceImspPath,
                intensityThreshold: ComputeSimulationImspThreshold(sliceScansForExport, 1e-4, 1.0));

            Console.WriteLine($"Simulated slice IMSP written: {simSliceImspPath}");
            Console.WriteLine($"Simulated slice models: {sliceFit.Models.Length}, scans: {sliceSimulation.Scans.Length}, peaks: {simSlicePeakCount}");

            var fullFit = FitProteoforms(qFilteredRecords, extractor, rtHalfWidth, useGlobalAbundanceRefit, globalRefitMaxModels);
            Assert.That(fullFit.Models, Is.Not.Empty, "No simulated models were fitted for full q<=0.01 run.");

            var fullScanTimes = allMs1Scans.Select(s => s.RetentionTime).ToArray();
            var fullSimulation = new Simulator().Simulate(
                fullFit.Models,
                fullFit.MinCharge,
                fullFit.MaxCharge,
                fullFit.SigmaMz,
                fullScanTimes,
                pointsPerSigma: 1,
                mzPaddingInSigmas: 4.0);

            var fullScansForExport = centroidOutput
                ? CentroidizeSimulatedScans(fullSimulation.Scans, relativeIntensityThreshold: 1e-4)
                : fullSimulation.Scans;

            int simFullPeakCount = WriteImspFile(
                fullScansForExport,
                simFullImspPath,
                intensityThreshold: ComputeSimulationImspThreshold(fullScansForExport, 1e-4, 1.0));

            Console.WriteLine($"Simulated full IMSP written: {simFullImspPath}");
            Console.WriteLine($"Simulated full models: {fullFit.Models.Length}, scans: {fullSimulation.Scans.Length}, peaks: {simFullPeakCount}");
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

        private static (ProteoformModel[] Models, int MinCharge, int MaxCharge, double SigmaMz) FitProteoforms(
            IReadOnlyList<MmResultRecord> records,
            GroundTruthExtractor extractor,
            double rtHalfWidth,
            bool useGlobalAbundanceRefit,
            int globalRefitMaxModels)
        {
            var fitter = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(fallbackSigmaMz: 0.012));
            var fitted = new List<FittedProteoform>(records.Count);
            var truths = new List<ProteoformGroundTruth>(records.Count);

            int minChargeGlobal = int.MaxValue;
            int maxChargeGlobal = int.MinValue;

            for (int i = 0; i < records.Count; i++)
            {
                var record = records[i];
                int minCharge = Math.Max(2, record.PrecursorCharge - 2);
                int maxCharge = Math.Min(80, record.PrecursorCharge + 2);
                if (minCharge > maxCharge)
                    continue;

                if (i % 25 == 0)
                    Console.WriteLine($"Fitting record {i + 1}/{records.Count}: {record.Identifier}");

                var truth = extractor.Extract(record.MonoisotopicMass, record.RetentionTime, rtHalfWidth, minCharge, maxCharge);
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
                truths.Add(truth);
                minChargeGlobal = Math.Min(minChargeGlobal, minCharge);
                maxChargeGlobal = Math.Max(maxChargeGlobal, maxCharge);
            }

            if (minChargeGlobal == int.MaxValue)
                minChargeGlobal = 2;
            if (maxChargeGlobal == int.MinValue)
                maxChargeGlobal = 80;

            var sigmaCandidates = fitted
                .Select(f => f.SigmaMz)
                .Where(s => !double.IsNaN(s) && !double.IsInfinity(s) && s > 0)
                .OrderBy(s => s)
                .ToArray();
            double sigmaMz = sigmaCandidates.Length == 0 ? 0.012 : sigmaCandidates[sigmaCandidates.Length / 2];

            if (useGlobalAbundanceRefit && fitted.Count > 1)
            {
                if (fitted.Count > globalRefitMaxModels)
                {
                    Console.WriteLine($"Skipping global abundance refit for {fitted.Count} models (max allowed: {globalRefitMaxModels}).");
                }
                else
                {
                var refitter = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(
                    MaxIterations: 8,
                    ConvergenceTolerance: 1e-3,
                    MinimumAbundance: 0,
                    Verbose: true));

                var refitResult = refitter.Refit(fitted, truths, minChargeGlobal, maxChargeGlobal, sigmaMz);
                fitted = refitResult.FittedProteoforms.ToList();
                Console.WriteLine($"Global abundance refit iterations: {refitResult.IterationsCompleted}, converged: {refitResult.Converged}");
                Console.WriteLine($"Global abundance refit residual fraction: {refitResult.InitialResidualFraction:G6} -> {refitResult.FinalResidualFraction:G6}");
                }
            }

            return (
                fitted.Select(f => f.Model).ToArray(),
                minChargeGlobal,
                maxChargeGlobal,
                sigmaMz);
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

        private static double ComputeSimulationImspThreshold(
            MsDataScan[] scans,
            double fractionOfMaxIntensity,
            double minimumThreshold)
        {
            double maxIntensity = 0;
            foreach (var scan in scans)
            {
                if (!scan.MassSpectrum.YArray.Any())
                    continue;

                double localMax = scan.MassSpectrum.YArray.Max();
                if (localMax > maxIntensity)
                    maxIntensity = localMax;
            }

            return Math.Max(minimumThreshold, maxIntensity * fractionOfMaxIntensity);
        }

        private static MsDataScan[] CentroidizeSimulatedScans(MsDataScan[] profileScans, double relativeIntensityThreshold)
        {
            var centroided = new MsDataScan[profileScans.Length];

            for (int s = 0; s < profileScans.Length; s++)
            {
                var scan = profileScans[s];
                var x = scan.MassSpectrum.XArray;
                var y = scan.MassSpectrum.YArray;

                if (x.Length == 0)
                {
                    centroided[s] = CloneScanWithSpectrum(scan, Array.Empty<double>(), Array.Empty<double>(), isCentroid: true);
                    continue;
                }

                double maxIntensity = y.Max();
                double floor = Math.Max(0.0, maxIntensity * relativeIntensityThreshold);

                var mzList = new List<double>();
                var intList = new List<double>();

                for (int i = 1; i < y.Length - 1; i++)
                {
                    double current = y[i];
                    if (current < floor)
                        continue;

                    bool isPeak = current >= y[i - 1] && current >= y[i + 1]
                                  && (current > y[i - 1] || current > y[i + 1]);
                    if (!isPeak)
                        continue;

                    mzList.Add(x[i]);
                    intList.Add(current);
                }

                centroided[s] = CloneScanWithSpectrum(scan, mzList.ToArray(), intList.ToArray(), isCentroid: true);
            }

            return centroided;
        }

        private static MsDataScan CloneScanWithSpectrum(MsDataScan source, double[] mz, double[] intensities, bool isCentroid)
        {
            double tic = intensities.Sum();
            MzRange scanWindowRange = mz.Length > 0
                ? new MzRange(mz[0], mz[^1])
                : source.ScanWindowRange;

            return new MsDataScan(
                massSpectrum: new MzSpectrum(mz, intensities, false),
                oneBasedScanNumber: source.OneBasedScanNumber,
                msnOrder: source.MsnOrder,
                isCentroid: isCentroid,
                polarity: source.Polarity,
                retentionTime: source.RetentionTime,
                scanWindowRange: scanWindowRange,
                scanFilter: source.ScanFilter,
                mzAnalyzer: source.MzAnalyzer,
                totalIonCurrent: tic,
                injectionTime: source.InjectionTime,
                noiseData: source.NoiseData,
                nativeId: source.NativeId,
                selectedIonMz: source.SelectedIonMZ,
                selectedIonChargeStateGuess: source.SelectedIonChargeStateGuess,
                selectedIonIntensity: source.SelectedIonIntensity,
                isolationMZ: source.IsolationMz,
                isolationWidth: source.IsolationWidth,
                dissociationType: source.DissociationType,
                oneBasedPrecursorScanNumber: source.OneBasedPrecursorScanNumber,
                selectedIonMonoisotopicGuessMz: source.SelectedIonMonoisotopicGuessMz,
                hcdEnergy: source.HcdEnergy,
                scanDescription: source.ScanDescription,
                compensationVoltage: source.CompensationVoltage);
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

        /// <summary>
        /// Generates fixture .imsp files for use in TypeScript parser tests.
        ///
        /// tiny-known  – 3 synthetic scans with exact known values; keep this one in the MsBrowser repo.
        /// small-real  – first 100 MS1 scans from the Jurkat raw (~300 KB–1 MB target).
        /// medium-real – first 500 MS1 scans                        (~5–15 MB target).
        /// large-real  – first 2000 MS1 scans                       (~50–100 MB target).
        ///
        /// Adjust the Take() counts if the resulting file sizes fall outside the desired ranges.
        /// </summary>
        [Test]
        public static void GenerateImspFixtures()
        {
            string fixtureDir = @"D:\JurkatTopdown\fixtures";
            Directory.CreateDirectory(fixtureDir);

            // ── tiny-known: synthetic data with exact, documentable expected values ──────
            // All peaks are kept (threshold = 0) so the file is fully deterministic.
            //
            // Expected header:  version=1, binsPerDalton=100, scanCount=3
            // Expected scans:
            //   [0] scanNumber=1, rt=1.0 min, tic=225000
            //   [1] scanNumber=2, rt=2.0 min, tic=250000
            //   [2] scanNumber=3, rt=3.0 min, tic=210000
            // Expected peaks (mzTenThousandths / intensity / scanIndex):
            //   bin 5001: (5001234, 50000, 0)  (5001234, 60000, 1)
            //   bin 6002: (6002345, 45000, 2)
            //   bin 7506: (7505678, 100000, 0)  (7505678, 110000, 1)
            //   bin 8507: (8506789, 95000, 2)
            //   bin 10009: (10009012, 75000, 0)  (10009012, 70000, 2)
            //   bin 12503: (12503456, 80000, 1)
            var tinyScans = new[]
            {
                MakeScan(1, 1.0, new[] { (500.1234, 50000.0),  (750.5678, 100000.0), (1000.9012, 75000.0)  }),
                MakeScan(2, 2.0, new[] { (500.1234, 60000.0),  (750.5678, 110000.0), (1250.3456, 80000.0)  }),
                MakeScan(3, 3.0, new[] { (600.2345, 45000.0),  (850.6789, 95000.0),  (1000.9012, 70000.0)  }),
            };
            WriteImspFile(tinyScans, Path.Combine(fixtureDir, "tiny-known.imsp"), intensityThreshold: 0);

            // ── real fixtures sliced from the Jurkat raw file ────────────────────────────
            string rawPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var allMs1 = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            WriteImspFile(allMs1.Take(100).ToArray(),  Path.Combine(fixtureDir, "small-real.imsp"));
            WriteImspFile(allMs1.Take(500).ToArray(),  Path.Combine(fixtureDir, "medium-real.imsp"));
            WriteImspFile(allMs1.Take(2000).ToArray(), Path.Combine(fixtureDir, "large-real.imsp"));
        }

        /// <summary>
        /// Core IMSP writer. Bins peaks from the supplied MS1 scans and serialises them
        /// to the given output path. Returns the total number of peaks written.
        /// </summary>
        public static int WriteImspFile(MsDataScan[] ms1Scans, string outputPath,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold)
        {
            return new ImspExportService().WriteFile(ms1Scans, outputPath, intensityThreshold: intensityThreshold);
        }

        /// <summary>Creates a minimal synthetic MS1 scan for fixture generation.</summary>
        private static MsDataScan MakeScan(int oneBasedScanNumber, double retentionTime,
            (double mz, double intensity)[] peaks)
        {
            double[] mzArr  = peaks.Select(p => p.mz).ToArray();
            double[] intArr = peaks.Select(p => p.intensity).ToArray();
            return new MsDataScan(
                massSpectrum:    new MzSpectrum(mzArr, intArr, false),
                oneBasedScanNumber: oneBasedScanNumber,
                msnOrder:        1,
                isCentroid:      true,
                polarity:        Polarity.Positive,
                retentionTime:   retentionTime,
                scanWindowRange: new MzRange(mzArr.Min(), mzArr.Max()),
                scanFilter:      "FTMS + p NSI Full ms",
                mzAnalyzer:      MZAnalyzerType.Orbitrap,
                totalIonCurrent: intArr.Sum(),
                injectionTime:   1.0,
                noiseData:       null,
                nativeId:        "scan=" + oneBasedScanNumber);
        }

        private const string ImspFormatMarkdown =
@"# IMSP Peak Index File Format (`.imsp`)

## Overview

An `.ind` file contains all MS1 peaks from a single mass-spec run, organised into
a sparse m/z bin index.  The binning strategy is identical to that used by the
FlashLFQ `PeakIndexingEngine`: every peak is placed in a bin whose index is

    BinIndex = round(mz × BinsPerDalton)

where `BinsPerDalton = 100`.  Peaks in bin 1000 therefore cover m/z values that
round to 10.00 Th (roughly 9.995–10.005 Th).

All numeric values are **little-endian**.

---

## File Layout

```
[ Section 1 ] File Header          –  24 bytes  (fixed)
[ Section 2 ] Scan Table           –  S × 16 bytes
[ Section 3 ] Bin Directory        –  N × 12 bytes
[ Section 4 ] Peak Array           –  T × 12 bytes
```

where:
- **S** = number of MS1 scans  (`ScanCount` in the header)
- **N** = number of non-empty bins  (`NonEmptyBinCount` in the header)
- **T** = total number of peaks  (`TotalPeakCount` in the header)

---

## Section 1 — File Header (offset 0, 24 bytes)

| Offset | Bytes | Type    | Field            | Description                          |
|-------:|------:|---------|------------------|--------------------------------------|
|      0 |     4 | ASCII   | Magic            | Always the four bytes `I M S P`      |
|      4 |     4 | uint32  | Version          | Format version; currently `1`        |
|      8 |     4 | uint32  | BinsPerDalton    | Bin resolution; currently `100`      |
|     12 |     4 | uint32  | NonEmptyBinCount | Number of non-empty bins (N)         |
|     16 |     4 | uint32  | TotalPeakCount   | Total peaks across all bins (T)      |
|     20 |     4 | uint32  | ScanCount        | Number of MS1 scans (S)              |

---

## Section 2 — Scan Table (offset 24, S × 16 bytes)

One entry per MS1 scan, ordered by ascending scan number.
The zero-based position of an entry is the `ScanIndex` referenced in the peak array.

| Relative offset | Bytes | Type    | Field               | Description                                              |
|----------------:|------:|---------|---------------------|----------------------------------------------------------|
|               0 |     4 | uint32   | OneBasedScanNumber  | Instrument scan number                                   |
|               4 |     8 | float64 | RetentionTime       | Retention time (minutes)                                 |
|              12 |     4 | float32 | Tic                 | Sum of filtered peak intensities; use for TIC chromatogram |

---

## Section 3 — Bin Directory (offset `24 + S×16`, N × 12 bytes)

One entry per non-empty bin, ordered by ascending `BinIndex`.
Use this section for fast m/z-range queries without scanning the full peak array.

| Relative offset | Bytes | Type    | Field       | Description                                          |
|----------------:|------:|---------|-------------|------------------------------------------------------|
|               0 |     4 | uint32  | BinIndex    | `round(representativeMz × BinsPerDalton)`            |
|               4 |     4 | uint32  | PeakOffset  | 0-based index of the first peak for this bin in §4   |
|               8 |     4 | uint32  | PeakCount   | Number of peaks belonging to this bin                |

---

## Section 4 — Peak Array (offset `24 + S×16 + N×12`, T × 12 bytes)

Peaks are stored contiguously, bin-by-bin in the same order as the bin directory.
Each record is exactly **12 bytes**, so the byte offset of peak `i` within this
section is `i × 12`.

Retention time is not stored per-peak; look it up via `scans[ZeroBasedScanIndex].retentionTime`.

| Relative offset | Bytes | Type    | Field              | Notes                                                           |
|----------------:|------:|---------|--------------------|------------------------------------------------------------------|
|               0 |     4 | uint32  | MzTenThousandths   | `round(mz × 10000)` — recover with `/ 10000`; 4 d.p. max       |
|               4 |     4 | float32 | Intensity          | Raw detector intensity (always positive; covers up to ~3.4e38)  |
|               8 |     4 | uint32  | ZeroBasedScanIndex | Index into the Scan Table (Section 2)                           |

---

## TypeScript Reader (DataView)

```typescript
const HEADER_BYTES = 24;
const SCAN_BYTES   = 16;
const BIN_BYTES    = 12;
const PEAK_BYTES   = 12;

interface ImspHeader {
  binsPerDalton:    number;
  nonEmptyBinCount: number;
  totalPeakCount:   number;
  scanCount:        number;
}

interface ImspScan {
  oneBasedScanNumber: number;
  retentionTime:      number;
  tic:                number;
}

interface ImspBinEntry {
  binIndex:   number;
  peakOffset: number;
  peakCount:  number;
}

interface ImspPeak {
  mz:        number;  // Th, 4 d.p.
  intensity: number;
  scanIndex: number;
}

function parseImsp(buffer: ArrayBuffer) {
  const v = new DataView(buffer);

  // ── Header ───────────────────────────────────────────────────────────────────
  const magic = String.fromCharCode(
    v.getUint8(0), v.getUint8(1), v.getUint8(2), v.getUint8(3)
  );
  if (magic !== 'IMSP') throw new Error('Not an IMSP file');

  const header: ImspHeader = {
    binsPerDalton:    v.getUint32( 8, true),
    nonEmptyBinCount: v.getUint32(12, true),
    totalPeakCount:   v.getUint32(16, true),
    scanCount:        v.getUint32(20, true),
  };

  // ── Derived offsets ───────────────────────────────────────────────────────────
  const scanTableStart = HEADER_BYTES;
  const binDirStart    = scanTableStart + header.scanCount        * SCAN_BYTES;
  const peakArrayStart = binDirStart    + header.nonEmptyBinCount * BIN_BYTES;

  // ── Scan Table ────────────────────────────────────────────────────────────────
  const scans: ImspScan[] = [];
  for (let i = 0; i < header.scanCount; i++) {
    const off = scanTableStart + i * SCAN_BYTES;
    scans.push({
      oneBasedScanNumber: v.getUint32 (off,      true),
      retentionTime:      v.getFloat64(off +  4, true),
      tic:                v.getFloat32(off + 12, true),
    });
  }

  // ── Bin Directory ─────────────────────────────────────────────────────────────
  const bins: ImspBinEntry[] = [];
  for (let i = 0; i < header.nonEmptyBinCount; i++) {
    const off = binDirStart + i * BIN_BYTES;
    bins.push({
      binIndex:   v.getUint32(off,     true),
      peakOffset: v.getUint32(off + 4, true),
      peakCount:  v.getUint32(off + 8, true),
    });
  }

  // ── Peak accessor ─────────────────────────────────────────────────────────────
  function readPeak(peakIndex: number): ImspPeak {
    const off = peakArrayStart + peakIndex * PEAK_BYTES;
    return {
      mz:        v.getUint32 (off,     true) / 10000,  // ten-thousandths Th → Th
      intensity: v.getFloat32(off + 4, true),
      scanIndex: v.getUint32 (off + 8, true),
    };
  }

  // ── Convenience: peaks for an m/z range ──────────────────────────────────────
  function peaksInMzRange(mzMin: number, mzMax: number): ImspPeak[] {
    const binMin = Math.floor(mzMin * header.binsPerDalton);
    const binMax = Math.ceil (mzMax * header.binsPerDalton);
    const result: ImspPeak[] = [];
    for (const bin of bins) {
      if (bin.binIndex < binMin) continue;
      if (bin.binIndex > binMax) break;
      for (let k = 0; k < bin.peakCount; k++)
        result.push(readPeak(bin.peakOffset + k));
    }
    return result;
  }

  return { header, scans, bins, readPeak, peaksInMzRange };
}
```

---

## Visualisation Notes

**TIC chromatogram** — render immediately from the scan table alone, no peak data required:

| Axis | Source field           |
|------|------------------------|
| X    | `scans[i].retentionTime` |
| Y    | `scans[i].tic`           |

**XIC (Extracted Ion Chromatogram)** — shows the intensity of a specific m/z across retention time.
Use `peaksInMzRange(mzMin, mzMax)` to retrieve peaks within a narrow m/z window, then plot:

| Axis | Source field                                    |
|------|-------------------------------------------------|
| X    | `scans[peak.scanIndex].retentionTime`           |
| Y    | `peak.intensity`                                |

The bin directory makes m/z-range extraction efficient — only the relevant bins are read,
leaving the rest of the peak array untouched.

---

## Version Notes

### v1 (2026-04-11)
- Initial format definition
- Magic bytes: `IMSP`
- Canonical file extension: `.imsp`
- Scan table entries: 16 bytes (uint32 scanNumber + float64 RT + float32 TIC)
- Bin directory entries: 12 bytes (uint32 BinIndex + uint32 PeakOffset + uint32 PeakCount)
- Peak records: 12 bytes (uint32 mzTenThousandths + float32 intensity + uint32 scanIndex)
- Intensity threshold applied during indexing: peaks below threshold are excluded from both the peak array and the TIC
- All values little-endian

### v1 Clarifications
- 2026-04-11: Confirmed scan table entry size is 16 bytes (earlier drafts incorrectly stated 12)
- 2026-04-11: Confirmed `OneBasedScanNumber` is written and read as uint32
";
    }
}
