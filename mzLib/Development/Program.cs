// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using Development.Dia;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Benchmarks;
using System;
using static MassSpectrometry.Dia.DiaLibraryQueryGenerator;

namespace Development
{
    /// <summary>
    /// Entry point for the Development project.
    /// 
    /// Runs all DIA engine benchmarks: index construction, parallel orchestration,
    /// and CPU vs GPU comparison (if GPU is available).
    /// 
    /// To run:
    ///   1. Right-click Development project → "Set as Startup Project"
    ///   2. Ctrl+F5 (Start Without Debugging — important for accurate timing)
    /// 
    /// Or from command line:
    ///   dotnet run --project Development --configuration Release
    /// </summary>
    public class Program
    {
        public static async Task Main(string[] args)
        {
            Console.WriteLine("mzLib Development Benchmarks");
            Console.WriteLine(new string('=', 60));
            Console.WriteLine();

            //DiaScanIndexBenchmark.RunAll();

            //DiaOrchestrationBenchmark.RunAll();

            //LibraryBridgeBenchmark.RunAll();

            //RtCalibrationBenchmark.RunAll();

            //RealDataBenchmark.Run(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv"
            //    );

            //var scoredResults = Phase10ClassifierBenchmark.RunAll(
            //    @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    mspLibraryPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv");

            //Phase10_5_FeatureRefinementBenchmark.RunAll(scoredResults,
            //    @"F:\DiaBenchmark\PXD005573\phase10_5_features_koina.tsv");

            //Phase11FdrBenchmark.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase11_fdr_results.tsv");

            ////DecoyKoinaTableGenerator.Generate(
            ////    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv",
            ////    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.tsv"
            ////    );

            ////KoinaTableGenerator.Generate(
            ////    @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            ////    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv",
            ////    collisionEnergy: 27);

            //Phase12PeakGroupBenchmark.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase12_peak_group_results.tsv");

            ////DiaClassifierBenchmark.RunAll();

            //Phase13BenchmarkRunner.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //   outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase13_results.tsv");

            //Phase14BenchmarkRunner.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase14_results.tsv");

            //Phase15BenchmarkRunner.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase15_results.tsv");


            //Phase23BenchmarkRunner.RunAll(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase23_results.tsv",
            //    runABComparison: false);

            //MetaMorpheusEquivalentBenchmark.Run(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
            //    groundTruthPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputFolder: @"F:\DiaBenchmark\PXD005573\metamorpheusBenchmarkInMzLib.tsv");

            //PeakSelectionAccuracyBenchmark.Run(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
            //    groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    outputTsvPath: @"F:\DiaBenchmark\PXD005573\peakSelectionAccuracyBenchmark.tsv");

            //DiaBootstrapRunner.Run(
            //    rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
            //    mspPaths: new string[] { @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp", @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp" },
            //    ppmTol: 20,
            //    filenameDecoySuffix: true);

            await DiaSearchRunner.Run(
                rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
                //targetMspPath: @"F:\DiaBenchmark\PXD005573\koina_input.msp",
                targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diannFig2helalib.msp",
                //decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
                outputTsvPath: @"F:\DiaBenchmark\PXD005573\Runner\dia_search_results.tsv",
                groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_report.tsv",
                classifierType: DiaClassifierType.GradientBoostedTree);

            Console.WriteLine(new string('=', 60));
            Console.WriteLine("All benchmarks complete.");
            Console.WriteLine();

            if (!Console.IsInputRedirected)
            {
                Console.WriteLine("Press any key to exit.");
                Console.ReadKey();
            }
        }
    }
}
