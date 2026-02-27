// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using Development.Dia;
using MassSpectrometry.Dia;
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
        public static void Main(string[] args)
        {
            Console.WriteLine("mzLib Development Benchmarks");
            Console.WriteLine(new string('=', 60));
            Console.WriteLine();

            DiaScanIndexBenchmark.RunAll();

            DiaOrchestrationBenchmark.RunAll();

            LibraryBridgeBenchmark.RunAll();

            RtCalibrationBenchmark.RunAll();

            RealDataBenchmark.Run(
                mzmlPath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
                groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv"
                );

            var scoredResults = Phase10ClassifierBenchmark.RunAll(
                @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
                mspLibraryPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
                groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv");

            Phase10_5_FeatureRefinementBenchmark.RunAll(scoredResults,
                @"F:\DiaBenchmark\PXD005573\phase10_5_features_koina.tsv");

            Phase11FdrBenchmark.RunAll(
                rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
                targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
                decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
                groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
                outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase11_fdr_results.tsv");

            //DecoyKoinaTableGenerator.Generate(
            //    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv",
            //    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.tsv"
            //    );

            //KoinaTableGenerator.Generate(
            //    @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
            //    @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv",
            //    collisionEnergy: 27);

            Phase12PeakGroupBenchmark.RunAll(
                rawFilePath: @"F:\DiaBenchmark\PXD005573\Fig2HeLa-0-5h_MHRM_R01_T0.raw",
                targetMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.msp",
                decoyMspPath: @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoys.msp",
                groundTruthTsvPath: @"F:\DiaBenchmark\PXD005573\DiannOut\diann_ground_truth.tsv",
                outputTsvPath: @"F:\DiaBenchmark\PXD005573\phase12_peak_group_results.tsv");

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
