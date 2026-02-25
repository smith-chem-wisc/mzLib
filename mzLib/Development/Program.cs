// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using Development.Dia;
using System;

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

            // Phase 1–2: Index construction and access benchmarks
            DiaScanIndexBenchmark.RunAll();

            // Phase 5–6: Orchestration + CPU vs GPU benchmarks
            DiaOrchestrationBenchmark.RunAll();

            LibraryBridgeBenchmark.RunAll();

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
