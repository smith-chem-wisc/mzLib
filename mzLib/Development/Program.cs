// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using Development.Dia;
using System;

namespace Development
{
    /// <summary>
    /// Entry point for the Development project.
    /// 
    /// This project is used for running benchmarks and development experiments
    /// that are too heavy or too informal for the unit test suite.
    /// 
    /// To run: right-click the Development project in Solution Explorer → 
    /// "Set as Startup Project", then press Ctrl+F5 (Start Without Debugging).
    /// Running without the debugger gives more accurate timing results.
    /// </summary>
    public class Program
    {
        public static void Main(string[] args)
        {
            Console.WriteLine("mzLib Development Benchmarks");
            Console.WriteLine(new string('=', 40));
            Console.WriteLine();

            DiaScanIndexBenchmark.RunAll();

            Console.WriteLine();
            Console.WriteLine("Done. Press any key to exit.");
            Console.ReadKey();
        }
    }
}
