using BenchmarkDotNet.Running;
using Benchmarks;

// Run all benchmarks in the assembly with command-line argument support
// Examples:
//   dotnet run -c Release                           # Run all benchmarks
//   dotnet run -c Release -- --filter "*SSRCalc*"   # Run only SSRCalc benchmarks
//   dotnet run -c Release -- --filter "*SinglePeptide*"  # Run single peptide benchmarks
//   dotnet run -c Release -- --list flat            # List all available benchmarks

var config = new MzLibBenchmarkConfig();

#if DEBUG
Console.ForegroundColor = ConsoleColor.Yellow;
Console.WriteLine("WARNING: Running benchmarks in DEBUG mode. Results will be unreliable.");
Console.WriteLine("Run with: dotnet run -c Release");
Console.ResetColor();
Console.WriteLine();
#endif

if (args.Length == 0)
{
    // Interactive mode - let user choose which benchmark to run
    BenchmarkSwitcher
        .FromAssembly(typeof(Program).Assembly)
        .Run(args, config);
}
else
{
    // Command-line mode with filters
    BenchmarkSwitcher
        .FromAssembly(typeof(Program).Assembly)
        .Run(args);
}