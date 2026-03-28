using BenchmarkDotNet.Running;

BenchmarkSwitcher
	.FromAssembly(typeof(Development.MSL.MslBenchmarks).Assembly)
	.Run(args);