```

BenchmarkDotNet v0.14.0, Windows 11 (10.0.26100.7623)
AMD Ryzen Threadripper PRO 5975WX 32-Cores, 1 CPU, 64 logical and 32 physical cores
.NET SDK 10.0.102
  [Host]   : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2
  ShortRun : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2

Job=ShortRun  IterationCount=5  LaunchCount=1  
WarmupCount=3  

```
| Method                | Mean     | Error    | StdDev   | Ratio        | RatioSD | Gen0   | Allocated | Alloc Ratio |
|---------------------- |---------:|---------:|---------:|-------------:|--------:|-------:|----------:|------------:|
| SSRCalc3_ShortPeptide | 14.37 μs | 0.070 μs | 0.011 μs |     baseline |         | 0.0458 |     840 B |             |
| SSRCalc3_LongPeptide  | 21.31 μs | 0.478 μs | 0.124 μs | 1.48x slower |   0.01x | 0.9460 |   16080 B | 19.14x more |
