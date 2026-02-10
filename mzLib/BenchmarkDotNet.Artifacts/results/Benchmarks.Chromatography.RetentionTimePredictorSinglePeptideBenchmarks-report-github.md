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
| SSRCalc3_ShortPeptide | 16.30 μs | 0.077 μs | 0.020 μs |     baseline |         | 0.0916 |   1.95 KB |             |
| SSRCalc3_LongPeptide  | 24.60 μs | 0.178 μs | 0.046 μs | 1.51x slower |   0.00x | 1.4954 |  24.56 KB | 12.63x more |
