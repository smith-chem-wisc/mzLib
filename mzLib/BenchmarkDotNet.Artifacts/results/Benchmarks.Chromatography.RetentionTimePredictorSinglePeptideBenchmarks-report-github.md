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
| SSRCalc3_ShortPeptide | 18.95 μs | 0.122 μs | 0.032 μs |     baseline |         | 0.1526 |    2.5 KB |             |
| SSRCalc3_LongPeptide  | 28.08 μs | 0.233 μs | 0.060 μs | 1.48x slower |   0.00x | 1.6174 |  26.71 KB | 10.68x more |
