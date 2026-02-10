```

BenchmarkDotNet v0.14.0, Windows 11 (10.0.26100.7623)
AMD Ryzen Threadripper PRO 5975WX 32-Cores, 1 CPU, 64 logical and 32 physical cores
.NET SDK 10.0.102
  [Host]   : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2
  ShortRun : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2

Job=ShortRun  IterationCount=5  LaunchCount=1  
WarmupCount=3  

```
| Method              | PeptideCount | Mean        | Error    | StdDev   | Ratio    | RatioSD | Gen0     | Allocated | Alloc Ratio |
|-------------------- |------------- |------------:|---------:|---------:|---------:|--------:|---------:|----------:|------------:|
| **SSRCalc3_PredictAll** | **10**           |    **100.6 μs** |  **0.76 μs** |  **0.20 μs** | **baseline** |        **** |   **1.4648** |  **25.67 KB** |            **** |
|                     |              |             |          |          |          |         |          |           |             |
| **SSRCalc3_PredictAll** | **100**          |  **1,008.7 μs** | **12.14 μs** |  **1.88 μs** | **baseline** |        **** |  **15.6250** | **256.72 KB** |            **** |
|                     |              |             |          |          |          |         |          |           |             |
| **SSRCalc3_PredictAll** | **1000**         | **10,110.3 μs** | **99.03 μs** | **15.33 μs** | **baseline** |        **** | **156.2500** | **2567.2 KB** |            **** |
