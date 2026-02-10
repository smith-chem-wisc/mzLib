```

BenchmarkDotNet v0.14.0, Windows 11 (10.0.26100.7623)
AMD Ryzen Threadripper PRO 5975WX 32-Cores, 1 CPU, 64 logical and 32 physical cores
.NET SDK 10.0.102
  [Host]   : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2
  ShortRun : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2

Job=ShortRun  IterationCount=5  LaunchCount=1  
WarmupCount=3  

```
| Method              | PeptideCount | Mean        | Error     | StdDev   | Ratio    | RatioSD | Gen0    | Allocated  | Alloc Ratio |
|-------------------- |------------- |------------:|----------:|---------:|---------:|--------:|--------:|-----------:|------------:|
| **SSRCalc3_PredictAll** | **10**           |    **75.47 μs** |  **1.032 μs** | **0.160 μs** | **baseline** |        **** |  **0.7324** |   **13.18 KB** |            **** |
|                     |              |             |           |          |          |         |         |            |             |
| **SSRCalc3_PredictAll** | **100**          |   **758.27 μs** | **19.270 μs** | **5.004 μs** | **baseline** |        **** |  **7.8125** |   **131.8 KB** |            **** |
|                     |              |             |           |          |          |         |         |            |             |
| **SSRCalc3_PredictAll** | **1000**         | **7,505.09 μs** | **29.840 μs** | **7.749 μs** | **baseline** |        **** | **78.1250** | **1317.97 KB** |            **** |
