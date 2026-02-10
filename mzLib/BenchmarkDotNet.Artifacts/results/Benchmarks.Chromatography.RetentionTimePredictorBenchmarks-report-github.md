```

BenchmarkDotNet v0.14.0, Windows 11 (10.0.26100.7623)
AMD Ryzen Threadripper PRO 5975WX 32-Cores, 1 CPU, 64 logical and 32 physical cores
.NET SDK 10.0.102
  [Host]   : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2
  ShortRun : .NET 8.0.23 (8.0.2325.60607), X64 RyuJIT AVX2

Job=ShortRun  IterationCount=5  LaunchCount=1  
WarmupCount=3  

```
| Method              | PeptideCount | Mean        | Error     | StdDev   | Ratio    | RatioSD | Gen0     | Allocated  | Alloc Ratio |
|-------------------- |------------- |------------:|----------:|---------:|---------:|--------:|---------:|-----------:|------------:|
| **SSRCalc3_PredictAll** | **10**           |    **86.87 μs** |  **1.615 μs** | **0.419 μs** | **baseline** |        **** |   **1.2207** |   **21.87 KB** |            **** |
|                     |              |             |           |          |          |         |          |            |             |
| **SSRCalc3_PredictAll** | **100**          |   **865.73 μs** |  **4.401 μs** | **1.143 μs** | **baseline** |        **** |  **12.6953** |  **218.67 KB** |            **** |
|                     |              |             |           |          |          |         |          |            |             |
| **SSRCalc3_PredictAll** | **1000**         | **8,659.87 μs** | **18.992 μs** | **4.932 μs** | **baseline** |        **** | **125.0000** | **2186.73 KB** |            **** |
