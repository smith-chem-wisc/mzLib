```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.8037/25H2/2025Update/HudsonValley2)
Intel Core i7-8550U CPU 1.80GHz (Max: 1.99GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.103
  [Host]   : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3
  .NET 8.0 : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3

Job=.NET 8.0  Runtime=.NET 8.0  

```
| Method                        | NPrecursors | AvgFragmentsPerPrecursor | Mean                | Error             | StdDev            | Gen0       | Gen1       | Gen2     | Allocated   |
|------------------------------ |------------ |------------------------- |--------------------:|------------------:|------------------:|-----------:|-----------:|---------:|------------:|
| **Write_MslLibrary**              | **1000**        | **10**                       |     **6,308,421.63 ns** |    **121,365.720 ns** |    **134,897.676 ns** |    **31.2500** |          **-** |        **-** |    **186858 B** |
| FullLoad_MslLibrary           | 1000        | 10                       |     1,471,923.70 ns |     29,030.150 ns |     22,664.831 ns |   277.3438 |   179.6875 |  82.0313 |   1686281 B |
| IndexOnlyLoad_MslLibrary      | 1000        | 10                       |                  NA |                NA |                NA |         NA |         NA |       NA |          NA |
| QueryMzWindow_SingleQuery     | 1000        | 10                       |           910.49 ns |         11.648 ns |         10.896 ns |          - |          - |        - |           - |
| QueryMzWindow_1000Queries     | 1000        | 10                       |        14,707.94 ns |        279.328 ns |        261.284 ns |          - |          - |        - |           - |
| DdaLookup_ExistingEntry       | 1000        | 10                       |            87.00 ns |          1.743 ns |          1.630 ns |     0.0095 |          - |        - |        40 B |
| DdaLookup_MissingEntry        | 1000        | 10                       |            78.79 ns |          1.099 ns |          1.028 ns |     0.0134 |          - |        - |        56 B |
| BuildIndex_FromEntries        | 1000        | 10                       |       126,461.73 ns |      2,515.784 ns |      3,089.609 ns |    48.5840 |    12.2070 |        - |    204264 B |
| RtCalibration_LinearTransform | 1000        | 10                       |       834,567.16 ns |     14,232.355 ns |     13,978.076 ns |   191.4063 |   190.4297 |        - |   1204360 B |
| **Write_MslLibrary**              | **50000**       | **10**                       |   **159,134,713.33 ns** |  **3,178,565.350 ns** |  **3,660,440.585 ns** |   **333.3333** |          **-** |        **-** |   **2538987 B** |
| FullLoad_MslLibrary           | 50000       | 10                       |   206,522,559.02 ns |  6,179,747.638 ns | 17,928,569.739 ns | 10250.0000 |  5250.0000 | 500.0000 |  84647692 B |
| IndexOnlyLoad_MslLibrary      | 50000       | 10                       |                  NA |                NA |                NA |         NA |         NA |       NA |          NA |
| QueryMzWindow_SingleQuery     | 50000       | 10                       |         2,106.49 ns |         40.579 ns |         87.351 ns |          - |          - |        - |           - |
| QueryMzWindow_1000Queries     | 50000       | 10                       |        68,421.93 ns |      1,895.456 ns |      5,588.794 ns |          - |          - |        - |           - |
| DdaLookup_ExistingEntry       | 50000       | 10                       |            91.97 ns |          0.679 ns |          0.567 ns |     0.0095 |          - |        - |        40 B |
| DdaLookup_MissingEntry        | 50000       | 10                       |            78.12 ns |          1.436 ns |          2.626 ns |     0.0134 |          - |        - |        56 B |
| BuildIndex_FromEntries        | 50000       | 10                       |    16,070,786.07 ns |    606,369.738 ns |  1,700,325.879 ns |   890.6250 |   390.6250 |  78.1250 |  10645260 B |
| RtCalibration_LinearTransform | 50000       | 10                       |   144,421,042.86 ns |  2,876,205.312 ns |  2,549,681.094 ns |  9000.0000 |  4750.0000 | 500.0000 |  60647292 B |
| **Write_MslLibrary**              | **500000**      | **10**                       | **1,437,155,880.00 ns** | **27,384,012.486 ns** | **25,615,022.282 ns** |  **3000.0000** |  **1000.0000** |        **-** |  **24139056 B** |
| FullLoad_MslLibrary           | 500000      | 10                       | 1,271,680,393.33 ns | 13,929,709.000 ns | 13,029,858.447 ns | 96000.0000 | 47000.0000 |        - | 840140320 B |
| IndexOnlyLoad_MslLibrary      | 500000      | 10                       |                  NA |                NA |                NA |         NA |         NA |       NA |          NA |
| QueryMzWindow_SingleQuery     | 500000      | 10                       |         1,797.82 ns |         27.973 ns |         23.359 ns |          - |          - |        - |           - |
| QueryMzWindow_1000Queries     | 500000      | 10                       |       799,423.09 ns |      7,850.750 ns |      6,959.486 ns |          - |          - |        - |           - |
| DdaLookup_ExistingEntry       | 500000      | 10                       |            90.31 ns |          1.357 ns |          1.133 ns |     0.0095 |          - |        - |        40 B |
| DdaLookup_MissingEntry        | 500000      | 10                       |            77.52 ns |          1.052 ns |          0.984 ns |     0.0134 |          - |        - |        56 B |
| BuildIndex_FromEntries        | 500000      | 10                       |   112,193,722.86 ns |  1,640,941.980 ns |  1,454,652.324 ns |  5000.0000 |   400.0000 |        - | 100138208 B |
| RtCalibration_LinearTransform | 500000      | 10                       | 1,038,034,938.46 ns | 13,940,074.934 ns | 11,640,593.922 ns | 83000.0000 | 41000.0000 |        - | 600138304 B |

Benchmarks with issues:
  MslBenchmarks.IndexOnlyLoad_MslLibrary: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=1000, AvgFragmentsPerPrecursor=10]
  MslBenchmarks.IndexOnlyLoad_MslLibrary: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=50000, AvgFragmentsPerPrecursor=10]
  MslBenchmarks.IndexOnlyLoad_MslLibrary: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=500000, AvgFragmentsPerPrecursor=10]
