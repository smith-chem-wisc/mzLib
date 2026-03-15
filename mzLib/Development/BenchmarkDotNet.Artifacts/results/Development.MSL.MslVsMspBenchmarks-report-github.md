```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.8037/25H2/2025Update/HudsonValley2)
Intel Core i7-8550U CPU 1.80GHz (Max: 1.99GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.103
  [Host]   : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3
  .NET 8.0 : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3

Job=.NET 8.0  Runtime=.NET 8.0  

```
| Method                                | NPrecursors | Mean          | Error        | StdDev       | Ratio | RatioSD | Gen0       | Gen1      | Gen2     | Allocated    | Alloc Ratio |
|-------------------------------------- |------------ |--------------:|-------------:|-------------:|------:|--------:|-----------:|----------:|---------:|-------------:|------------:|
| **MspLoad_FullIndex**                     | **1000**        |   **5,624.30 μs** |   **102.366 μs** |   **143.503 μs** |  **1.00** |    **0.04** |   **781.2500** |         **-** |        **-** |    **3242.9 KB** |        **1.00** |
| MslLoad_Full                          | 1000        |   1,714.36 μs |    30.687 μs |    46.862 μs |  0.31 |    0.01 |   304.6875 |  222.6563 |  50.7813 |   1646.76 KB |        0.51 |
| MslLoad_IndexOnly                     | 1000        |            NA |           NA |           NA |     ? |       ? |         NA |        NA |       NA |           NA |           ? |
| Msp_TryGetSpectrum_1000Lookups        | 1000        |      63.35 μs |     1.251 μs |     1.948 μs |  0.01 |    0.00 |    16.2354 |         - |        - |     66.41 KB |        0.02 |
| Msl_TryGetLibrarySpectrum_1000Lookups | 1000        |     309.62 μs |     5.750 μs |     5.379 μs |  0.06 |    0.00 |   230.4688 |         - |        - |    941.41 KB |        0.29 |
|                                       |             |               |              |              |       |         |            |           |          |              |             |
| **MspLoad_FullIndex**                     | **10000**       |  **55,747.54 μs** |   **823.150 μs** |   **769.975 μs** | **1.000** |    **0.02** |  **7900.0000** |         **-** |        **-** |  **32458.23 KB** |       **1.000** |
| MslLoad_Full                          | 10000       |  32,098.65 μs |   638.699 μs |   709.912 μs | 0.576 |    0.01 |  2187.5000 | 1187.5000 | 187.5000 |  16236.42 KB |       0.500 |
| MslLoad_IndexOnly                     | 10000       |            NA |           NA |           NA |     ? |       ? |         NA |        NA |       NA |           NA |           ? |
| Msp_TryGetSpectrum_1000Lookups        | 10000       |      63.66 μs |     1.244 μs |     1.382 μs | 0.001 |    0.00 |    16.2354 |         - |        - |     66.41 KB |       0.002 |
| Msl_TryGetLibrarySpectrum_1000Lookups | 10000       |     317.12 μs |     6.137 μs |     6.821 μs | 0.006 |    0.00 |   230.4688 |         - |        - |    941.41 KB |       0.029 |
|                                       |             |               |              |              |       |         |            |           |          |              |             |
| **MspLoad_FullIndex**                     | **50000**       | **281,574.66 μs** | **4,601.400 μs** | **4,079.021 μs** | **1.000** |    **0.02** | **39500.0000** |         **-** |        **-** | **162719.36 KB** |       **1.000** |
| MslLoad_Full                          | 50000       | 174,007.68 μs | 3,465.631 μs | 9,007.637 μs | 0.618 |    0.03 | 10000.0000 | 5000.0000 | 250.0000 |  82663.68 KB |       0.508 |
| MslLoad_IndexOnly                     | 50000       |            NA |           NA |           NA |     ? |       ? |         NA |        NA |       NA |           NA |           ? |
| Msp_TryGetSpectrum_1000Lookups        | 50000       |      60.40 μs |     0.835 μs |     0.781 μs | 0.000 |    0.00 |    16.2354 |         - |        - |     66.41 KB |       0.000 |
| Msl_TryGetLibrarySpectrum_1000Lookups | 50000       |     312.33 μs |     5.940 μs |     5.266 μs | 0.001 |    0.00 |   230.4688 |         - |        - |    941.41 KB |       0.006 |

Benchmarks with issues:
  MslVsMspBenchmarks.MslLoad_IndexOnly: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=1000]
  MslVsMspBenchmarks.MslLoad_IndexOnly: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=10000]
  MslVsMspBenchmarks.MslLoad_IndexOnly: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=50000]
