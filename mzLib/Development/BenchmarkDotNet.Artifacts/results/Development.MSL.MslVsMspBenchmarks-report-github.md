```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.8037/25H2/2025Update/HudsonValley2)
Intel Core i7-8550U CPU 1.80GHz (Max: 1.99GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.103
  [Host]   : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3
  .NET 8.0 : .NET 8.0.21 (8.0.21, 8.0.2125.47513), X64 RyuJIT x86-64-v3

Job=.NET 8.0  Runtime=.NET 8.0  LaunchCount=Default  
WarmupCount=Default  

```
| Method                     | NPrecursors | AvgFrag | Mean              | Error            | StdDev            | Gen0       | Gen1      | Gen2     | Allocated   |
|--------------------------- |------------ |-------- |------------------:|-----------------:|------------------:|-----------:|----------:|---------:|------------:|
| **MslWrite_Full**              | **1000**        | **10**      |   **6,571,168.53 ns** |   **105,264.039 ns** |     **93,313.829 ns** |    **46.8750** |         **-** |        **-** |    **211925 B** |
| MspWrite_Full              | 1000        | 10      |  10,417,464.73 ns |   194,835.656 ns |    172,716.734 ns |   937.5000 |   62.5000 |  62.5000 |   4038008 B |
| MslLoad_Full               | 1000        | 10      |   2,200,724.79 ns |    43,201.530 ns |     69,762.424 ns |   355.4688 |  242.1875 |  74.2188 |   1717937 B |
| MspLoad_Full               | 1000        | 10      |   5,816,493.25 ns |   113,087.944 ns |    138,882.163 ns |   765.6250 |  109.3750 |        - |   3237847 B |
| MslLoad_IndexOnly          | 1000        | 10      |   2,246,101.35 ns |    25,985.795 ns |     21,699.316 ns |   300.7813 |  187.5000 |  74.2188 |   1754202 B |
| MslLookup_BySequenceCharge | 1000        | 10      |          97.85 ns |         1.971 ns |          2.562 ns |     0.0095 |         - |        - |        40 B |
| MspLookup_BySequenceCharge | 1000        | 10      |                NA |               NA |                NA |         NA |        NA |       NA |          NA |
| **MslWrite_Full**              | **10000**       | **10**      |  **31,789,980.00 ns** |   **344,468.974 ns** |    **287,647.194 ns** |    **66.6667** |         **-** |        **-** |    **643940 B** |
| MspWrite_Full              | 10000       | 10      |  88,237,038.46 ns | 1,399,003.802 ns |  1,168,231.536 ns |  8666.6667 |         - |        - |  37419232 B |
| MslLoad_Full               | 10000       | 10      |  35,367,589.73 ns |   704,275.720 ns |  1,287,808.080 ns |  2375.0000 | 1375.0000 | 312.5000 |  17190196 B |
| MspLoad_Full               | 10000       | 10      |  57,850,788.89 ns | 1,145,902.198 ns |  1,226,102.515 ns |  7555.5556 |  111.1111 |        - |  32047400 B |
| MslLoad_IndexOnly          | 10000       | 10      |  82,690,752.14 ns | 1,607,092.608 ns |  1,850,730.238 ns |  2428.5714 | 1285.7143 | 285.7143 |  17514480 B |
| MslLookup_BySequenceCharge | 10000       | 10      |         104.95 ns |         1.974 ns |          2.112 ns |     0.0095 |         - |        - |        40 B |
| MspLookup_BySequenceCharge | 10000       | 10      |                NA |               NA |                NA |         NA |        NA |       NA |          NA |
| **MslWrite_Full**              | **50000**       | **10**      | **152,633,264.29 ns** | **2,732,450.693 ns** |  **2,422,246.368 ns** |   **250.0000** |         **-** |        **-** |   **2564016 B** |
| MspWrite_Full              | 50000       | 10      | 427,246,540.00 ns | 3,903,596.768 ns |  3,651,426.840 ns | 44000.0000 |         - |        - | 185780896 B |
| MslLoad_Full               | 50000       | 10      | 162,079,510.25 ns | 3,388,140.986 ns |  9,990,010.868 ns | 10250.0000 | 5250.0000 | 250.0000 |  84996708 B |
| MspLoad_Full               | 50000       | 10      | 280,469,875.00 ns | 4,748,718.936 ns |  4,209,615.648 ns | 38000.0000 |         - |        - | 160089336 B |
| MslLoad_IndexOnly          | 50000       | 10      | 426,574,461.76 ns | 8,305,729.876 ns | 13,412,206.634 ns | 10000.0000 | 5000.0000 |        - |  86600872 B |
| MslLookup_BySequenceCharge | 50000       | 10      |         106.83 ns |         2.655 ns |          7.788 ns |     0.0095 |         - |        - |        40 B |
| MspLookup_BySequenceCharge | 50000       | 10      |                NA |               NA |                NA |         NA |        NA |       NA |          NA |

Benchmarks with issues:
  MslVsMspBenchmarks.MspLookup_BySequenceCharge: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=1000, AvgFrag=10]
  MslVsMspBenchmarks.MspLookup_BySequenceCharge: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=10000, AvgFrag=10]
  MslVsMspBenchmarks.MspLookup_BySequenceCharge: .NET 8.0(Runtime=.NET 8.0) [NPrecursors=50000, AvgFrag=10]
