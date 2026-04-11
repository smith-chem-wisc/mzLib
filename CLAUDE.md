# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository layout

The git repo root contains `README.md` and a single `mzLib/` subdirectory that holds the Visual Studio solution (`mzLib/mzLib.sln`) and all project folders. **All build/test commands below must be run from `mzLib/`**, not the repo root. CI does this with `cd mzLib && dotnet ...`.

## Build, test, and run

Target framework is `net8.0` (most libraries) and `net8.0-windows` for `Test` and `mzPlot`. Platform is `x64` everywhere. CI runs on `windows-latest` — the Test project references `mzPlot` (OxyPlot.Wpf) and Windows-only Thermo/Bruker native DLLs, so the full test suite is only reliably buildable on Windows. On macOS/Linux, non-Windows projects build fine individually but `Test/Test.csproj` will fail to restore.

```bash
cd mzLib
dotnet restore
dotnet build --no-restore                               # Debug
dotnet build --no-restore --configuration Release       # Release (used by the integration job + nuspec)
dotnet build --no-restore ./Test/Test.csproj            # build test project explicitly
dotnet test  --no-build   ./Test/Test.csproj            # run all tests
```

Run a single test or fixture by fully qualified name:

```bash
dotnet test ./Test/Test.csproj --filter "FullyQualifiedName~TestFlashLFQ.SomeTest"
dotnet test ./Test/Test.csproj --filter "FullyQualifiedName~IndexingEngineTests"
```

Coverage (matches CI):

```bash
dotnet add Test/Test.csproj package coverlet.collector -v 6.0.2
dotnet test --no-build --collect:"XPlat Code Coverage" /p:CoverletOutputFormat=cobertura ./Test/Test.csproj
```

mzLib is a library and has no runnable entry point. `Development/` is a scratch project for experimentation; `Test/` is the NUnit 4 test project (also excluded from code coverage via an assembly-level `[ExcludeFromCodeCoverage]` attribute declared in its csproj).

## Packaging and the MetaMorpheus contract

`mzLib` ships as a NuGet package (`mzLib.nuspec`) consumed primarily by [MetaMorpheus](https://github.com/smith-chem-wisc/MetaMorpheus). The CI `integration` job packs mzLib with a stamped `9.9.9` version, adds it as a local NuGet source, swaps MetaMorpheus's `mzLib` package reference to `9.9.9`, and builds + tests MetaMorpheus against the PR. **Any breaking API change will fail the integration job even if all mzLib unit tests pass.** When renaming or removing public APIs, grep MetaMorpheus (or expect integration CI to catch it) before merging.

The `nuspec` `<files>` section hand-lists every DLL (and Thermo/Bruker native DLLs under `Readers/bin/.../Release/net8.0/`) that ships in the package. New top-level projects or new native dependencies must be added there or they won't make it into the NuGet package.

## Architectural big picture

mzLib is a layered mass-spectrometry library. Dependencies flow strictly upward from small, foundational projects to higher-level engines:

1. **`MzLibUtil`** — shared primitives: `Tolerance` (PPM/absolute), `DoubleRange`, extension methods, exceptions. Depended on by everything.
2. **`Chemistry`** — elements, isotopes, `ChemicalFormula`, `ClassExtensions.ToMass/ToMz`. The periodic table is loaded statically via `UsefulProteomicsDatabases.Loaders.LoadElements()`.
3. **`MassSpectrometry`** — core spectra types: `MsDataFile`, `MsDataScan`, `MzSpectrum`, isotopic envelopes, deconvolution algorithms (`ClassicDeconvolutionAlgorithm`, `IsoDecAlgorithm`), and the generic **peak-indexing subsystem** under `MassSpectrometry/PeakIndexing/` (`IIndexedPeak`, `IndexedMassSpectralPeak`, `IndexedMass`, abstract `IndexingEngine<T>`, concrete `MassIndexingEngine`, plus `ExtractedIonChromatogram` with XIC smoothing/boundary detection).
4. **`Omics`**, **`Proteomics`**, **`Transcriptomics`** — biopolymer types. `Omics` defines shared interfaces (`IBioPolymer`, modifications, digestion, sequence conversion); `Proteomics` and `Transcriptomics` are parallel concrete implementations. Treat entrapment/contaminant/decoy flags uniformly — see recent commit "Treat Entrapment like Contaminants".
5. **`UsefulProteomicsDatabases`** — loaders for UniProt XML, FASTA, PTM lists, Unimod, element data. Also handles protein XML **writing** (see recent "Xml protein writer fix").
6. **`Readers`** — file-format adapters (`Mzml`, `ThermoRawFileReader`, Bruker TIMS, MGF, etc.) that produce `MsDataFile`. Thermo and Bruker readers pull in native DLLs copied to the output directory; these are the reason Test is `net8.0-windows`.
7. **`MzIdentML`**, **`pepXML`** — serializers for identification result formats.
8. **`Chromatography`**, **`SpectralAveraging`**, **`BayesianEstimation`** — analytical utilities layered on top of `MassSpectrometry`.
9. **`FlashLFQ`** + **`Quantification`** — the label-free quantification engine. `FlashLfqEngine` is the top-level orchestrator; it owns a per-`SpectraFileInfo` dictionary of `PeakIndexingEngine` instances (`FlashLFQ/PeakIndexingEngine/PeakIndexingEngine.cs`, which is `IndexingEngine<IndexedMassSpectralPeak>` from `MassSpectrometry`) and uses them to build XICs, trace isotopic envelopes, do match-between-runs (`FlashLFQ/MBR/`) and isotope tracking (`FlashLFQ/IsoTracker/`).
10. **`PredictionClients`** — gRPC/protobuf clients for external ML inference services (Koina). Recently refactored.

Two points that are easy to miss:

- **`InternalsVisibleTo`**: several library csprojs grant `internal` access to `Development` and `Test`. When adding tests that poke at internals, check the target project's csproj — not all of them expose internals.
- **Peak indexing is generic**: `IndexingEngine<T>` in `MassSpectrometry/PeakIndexing/` is the reusable base. `PeakIndexingEngine` (in `FlashLFQ/PeakIndexingEngine/`) is the concrete m/z-based engine used by FlashLFQ; `MassIndexingEngine` (in `MassSpectrometry/PeakIndexing/`) indexes deconvoluted masses instead. New indexing variants should subclass `IndexingEngine<T>` rather than duplicating the binary-search/bin logic.

## Test-data conventions

Test data files live under `Test/` in area-specific subfolders (`Test/AveragingTests/TestData/`, `Test/DatabaseTests/`, `Test/FlashLFQ/...`, etc.) and are shipped to the test bin directory via `<None Update="..." CopyToOutputDirectory="..." />` entries in `Test/Test.csproj`. **When adding a new test data file, you must also add a matching entry to `Test/Test.csproj`** or the file will not be found at runtime.
