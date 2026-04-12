# Top-Down MS1 Simulator — Session Status

**Last updated:** 2026-04-11  
**Branch:** `TopDownSimulator`  
**Plan doc:** `agent_info/TopDown-Simulator-Plan.md`  
**Project root:** `mzLib/TopDownSimulator/`  
**Tests:** `mzLib/Test/TopDownSimulator/`

---

## Codebase summary

Every folder in the planned layout exists and has code. The **latest commit `8288a664` ("Changes partially implemented")** added the Extraction, Simulation, and Comparison layers all at once, but may not be fully building or passing tests yet. Verify with a build before continuing.

---

## Implementation status by component

### Model layer — complete

| File | What it does |
|---|---|
| `Model/ProteoformModel.cs` | Parameter record `θ_p = {MonoisotopicMass, Abundance, RtProfile, ChargeDistribution, Identifier}` |
| `Model/ForwardModel.cs` | Evaluates `I(t, mz) = Σ_p A_p · g_p(t) · Σ_z f_p(z) · φ(b; M_p, z, σ_m)`; has scalar `Evaluate(t, mz)` and grid `Rasterize(scanTimes, mzGrid)` |
| `Model/EmgProfile.cs` | Exponentially-modified Gaussian `g_p(t)` with params `(Mu, Sigma, Tau)` |
| `Model/ChargeStateDistribution.cs` | `IChargeStateDistribution` interface + `GaussianChargeDistribution(MuZ, SigmaZ)` |
| `Model/IsotopeEnvelopeKernel.cs` | Averagine-based isotopologue envelope; `CentroidMzs(z)` and `Evaluate(mz, z, sigmaMz)` |

### Fitting layer — complete (Phase 1 done, 19 tests pass)

The fitters are chained in sequence by `ParameterFitter`: `σ_m → g_p → f_p(z) → A_p`.

| File | Algorithm |
|---|---|
| `Fitting/EnvelopeWidthFitter.cs` | Pooled weighted 2nd-central-moment across isotopologues when window has ≥2 peaks; falls back to user-supplied σ_m for centroided data |
| `Fitting/RtProfileFitter.cs` | Closed-form EMG method of moments: `τ=(m3/2)^(1/3)`, `σ²=m2−τ²`, `μ=m1−τ`; Gaussian fallback when `m3 ≤ 0` |
| `Fitting/ChargeDistributionFitter.cs` | Weighted-moment Gaussian on z, derived from per-charge XIC apex intensities |
| `Fitting/AbundanceFitter.cs` | Closed-form linear least squares through the origin on `A·g(t)·f(z)·w_i/(σ_m·√(2π))` |
| `Fitting/ParameterFitter.cs` | Top-level orchestrator; chains all four fitters in sequence and returns a `ProteoformModel` |

### Extraction layer — code present, tests present, verify build

| File | What it does |
|---|---|
| `Extraction/PeakSample.cs` | `readonly record struct PeakSample(double Mz, double Intensity)` — leaf value in the peak-window tensor |
| `Extraction/ProteoformGroundTruth.cs` | Raw MS1 tensor around one proteoform: `IsotopologueIntensities[charge][iso][scan]`, `IsotopologuePeakWindows[charge][iso][scan][]`, `ChargeXics[charge][scan]` |
| `Extraction/GroundTruthExtractor.cs` | Wraps a `PeakIndexingEngine`, pulls the `ProteoformGroundTruth` tensor for a given `(mass, rtCenter, rtHalfWidth, minCharge, maxCharge)` |
| `Extraction/MmResultLoader.cs` | Parses a MetaMorpheus `.psmtsv` via `PsmFromTsvFile`; returns `IReadOnlyList<MmResultRecord>` |
| `Extraction/CoeluterFinder.cs` | Finds `MmResultRecord`s within `(rtHalfWidth, massHalfWidthDa)` of an anchor; Phase 3 helper |

### Simulation layer — code present, tests present, verify build

| File | What it does |
|---|---|
| `Simulation/GridRasterizer.cs` | Builds an m/z grid from the proteoform set (padding in σ_m units), rasterizes `ForwardModel` onto `(scanTimes × mzGrid)` → `RasterizedScanGrid` |
| `Simulation/ScanBuilder.cs` | Converts `RasterizedScanGrid` → `MsDataScan[]` → `GenericMsDataFile` |
| `Simulation/Simulator.cs` | High-level entry point: `Simulate(...)` → `SimulationResult`; `WriteMzml(...)` calls `MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra` |

### Comparison layer — code present, tests present, verify build

| File | What it does |
|---|---|
| `Comparison/ModelProjection.cs` | `internal` helper: projects `ForwardModel` onto the `ProteoformGroundTruth` centroid grid to get `predicted[charge][iso][scan]` |
| `Comparison/SpectralAngle.cs` | Cosine similarity + normalized spectral angle per scan; `ComputePerScan(truth, proteoforms, sigmaMz)` |
| `Comparison/XicCorrelation.cs` | Pearson correlation per charge between observed and predicted XICs |
| `Comparison/ResidualAnalyzer.cs` | Residual energy fraction `||y − ŷ||² / ||y||²` over the full tensor |
| `Comparison/ComparisonReport.cs` | `ComparisonReportBuilder.Create(truth, proteoforms, sigmaMz)` → `ComparisonReport` record with per-scan angles, per-charge correlations, residual summary, and convenience mean properties |

---

## Test files

All test classes live in `mzLib/Test/TopDownSimulator/`.

| Test file | Status | Notes |
|---|---|---|
| `AbundanceFitterTests.cs` | Passing (Phase 1) | |
| `ChargeDistributionFitterTests.cs` | Passing (Phase 1) | |
| `EmgProfileTests.cs` | Passing (Phase 1) | |
| `EnvelopeWidthFitterTests.cs` | Passing (Phase 1) | |
| `ForwardModelTests.cs` | Passing (Phase 1) | |
| `GroundTruthExtractorTests.cs` | Passing (Phase 1) | |
| `IsotopeEnvelopeKernelTests.cs` | Passing (Phase 1) | |
| `ParameterFitterTests.cs` | Passing (Phase 1) | |
| `RtProfileFitterTests.cs` | Passing (Phase 1) | |
| `MmResultLoaderTests.cs` | Added in latest commit — verify | Writes a temp `.psmtsv` and round-trips it through `MmResultLoader` |
| `CoeluterFinderTests.cs` | Added in latest commit — verify | Tests `FindCoeluters` with same-file/other-file/out-of-window filtering |
| `SimulationAndComparisonTests.cs` | Added in latest commit — verify | Two tests: `SimulatorBuildsConsistentMsDataFile` and `ComparisonMetricsAreNearPerfectForMatchingModel` (self-consistency check: fit a synthetic model, extract from it, compare — all metrics expected ≥ 0.999) |

Run the test suite:
```
cd mzLib && dotnet test ./Test/Test.csproj --filter "FullyQualifiedName~TopDownSimulator"
```

---

## What remains to be done

The code structure from the plan is fully scaffolded. After verifying the build passes:

1. **Fix any compilation/test failures** introduced in commit `8288a664`.
2. **Phase 2** — resolution-scaled `σ_m(m/z) ∝ m/z^1.5 / R` model in `EnvelopeWidthFitter`; currently uses a constant σ_m.
3. **Phase 3** — wire `CoeluterFinder` into a joint-NNLS abundance fit across co-eluting proteoforms.
4. **Phase 4** — noise model; population priors across all MM results in a file.
5. **Phase 5** — already structurally complete in `Simulator.cs`/`WriteMzml`; needs integration with a real data file and the `TopDownEngine` test harness (task T20.1).

---

## Key design decisions (locked)

1. **Parametric forward model**, not empirical re-synthesis.
2. **Sequential per-factor fitting** (`σ_m → g_p → f_p(z) → A_p`) — more debuggable than joint nonlinear least squares.
3. **Phase 1 uses a constant σ_m**; Phase 2 upgrades to `σ_m(m/z)` for Orbitrap.
4. **Log-intensity for shape fits, linear for abundance** — log handles dynamic range for shapes, NNLS stays clean.
5. **Start Gaussian for `f_p(z)`**; revisit skew-normal only if residuals show systematic asymmetry.
6. **All fitting helpers are `internal`** unless a cross-project consumer genuinely needs them.
7. **Project ships in the mzLib NuGet** — treat every public API as a potential MetaMorpheus breaking change.

---

## mzLib integration notes

- `InternalsVisibleTo` for `Test` project must be present in `TopDownSimulator.csproj` if tests need to access `internal` members.
- Do **not** add `TopDownSimulator` to `mzLib.nuspec` until the public API has settled.
- `PeakIndexingEngine.InitializeIndexingEngine(MsDataScan[])` (in `FlashLFQ`) is the factory used by `GroundTruthExtractorTests` — this is a FlashLFQ internal, so the Test project must reference FlashLFQ.
