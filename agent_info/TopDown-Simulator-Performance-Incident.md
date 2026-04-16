# TopDown Simulator Performance Incident (Rep2 Fract7)

## What looked like a bug

`AnalysisExample.SimulateJurkatRep2AndWriteImsp` appeared to hang even when using a very small number of MM records (`Take = 1..3`).

Observed symptom: long periods with no output and runs exceeding 5 minutes.

## Root cause

This was not an infinite loop or failed numeric convergence.

Runtime was dominated by a combination of expensive steps:

1. Thermo RAW load (`LoadAllStaticData`) was very expensive on this dataset.
2. Peak indexing over all MS1 scans was moderately expensive.
3. Simulation generated dense profile grids.
4. IMSP writing used `intensityThreshold: 0`, causing excessive peak retention and write volume.

Additionally, forward-model evaluation originally performed broad summation work per grid point.

## Evidence (timed run with `Take = 3`)

From a detailed run of `AnalysisExample.SimulateJurkatRep2AndWriteImsp`:

- RAW load: `~1m58s`
- MS1 scan collection: `~0.003s`
- Peak index build: `~23.6s`
- MM load/filter (`3` records): `~1.8s`
- Fit loop (`3` records): seconds
- Simulation: `~11.0s`
- IMSP write: `~6.1s`
- Total: `~2m43s`

Conclusion: no non-terminating loop; the workload was simply heavy and front-loaded by file I/O.

## Fixes implemented

### 1) Kernel evaluation optimization

File: `mzLib/TopDownSimulator/Model/IsotopeEnvelopeKernel.cs`

- Added centroid caching by charge.
- Added sorted centroid/intensity evaluation cache by charge.
- Restricted Gaussian summation to isotopologues in a bounded m/z window (`+/- 6 sigma`).
- Added binary-search helpers (`LowerBound`, `UpperBound`) for fast window lookup.
- Added per-charge m/z bounds cache (`GetMzBounds`) used by rasterization.

### 2) Rasterization optimization

File: `mzLib/TopDownSimulator/Model/ForwardModel.cs`

- Reworked `Rasterize` to accumulate contributions species-by-species and charge-by-charge.
- For each charge, only evaluates m/z bins inside that envelope's bounded support (`centroid min/max +/- padding`).
- Skips tiny contributions via a minimum weight cutoff.
- Added binary-search helpers for fast m/z index range selection.

### 3) Simulation pipeline and output throttling

File: `mzLib/Test/FileReadingTests/AnalysisExample.cs`

- Added stage-level stopwatch logging and progress messages.
- Added run profiles (quick vs full-fidelity) selected by env var.
- Replaced pathological `intensityThreshold: 0` during simulated IMSP write with a dynamic threshold:
  - `max(MinImspThreshold, maxSimIntensity * ImspThresholdFraction)`

## New run profiles

In `SimulateJurkatRep2AndWriteImsp`:

- Default mode: `QuickDev`
- Full mode: set `MZLIB_TOPDOWN_SIM_MODE=full`

Profile settings:

- `QuickDev`
  - `MaxRecords = 25`
  - `RtHalfWidth = 0.25`
  - `PointsPerSigma = 1`
  - `MzPaddingInSigmas = 4.0`
  - `ImspThresholdFraction = 1e-4`
  - `MinImspThreshold = 1.0`

- `FullFidelity`
  - `MaxRecords = 200`
  - `RtHalfWidth = 0.40`
  - `PointsPerSigma = 3`
  - `MzPaddingInSigmas = 6.0`
  - `ImspThresholdFraction = 1e-5`
  - `MinImspThreshold = 0.1`

## Validation status

- `TopDownSimulator` test group passed (`23/23`).
- `TopDownEngine` test group previously validated after path updates (`30/30`).
- `AnalysisExample.SimulateJurkatRep2AndWriteImsp` completed successfully in quick mode and produced:
  - `D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.simulated.imsp`

## Operational notes

- If runs feel stalled, use console logger verbosity to surface progress output.
- For rapid iteration, keep default quick mode.
- Use full mode for more complete simulation once parameters look good.
