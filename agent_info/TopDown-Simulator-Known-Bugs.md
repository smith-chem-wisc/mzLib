# Top-Down Simulator — Known Bugs and Fix Plan

**Written:** 2026-07-26
**Updated:** 2026-07-26 — **both bugs are now fixed.** See "Outcome" under each.
**Branch:** `TopDownSimulator` (bugs confirmed at `bf36b4bc`)

Build and test from `mzLib/`, not the repo root:

```
cd mzLib
dotnet build ./Test/Test.csproj
dotnet test  ./Test/Test.csproj --filter "FullyQualifiedName~TopDownSimulator"
```

Baseline before the fixes was **78 passed, 0 failed** under that filter. It is now **82 passed, 0 failed**
(four new tests, described below).

⚠️ **Build gotcha discovered while validating this work.** `Copy-Item` on Windows preserves the source
file's last-write time. Restoring a file that way can leave it *older* than the compiled output, so
MSBuild skips recompiling it and `dotnet test --no-build` silently runs stale code. Two measurement
runs during this work reported "no change" for exactly that reason. If a change appears to have no
effect, check `(Get-Item <file>).LastWriteTime` against the DLL before believing it.

---

## Bug 1 — `CentroidMzs` returned isotopologues in intensity order, and two callers assumed m/z order

**Severity: high. Silently wrong rather than merely slow, and live on the current dataset.**

### Root cause

`Averagine`'s static constructor sorts each theoretical envelope by *intensity* and reverses it, so
`AllMasses[i]` is ordered **most-intense first**, not by mass:

- `mzLib/MassSpectrometry/Deconvolution/AverageResidue/Averagine.cs:45-47`
  ```csharp
  Array.Sort(intensities, masses);
  Array.Reverse(intensities);
  Array.Reverse(masses);
  ```

That ordering is deliberate and **was not changed** — the very next lines depend on it:

```csharp
MostIntenseMasses[i] = masses[0];
DiffToMonoisotopic[i] = masses[0] - chemicalFormulaReg.MonoisotopicMass;
```

`IsotopeEnvelopeKernel` copied that ordering into `_neutralMasses` / `_normalizedIntensities`, so
`CentroidMzs(z)`, `NeutralMass(i)` and `Intensity(i)` were all intensity-descending.

### The two broken callers

**1. `EnvelopeWidthFitter` — live, corrupted fitted σ_m**
`mzLib/TopDownSimulator/Fitting/EnvelopeWidthFitter.cs:74-80` computes midpoint boundaries between
isotopologue `i` and its index-neighbours `i-1` / `i+1` to stop a neighbouring isotopologue's peaks
leaking into this isotopologue's second-moment estimate. With intensity ordering those neighbours
were arbitrary isotopologues, so the boundaries were arbitrary — frequently *inverted*
(`rightBoundary < leftBoundary`), which rejects every sample and drops the isotopologue entirely.

**2. `GridRasterizer.BuildMzGrid` — dormant but wrong**
`mzLib/TopDownSimulator/Simulation/GridRasterizer.cs:114-119` took `centroids[0]` and `centroids[^1]`
as the m/z extremes; they were the most- and least-intense isotopologues.

### Fix applied — option A, sort inside the kernel

`mzLib/TopDownSimulator/Model/IsotopeEnvelopeKernel.cs`:

1. A single paired `Array.Sort(_neutralMasses, _normalizedIntensities)` at the end of the constructor
   puts both arrays in ascending-mass order. `Averagine`'s static tables are untouched, so nothing
   outside `TopDownSimulator` is affected.
2. The monoisotopic-shift computation still reads `theorMasses[0]` from `Averagine`'s array *before*
   the local arrays are reordered, so it keeps working. Both facts now carry comments saying so.
3. The class remark states the contract: index 0 is the lowest-mass isotopologue, not the most
   intense.
4. `GetSortedEvaluationArrays` and its per-charge cache are gone — `Evaluate` reads `CentroidMzs`
   and `_normalizedIntensities` directly. `GetMzBounds` takes the ends instead of scanning, and its
   cache is gone too. Both callers above became correct without being edited.

Two of the three mutable dictionaries on the read path are therefore gone; the remaining
`_centroidCacheByCharge` still mutates on read, so **the "do not share kernels across threads" rule
still stands.**

### Outcome — measured

The regression test below recovers a known σ_m of **0.004** from a synthetic z=20 profile envelope.
Against the old intensity-ordered kernel it returned **0.017836** — 4.5x too wide.

On real data (`SimulateJurkatRep2AndWriteMzml`, quick mode, 25 records):

| | before | after |
| --- | --- | --- |
| median fitted σ_m | 0.007869 | **0.007212** |

Individual records moved more than the median suggests — e.g. 0.005535 → 0.004833, 0.009231 →
0.007955, 0.004222 → 0.004895 — and `PeaksUsed` rose across the board (11 → 17, 16 → 22, 15 → 29),
which is the inverted-boundary case no longer discarding whole isotopologues.

### Tests

Updated:

- `ForwardModelTests.cs` — `CentroidMzs(8)[0]` was being used as `apexMz`. It is now the monoisotopic
  peak, which for a 10 kDa proteoform is ~0.4% of the envelope apex and would have broken the
  assertions. Added a `MostIntenseMz` helper that searches for the brightest isotopologue, preserving
  the tests' original intent.
- `MzmlExportTests.cs:259,261` — the comment calling `CentroidMzs(8)[0]` "the monoisotopic m/z" was
  wrong before and is right now. Added a note saying why.

New:

- `IsotopeEnvelopeKernelTests.IsotopologuesAreOrderedByAscendingMassAtEveryCharge` — strict ascent of
  masses and m/z over four masses and four charges, plus `GetMzBounds` equalling the ends.
- `IsotopeEnvelopeKernelTests.ReorderingKeepsEachMassPairedWithItsOwnIntensity` — rebuilds the
  expected (mass, weight) pairs from `Averagine`'s tables and checks the kernel against them, so the
  pairing is verified against the source rather than against itself. Also pins that the brightest
  isotopologue is still `MostIntenseMasses[idx] + shift` and is no longer at index 0.
- `EnvelopeWidthFitterTests.ProfileInputRecoversSigmaWhenNeighbouringIsotopologuesShareTheWindow` —
  the regression test for the real defect. z=20, where isotopologue spacing (~0.050 m/z) is inside
  the extractor's ±0.05 window. **Verified to fail before the fix and pass after.**

The pre-existing profile fixture was refactored to share a `BuildProfileScans` helper with the new
one; its behaviour is unchanged.

---

## Bug 2 — the abundance refit used stale totals within a sweep

**Severity: medium. Wrong algorithm, not wrong output shape — it converged to a worse answer, more slowly.**

### Root cause

`mzLib/TopDownSimulator/Fitting/GlobalAbundanceRefitter.cs` computed `predictedTotals` once before
the model loop but mutated abundances inside it, so coordinate `p` reasoned about totals that still
reflected the *old* abundances of models `0 … p-1`. That is a Jacobi sweep, despite the class
documenting itself as coordinate descent. True coordinate descent is Gauss-Seidel: each update sees
every earlier update in the same sweep.

### Fix applied — maintain the residual incrementally

1. New `BasisTranspose` record: the transpose of the existing `SparseBasis`, grouping entries by
   model instead of by sample, with samples addressed by a flat index (`SampleSetOffsets[p] + k`) so
   one `liveTotals` array covers every sample set at once.
2. It is built **from the entries of the forward matrices**, not re-derived, so the two cannot
   disagree about which terms survived `BasisSparsityThreshold`.
3. Each sweep starts by rebuilding `liveTotals` from the current abundances (O(nnz)), so rounding
   error from incremental updates cannot accumulate across sweeps. Filling model-major in ascending
   model order reproduces the forward pass's summation order exactly.
4. After each coordinate update, `ApplyAbundanceDelta` publishes `delta * basis` to every sample that
   model touches. Cost per sweep stays O(nnz).
5. The uncached fallback (basis exceeds `MaxBasisCacheBytes`) recomputes one sample set's totals
   immediately before use, which is Gauss-Seidel too and costs the same per sweep as the single
   up-front pass it replaces.

`GlobalAbundanceRefitResult` gained `ResidualFractionByIteration` so convergence behaviour is
observable; `Verbose` output now prints the residual alongside the max relative change.

### Outcome — measured

On a synthetic fixture of four heavily overlapping proteoforms (masses within ~1 Da, same elution
window, same charge envelope):

| tolerance | Jacobi (before) | Gauss-Seidel (after) |
| --- | --- | --- |
| 1e-3 | converged in **45** sweeps | converged in **10** |
| 1e-5 | **not converged** after 60 | converged in **16** |

On the real rep2 31–35 min slice (77 models), with Bug 1 also fixed:

| | before | after |
| --- | --- | --- |
| residual energy fraction | 2.93 → 1.19 | **1.7157 → 0.3098** |
| iterations / converged | 8 / no | 8 / no |

### ⚠️ This changes fitted abundances

It is a **behavioural change**, not a speedup, and should be committed as one.

### Tests

- `GlobalAbundanceRefitterTests.ResidualFractionFallsMonotonicallyAndConvergesQuickly` — asserts the
  residual falls monotonically across sweeps (the property Gauss-Seidel buys and Jacobi does not
  guarantee), that it converges at 1e-5 well within the 60-sweep budget, and that all four
  abundances land within 1e-4 relative of truth.
- The three existing tests still pass unchanged, including
  `SparseBasisMatchesTheExactBasisToFarBetterThanFitPrecision`, which is the guard that the sparsity
  cutoff is applied identically on both paths.

---

## Open concern — the refit residual — largely resolved

The original concern was that the slice refit ended at a residual energy fraction of **1.19**, i.e.
`||y − ŷ||² > ||y||²`: the model mispredicting by more than the entire signal it is trying to explain.

With both bugs fixed that number is **0.3098**, comfortably below 1, and the starting point dropped
from 2.93 to 1.72 as well. Bug 1 explains the better start (a corrupted σ_m made every envelope the
wrong width, so no set of abundances could fit); Bug 2 explains the better finish.

Two things remain worth noting, neither a bug in the two fixed here:

1. **The slice still reports `converged: false` after 8 sweeps**, but the residual has effectively
   plateaued (0.31034 → 0.31003 → 0.30981 over the last three). The convergence *criterion* is what
   is failing to trip, not the fit: `maxRelativeChange` is set to a hard `1` whenever a model sitting
   at zero abundance takes any positive value, so a single model flickering on and off pins the
   metric at 1 forever. A relative-change criterion that tolerates that, or an absolute-change floor,
   would let this terminate honestly. Raising `MaxIterations` will not help.
2. **No noise model.** `Model/NoiseModel.cs` was never written, so baseline and chemical background
   remain structurally unexplainable. That is the likeliest source of the residual 0.31.

Also still open, and unchanged by this work: the refit evaluates every model across the *global*
charge span rather than each model's own fitted range.

---

## Context the next agent needs — do not redo these

- **Output is centroided mzML only.** IMSP is retired. mzLib's mzML reader throws on profile-mode
  files. See `agent_info/TopDown-Simulator-Status.md`.
- **Centroids are generated directly** at theoretical isotopologue m/z
  (`GridRasterizer.RasterizeAtCentroids`). There is no profile-grid-then-peak-pick step any more.
- **The fitting path was optimized in `ea8196a8` and `bf36b4bc`.** `ForwardModelTests` pins the
  `Rasterize` optimization as bit-identical against a reference implementation — **keep that test
  passing**; it is the guard on that rewrite.
- **Do not cache `IsotopeEnvelopeKernel` instances across threads.** Still true after the Bug 1 fix:
  `_centroidCacheByCharge` is a plain `Dictionary` mutated on the read path. Kernels are deliberately
  constructed per record so the parallel fit loop stays safe.
- **The 200-model cap on the global refit is still needed.** *Building* the basis still evaluates
  every model at every sample, so it remains quadratic in model count. The full q≤0.01 run (1087
  models) still skips the refit entirely. Lifting the cap needs spatial pruning at build time — only
  consider models whose RT and m/z ranges overlap the sample.
- **Three tests fail on this machine for unrelated reasons**: `ControlXIC`, `Ms1Example` and
  `Ms2LogExample` hardcode `D:\Human_Ecoli_TwoProteome_60minGradient\...`, which does not exist here,
  and call Plotly `.Show()`. They are not `[Explicit]` but they are not runnable without that data.
