# Top-Down Simulator — Known Issues and Deferred Work

**Written:** 2026-07-26
**Branch:** `TopDownSimulator`
**Context:** raised by two code reviews over the forward-model-fidelity change set
(`agent_info/TopDown-Simulator-Forward-Model-Fidelity.md`) and the two commits beneath it,
`b5216957` (Gauss-Seidel abundance refit) and `c461b6eb` (isotopologue ordering).

Everything here was **verified against the code, not just reported**. Items the reviews raised that
were fixed are listed in the "Fixed" section of the fidelity document and are not repeated.

Ordered by what I would act on first.

---

## 1. `BuildTranspose` allocates outside the memory cap it exists under

**Where:** `mzLib/TopDownSimulator/Fitting/GlobalAbundanceRefitter.cs`, `BuildTranspose` and its call
site in `Refit`.
**Severity:** high on large runs; this is the one I would fix first.

`TryBuildBasisMatrices` budgets carefully against `GlobalAbundanceRefitOptions.MaxBasisCacheBytes`
(default 1 GiB) and degrades gracefully — it returns `null` and the sweep falls back to recomputing
basis values on demand. `BuildTranspose` then unconditionally allocates `new int[nnz]` plus
`new double[nnz]` over the *same* entry count, with **no budget check and no fallback**, plus
`liveTotals` at 8 bytes per sample across all sets.

So a basis that just fits under the cap peaks at roughly twice it, and there is no
`OutOfMemoryException` path back to the uncached branch. The budget also ignores the per-set
`rowStart` arrays and the `List<int>`/`List<double>` backing capacity during construction, which can
add roughly another 2× of the largest single sample set.

**Why it matters here:** the full-fidelity rep2 run fits 200 models and
`MZLIB_TOPDOWN_SIM_GLOBAL_REFIT_MAX_MODELS` defaults to 200, so this path is exercised by the
headline workload.

**Fix sketch:** fold the transpose and `liveTotals` into the same budget — the entry count is known
once the matrices are built, so the check is `2 * nnz * bytesPerEntry + totalSamples * 8 <= maxBytes`
before calling `BuildTranspose`, falling back to the uncached branch when it fails. Cheap, but it
changes committed behaviour and deserves a test that actually sets `MaxBasisCacheBytes` small —
today nothing does, so the whole fallback branch has zero coverage.

## 2. The coordinate update does not minimise the residual it is scored against

**Where:** `GlobalAbundanceRefitter.SolveCoordinateUpdate`, and the residual reported by
`ComputeResidualEnergy`.
**Severity:** medium. Affects a diagnostic and a test's justification, not the emitted mzML.

`SolveCoordinateUpdate` accumulates its numerator and denominator **only over proteoform `p`'s own
sample set**. But model `p` has nonzero basis at samples belonging to other sample sets — that is
precisely why the sparse matrix for set `q` contains model `p`, and why `ApplyAbundanceDelta`
publishes `p`'s delta into other sets' totals.

Publishing the delta globally (what `b5216957` changed) makes each coordinate *see* its neighbours'
updates, but does not make the update *account for* its own effect on their samples. The fixed point
is the non-symmetric system "`b_p` orthogonal to the residual, restricted to rows of set `p`" — not
the normal equations of any least-squares problem. Gauss-Seidel on that carries no descent or
convergence guarantee.

**Consequence:** `GlobalAbundanceRefitterTests.ResidualFractionFallsMonotonicallyAndConvergesQuickly`
asserts monotone descent, which holds **empirically on that fixture**, not by construction. The
commit message's claim that "monotone descent is the property Gauss-Seidel buys" is wrong for this
formulation.

**Deliberately deferred:** per direction on 2026-07-26, making the update globally exact is not a
priority. Recorded so nobody later treats the monotonicity test as proof of a property the algorithm
has.

## 3. `ResidualFractionByIteration` differs between the cached and uncached paths

**Where:** `GlobalAbundanceRefitter.Refit`, the per-sweep residual.
**Severity:** medium-low. A latent test flake and a documentation inaccuracy.

The cached path reads `liveTotals` after a sweep's worth of incremental deltas; the uncached path
recomputes from scratch. The two therefore differ by more than the "same arithmetic" the class docs
claim, and `ResidualFractionByIteration[^1] != FinalResidualFraction` (the final value is always a
fresh recompute) even when nothing moved.

**Flake scenario:** feed `Refit` proteoforms already at the optimum. `InitialResidualFraction` is R;
the first sweep refreshes to exactly R, applies ~0 deltas, and reports R ± rounding. The assertion
`ResidualFractionByIteration[0] < InitialResidualFraction` then fails on a coin flip. The current
fixtures start far from the optimum, so this does not bite today.

## 4. The uncached fallback now costs two dense passes per sweep

**Where:** `GlobalAbundanceRefitter.Refit`, the `transpose is null` branch.
**Severity:** low, but it contradicts its own comment.

The comment on the per-coordinate `ComputeSampleSetTotals` says it "costs no more than the single
up-front pass it replaces", which is true in isolation — summed over `p` it is one pass. But the
sweep then adds a **full** `ComputePredictedTotals` purely to populate `ResidualFractionByIteration`.
Net: the branch used precisely when the problem is too big to cache went from one dense pass per
sweep to two, for diagnostics. Skipping or sampling the per-sweep diagnostic on that branch would
recover it.

## 5. `OverpredictedFraction` is not deduplicated

**Where:** `GlobalAbundanceRefitter.ComputeResidualEnergy`.
**Severity:** medium-low. Documented accurately as of this change set; not fixed.

`UnexplainedFraction` is computed over **distinct experimental peaks**, keyed on
`(zeroBasedScanIndex, peakIndex)`, which is the whole point of Issue 3 in the fidelity plan.
`OverpredictedFraction` is not: a sample that matched no peak has no peak index to key on, so two
proteoforms probing the same empty neighbourhood each contribute their own `predicted²`. The
numerator therefore still scales with proteoform crowding — the same incomparability defect, moved
to the second number. `UnmatchedSamples` is likewise a sample count, not a count of distinct empty
positions.

The XML docs on `ResidualEnergySummary` now state this explicitly, so it is a known limitation rather
than a false claim. It should be read as a within-run diagnostic and **not compared across runs of
differing proteoform density**.

**Fix sketch:** key unmatched samples on `(zeroBasedScanIndex, quantized m/z)` — quantized at
something principled like the extractor's `MzWindowHalfWidth`. Needs care at quantization
boundaries.

## 6. `AbundanceFitter` is not numerically identical on the `ConstantPeakWidth` path

**Where:** `mzLib/TopDownSimulator/Fitting/AbundanceFitter.cs`.
**Severity:** low. Deliberate, but worth knowing before comparing against old numbers.

The prediction is now `kernel.Evaluate(centroid, z, widthModel)` — which carries the tails of
neighbouring isotopologues — rather than `kernel.Intensity(i) / (σ√2π)`. This is what makes an
abundance fitted under a width model reproduce the peak height rendered under that same model, and it
matches what `GlobalAbundanceRefitter.EvaluateUnitContribution` already did.

But it means `Abundance`, `Residual` and `SamplesUsed` all move for pre-existing constant-σ fixtures.
`SamplesUsed` moves most: `Evaluate` is > 0 wherever any isotopologue is in window, whereas
`Intensity(i)` could be exactly 0 for a negligible isotopologue.

**The bit-identity claim for `ConstantPeakWidth` applies to `IsotopeEnvelopeKernel.Evaluate` and
`ForwardModel.Rasterize`, not to the fitter.**

## 7. `IsotopeEnvelopeKernel` is not thread-safe

**Where:** `mzLib/TopDownSimulator/Model/IsotopeEnvelopeKernel.cs`.
**Severity:** low, latent. Now documented on the type.

Centroid positions and per-width-model Gaussian constants are cached lazily into plain dictionaries
on every `Evaluate`. No current call site shares a kernel across threads —
`AnalysisExample`'s `Parallel.For` builds fresh fitters per iteration and `GroundTruthExtractor.Extract`
builds its own kernel — but `ForwardModel` holds kernels as instance state, so **parallelising
`ForwardModel.Rasterize`, or sharing a `ForwardModel` across threads, would corrupt the caches.**

Related: the width cache keys on `IPeakWidthModel` value equality. A custom `IPeakWidthModel`
implemented as a plain class (reference equality) would flush the cache on every `Evaluate` — correct
results, but O(nIso) rebuild per evaluated point.

## 8. The dedup key carries no extractor identity

**Where:** `GlobalAbundanceRefitter.BuildSampleSets`.
**Severity:** very low; no caller does this.

The key is `((long)zeroBasedScanIndex << 32) | (uint)peakIndex`. That index is a position within one
`GroundTruthExtractor`'s MS1 scan list. Truths extracted from **two different extractor instances**
(i.e. two different files) passed to a single `Refit` call would alias silently and be deduplicated
as if they were the same peak. Every current caller builds one extractor per run.

---

## Non-obvious library facts worth keeping

These are properties of mzLib itself, discovered by measurement during this work. They are not
specific to the simulator.

### Averagine drops the monoisotopic peak above ~31 kDa

`AverageResidue.MinRes = 1e-8` is an **absolute** probability cutoff, and
`IsotopicDistribution.CalculateFineGrain` drops terms with `prob <= minProbability`. The all-light
composition's absolute probability falls below that around 31 kDa, so above it the monoisotopic peak
is simply **absent from the averagine tables**.

Measured via `IsotopeEnvelopeKernel`, offset of the lowest retained isotopologue from the
monoisotopic mass:

| mass (Da) | offset | note |
| --- | --- | --- |
| 5,000 – 31,000 | 0.0000 | monoisotopic retained |
| 32,000 | +1.0029 | first isotopologue lost |
| 35,000 | +1.0034 | |
| 40,000 | +2.0068 | |
| 50,000 | +5.0151 | |

**Consequence:** `IsotopeEnvelopeKernel.CentroidMzs(z)[0]` is the lowest *surviving* isotopologue,
**not** the monoisotopic one, above ~31 kDa. Code needing the monoisotopic m/z must compute it from
the mass (`mass.ToMz(z)`). This caused a real bug in `FeatureGroundTruth`, now fixed and pinned by
`MonoisotopicMzStaysCorrectAboveTheMassWhereAveragineDropsTheMonoisotopicPeak`.

The envelope's *anchoring* is unaffected — the shift is derived from `Averagine.DiffToMonoisotopic`,
not from array position.

### σ_m is not measurable from a centroided acquisition

Peak width cannot be recovered from a file whose peaks have already been reduced to positions. Any
estimator run over such a file measures its own extraction window instead: the window is half the
isotopologue spacing, `0.5/z`, and since m/z ≈ M/z, centroid scatter fills a window **proportional to
m/z** and fits σ ∝ (m/z)^1.0 with high apparent precision.

`D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw` is such a file. Over 200 records, 4084 of 4097
windows are too coarsely sampled to be peak shapes, and the unguarded fit returns σ ∝ (m/z)^0.942 ±
0.081. Fitting a width law from this data set requires a **profile-mode acquisition**; otherwise `k`
must come from the instrument's known resolving power (R = 60,000 at m/z 400 gives k ≈ 3.5e-7, which
reproduces the expected σ at m/z 600/1000/1500 essentially exactly).

Reproduce the artifact with `MZLIB_TOPDOWN_SIM_MIN_SAMPLES_PER_SIGMA=0`, which disables the guard.

---

## Test-suite baseline (2026-07-26)

The full suite reports **86 failures out of 3028**, all environmental on this machine and unrelated to
the simulator:

- 51 cascade from a `OneTimeSetUp` calling `UsefulProteomicsDatabases.Loaders.DownloadContent` —
  UniProt/PTM-list downloads with no network. These present as protein-variant and decoy-database
  failures, which is misleading until you read the stack.
- 31 direct network failures.
- 3 missing files, 1 file lock, on `D:\`.

**Zero** are in any simulator-related test. `Test.TopDownSimulator` is 113/113 green. Worth knowing
before chasing a "regression" that is really an offline machine.
