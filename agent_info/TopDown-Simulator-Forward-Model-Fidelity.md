# Top-Down Simulator — Forward-Model Fidelity for Feature-Finder Benchmarking

**Written:** 2026-07-26
**Branch:** `TopDownSimulator` (HEAD `b5216957`)
**Status:** analysis and proposals. No code changed for this document.

## Purpose and scope

The intended use of the simulated mzML is **benchmarking feature-finding algorithms** — measuring
recall and precision of feature detection against known truth — **not** precise quantification. That
reframing changes which defects matter, and this document is ordered accordingly.

For a feature-finder benchmark, the properties that matter are:

1. **Peak positions** (m/z and RT) — what the finder is trying to locate.
2. **Peak width relative to isotopologue spacing** — this is the *difficulty knob*. Whether adjacent
   isotopologues are resolved or merged is what separates easy cases from hard ones.
3. **Relative intensities within an envelope** — isotope-pattern matching depends on them.
4. **Machine-readable per-feature truth** — without it you cannot score anything.
5. **Noise and background** — sets the false-positive rate, and therefore precision.

Absolute abundance accuracy matters comparatively little: a finder that locates a feature does not
care whether its intensity is 2× off, provided the envelope shape and detectability are right.

⚠️ **This revises the priority order I gave earlier in review.** I previously ranked the residual-metric
defect first, on quant-accuracy grounds. Under a benchmarking goal it drops to third: it is a
diagnostic that is currently misleading, not something that corrupts the emitted data. σ_m rises to
first, and a gap I had not flagged — the ground-truth format — is arguably ahead of both.

---

## Issue 1 — σ_m is a single global constant, so the benchmark is systematically too easy at high m/z

**Severity: high for this use case. It misrepresents difficulty in exactly the regime that
discriminates between feature finders.**

### Root cause

σ_m enters the forward model as a scalar and stays one all the way to the output:

- `ForwardModel` stores a single `_sigmaMz` (`Model/ForwardModel.cs:20,32`) and passes it to every
  `kernel.Evaluate(mz, z, _sigmaMz)` call regardless of m/z (`:51,107`).
- `Simulator.WriteMzml` takes one `sigmaMz` and threads it through to the rasterizer and the ground
  truth sidecar (`Simulation/Simulator.cs:112,126`).
- `AnalysisExample` fits σ per proteoform, then collapses all of them to a single median
  (`Test/FileReadingTests/AnalysisExample.cs:441-446`) which becomes that scalar.

### Why a constant is wrong, quantitatively

On an Orbitrap, resolving power falls as R ∝ (m/z)^−1/2, and R ≡ (m/z)/FWHM, so

```
FWHM ∝ (m/z)^1.5        σ_m = FWHM / 2.3548 = k · (m/z)^1.5
```

The rep2 records cluster around m/z 700–880 and produced a median σ of 0.007212, which implies
k ≈ 3.5 × 10⁻⁷. Extrapolating that one-parameter model across a realistic acquisition range:

| m/z | σ_m implied by physics | σ_m the simulator actually uses |
| --- | --- | --- |
| 600 | 0.0051 | 0.0072 |
| 1000 | 0.0111 | 0.0072 |
| 1500 | 0.0203 | 0.0072 |

So the simulator is ~40% too wide at the bottom of the range and ~3× too narrow at the top.

### Why that specifically breaks a feature-finder benchmark

Isotopologue spacing is `1.00235 / z` m/z. What determines whether a finder can resolve an envelope
is the ratio of that spacing to peak width. Taking a 25 kDa proteoform at z = 25, landing near
m/z 1000, and a larger species reaching m/z 1500:

| case | spacing | FWHM (true σ) | resolved? | FWHM (constant σ) | resolved? |
| --- | --- | --- | --- | --- | --- |
| z=25 at m/z 1000 | 0.040 | 0.026 | yes, comfortably | 0.017 | yes, very comfortably |
| z=25 at m/z 1500 | 0.040 | 0.048 | **no — peaks merge** | 0.017 | yes |

The bottom-right cell is the problem. Real high-mass, high-charge, high-m/z envelopes partially
merge, and handling that is precisely what distinguishes a good deconvolution/feature-finding
algorithm from a naive one. **The current simulation never produces that regime**, so a benchmark
built on it will score all finders too favourably and will not rank them the way real data would.

The converse error at low m/z is milder but real: peaks are simulated wider than they should be.

### Secondary defect — the σ used for fitting is not the σ used for rendering

`ParameterFitter.Fit` fits abundance under that proteoform's *own* σ
(`Fitting/ParameterFitter.cs:43` → `AbundanceFitter.Fit(truth, width.SigmaMz, ...)`), and
`AbundanceFitter` scales the envelope by `peakScale = 1/(σ√2π)` (`Fitting/AbundanceFitter.cs:30`).
The simulation then renders that abundance using the **global median** σ instead. Because the kernel
is normalised to unit area, peak height goes as 1/σ, so every proteoform is rendered at a height
scaled by σ_fitted/σ_median — a factor of up to ~1.7 either way given the observed 0.00378–0.01239
spread. Harmless for feature *detection*, but it means simulated peak heights do not correspond to
the fitted abundances, which will confuse anyone who checks.

### Third defect — the per-proteoform σ estimate is noisy by construction

`EnvelopeWidthFitter.Fit` uses a single (charge, scan) cell — the apex — because
`GroundTruthExtractor` deliberately retains peak windows only for that cell to control memory
(`Extraction/GroundTruthExtractor.cs:161-198`). Each σ therefore rests on a few dozen peaks. Across
25 rep2 records the fitted values span 0.00378 to 0.01239 (3.3×) with no clean m/z trend, which is
the signature of estimator variance rather than physics.

### Proposed fix

Replace the scalar with a width *model*, and fit its one parameter by pooling measurements instead
of taking a median of noisy per-record estimates.

1. **Introduce an abstraction.**
   ```csharp
   public interface IPeakWidthModel { double SigmaAt(double mz); }
   public sealed record ConstantPeakWidth(double SigmaMz) : IPeakWidthModel;      // today's behaviour
   public sealed record OrbitrapPeakWidth(double K) : IPeakWidthModel;            // σ = K · (m/z)^1.5
   ```
   Keeping `ConstantPeakWidth` means every existing test can be migrated mechanically and
   `ForwardModelTests`' bit-identical reference assertion stays meaningful.

2. **Thread it through.** `ForwardModel`, `GridRasterizer`, `Simulator.WriteMzml`,
   `IsotopeEnvelopeKernel.Evaluate`, `AbundanceFitter` and `GlobalAbundanceRefitter` all currently
   take `double sigmaMz`. This is the bulk of the work and it is mostly mechanical. Note
   `ForwardModel.Rasterize` hoists `mzPadding = 6σ` out of its loops (`:72`) and
   `IsotopeEnvelopeKernel.Evaluate` derives its search window from σ — with σ varying by m/z these
   must be computed from the σ at the m/z being evaluated, or conservatively from the maximum σ over
   the span. Getting this wrong silently truncates envelope tails.

3. **Widen what the extractor retains.** To fit k you need width measurements across a *range* of
   m/z. Retaining the apex scan for **every charge** rather than only the apex charge gives
   nCharges × nIso windows instead of nIso — a modest memory increase that yields, for each
   proteoform, several m/z points spanning its charge envelope.

4. **Fit k by pooled weighted regression.** Have the width fitter return the individual
   (m/z, σ²ᵢ, weight) measurements rather than only the pooled scalar. Regress log σ on log(m/z)
   across every measurement from every record. **Fit both slope and intercept, then report the
   slope as a diagnostic** — if it does not come out near 1.5, the data is telling you something
   about the instrument or the estimator, and you want to see that rather than have it assumed away.
   Ship the fit with the slope fixed at 1.5 and only k free.

5. **Use one width model everywhere.** Whatever model is fitted must be the same one passed to the
   simulation, which removes the fitting/rendering mismatch above as a side effect.

6. **Add guard rails.** Reject σ measurements from windows with fewer than N peaks, and clamp the
   final model to a plausible range so one bad record cannot distort k.

### Validation

- Assert the fitted slope is near 1.5 on rep2. If it is not, stop and investigate before shipping.
- Simulate a high-mass proteoform at z ≈ 25 near m/z 1500 and confirm the isotopologue peaks now
  merge — i.e. that the hard regime is reachable at all. This is the whole point of the change.
- Re-run `SimulateJurkatRep2AndWriteMzml` and record peak counts before and after. Wider peaks at
  high m/z will merge centroids on the shared axis (`GridRasterizer.BuildCentroidMzAxis` collapses
  positions closer than σ/100, so the axis itself shifts), so expect the count to move.
- `ForwardModelTests`' reference-implementation equality must still pass under `ConstantPeakWidth`.

### Tests worth adding

1. `OrbitrapPeakWidth.SigmaAt` reproduces the closed form and is monotonically increasing.
2. A forward-model test that peak FWHM at m/z 1500 is measurably greater than at m/z 600 — the
   property a constant σ cannot have.
3. A resolvability test: at fixed charge, an envelope simulated at high m/z has fewer distinct local
   maxima than the same envelope at low m/z.
4. A round-trip test that abundance fitted under a width model and then rendered under the *same*
   model reproduces the expected peak height, pinning the secondary defect closed.

---

## Issue 2 — the ground-truth sidecar records model parameters, not features

**Severity: high for this use case, and not previously flagged. Arguably the actual blocker.**

`Simulator.WriteGroundTruth` (`Simulation/Simulator.cs:122-150`) writes one row per proteoform with
columns `Identifier, MonoisotopicMass, Abundance, RtMu, RtSigma, RtTau, ChargeMu, ChargeSigma,
MinCharge, MaxCharge, SigmaMz`.

That is the *parameter vector that generated* the data. A feature finder emits features — typically
(monoisotopic mass, charge, RT apex, RT start/end, intensity), or (m/z, charge, RT range). To score
recall and precision you must compare like with like, which means a consumer of this file has to
re-derive, from the parameters alone:

- which charge states are actually present above the intensity threshold (a function of the Gaussian
  charge distribution *and* `ScanReductionOptions`, which the sidecar does not record),
- the m/z of each isotopologue at each charge (requires reimplementing the averagine kernel and its
  monoisotopic shift),
- the RT boundaries implied by the EMG profile at some detectability cutoff,
- which scans the feature is above threshold in.

In other words, scoring requires reimplementing the simulator. That is a lot of work per benchmark
consumer and a rich source of disagreement about what counted as a true positive.

### Proposed fix

Emit a **second, feature-level sidecar** alongside the parameter one, written by the simulator that
already knows all of this. One row per (proteoform × charge) that survived reduction:

```
Identifier  MonoisotopicMass  Charge  MonoisotopicMz  ApexMz  ApexRt  RtStart  RtEnd
ApexScanNumber  FirstScanNumber  LastScanNumber  SummedIntensity  ApexIntensity  NumIsotopologues
```

Rules worth fixing explicitly in the writer, because they are the scoring contract:

- RT boundaries defined by where the feature's simulated intensity crosses the same threshold used
  by `SimulatedScanReducer`, so the truth file describes what is actually *in* the mzML rather than
  what the model would have produced without reduction.
- A charge state that reduction removed entirely must not appear — otherwise every benchmark starts
  with unfixable false negatives.
- Scan numbers, not just retention times, so scoring does not depend on float comparison.

Keep the existing parameter file as-is; it is the right artifact for reproducing a run, just not for
scoring one.

---

## Issue 3 — the reported residual fraction double-counts shared signal

**Severity: low for this use case. It is a misleading diagnostic, not a defect in the emitted data.**
Included because it was requested and because the number is currently quoted in run output as though
it meant something specific.

### Root cause

`GlobalAbundanceRefitter.BuildSampleSets` (`Fitting/GlobalAbundanceRefitter.cs:327-372`) builds one
sample set **per proteoform**, drawn from that proteoform's own centroid grid and RT window.
`ComputeResidualEnergyFraction` (`:553-578`) then sums over every sample set independently:

```csharp
for (int p = 0; p < sampleSets.Count; p++)
    for (int k = 0; k < set.Observed.Length; k++)
    {
        observedEnergy += observed * observed;
        residualEnergy += residual * residual;
    }
```

When two proteoforms overlap in m/z and RT, the extractor matches the **same experimental peak** into
both sets, so that peak's intensity enters `observedEnergy` once per proteoform that claims it. With
77 co-eluting models in the rep2 slice the denominator is inflated by an unknown factor.

Consequently the reported `0.3098` is a weighted average over per-proteoform neighbourhoods, where
crowded regions carry more weight — **not** "31% of the file's energy is unexplained." It is not
comparable across runs with different degrees of proteoform overlap, which is exactly the comparison
anyone would want to make with it.

### Why deduplication is not simply "group by m/z"

The same physical peak enters set p at *p's* theoretical centroid and set q at *q's* theoretical
centroid. Those differ by Δmass/z, so the recorded m/z values are not equal even though the source
peak is identical. `BuildSampleSets` stores the theoretical centroid (`mzs.Add(mz)` at `:361`), while
the intensity comes from whichever experimental peak the extractor matched. Any dedup keyed on the
stored m/z will therefore fail to merge them.

Note there are now two `ComputeResidualEnergyFraction` overloads (`:553` jagged, `:580` flat over the
Gauss-Seidel live totals). Both have this property and both need the same fix.

### Proposed fix

1. **Record the source peak identity during extraction.** `GroundTruthExtractor` already knows the
   index of the peak it selected (`idx` in the pass-1 loop, `Extraction/GroundTruthExtractor.cs:143`)
   but keeps only the intensity. Add a parallel array carrying `(zeroBasedScanIndex, peakIndex)` — or
   the matched peak's actual m/z, which is equally unique per scan — to `ProteoformGroundTruth`.
2. **Deduplicate on that key** when computing the residual: accumulate observed and predicted into a
   map keyed by source peak, then compute the fraction once over distinct peaks.
3. **Report the two failure modes separately.** Samples where no peak matched have `observed = 0` but
   generally non-zero prediction. Folding those into one ratio conflates *signal we failed to
   explain* with *signal we predicted that is not there* — opposite problems with opposite fixes.
   Emit them as two numbers.

### Validation

- On a fixture with deliberately non-overlapping proteoforms the deduplicated and current metrics
  must agree, since there is nothing to double-count. That is the sharpest available check that the
  dedup key is right.
- On the rep2 slice, expect the reported fraction to **change**; that is the fix working. Record the
  new value.

### Caveat on the existing monotonicity test

`GlobalAbundanceRefitterTests.ResidualFractionFallsMonotonicallyAndConvergesQuickly` asserts monotone
descent of this metric. `SolveCoordinateUpdate` sums only over proteoform p's own sample set, so the
update is not the exact minimiser of the global objective and monotonicity is not guaranteed by
construction — it held empirically on that fixture. If the metric changes as proposed here, re-check
that this assertion still holds and loosen it to "the final residual is below the initial" if not.
Per direction on 2026-07-26, making the update globally exact is **not** a current priority.

---

## Suggested order

1. **Feature-level ground truth (Issue 2).** Nothing can be scored without it, and it is
   self-contained — no changes to the model, only to what is written out.
2. **Peak width model (Issue 1).** This is what makes the benchmark represent real difficulty. It is
   the largest change, mostly mechanical signature threading, with two places (evaluation windows and
   m/z padding) where care is needed.
3. **Residual metric (Issue 3).** Cheap, and worth doing before anyone tunes a noise model against
   the number, but it does not affect the emitted data.

Noise modelling sits naturally after 1 and 2: a noise term is flexible enough to absorb mis-shaped
peaks, so fitting one while the peak width model is known to be wrong at high m/z would produce a
noise model that looks well-fitted and is compensating for the wrong thing.

## Implementation outcome (2026-07-26)

All three issues are implemented. One finding from the Issue 1 validation changed the shipped
behaviour and is recorded here because it is not something the analysis above anticipated.

### σ_m is not measurable from the rep2 file, and the estimator did not say so

The plan's validation step — "assert the fitted slope is near 1.5 on rep2; if it is not, stop and
investigate" — failed, and investigating it was worth the trouble.

Over the 200-record full-fidelity slice the free slope came out **0.942 ± 0.081**, which is 6.9
standard errors from 1.5. That is not scatter; it is a different law. The cause is that the rep2 raw
is **centroided**, so the "profile peak shape" the width fitter was measuring was really a handful of
isolated centroids sitting inside the isotopologue's midpoint boundaries. Those boundaries are half
the isotopologue spacing, 0.5/z, and since m/z ≈ M/z the boundary width is directly proportional to
m/z — so centroid scatter fills a window proportional to m/z and fits σ ∝ (m/z)^1.0 with high
apparent precision. The estimator was confidently measuring the extraction window.

Three guards now prevent that number from being shipped, each of which the rep2 data trips:

1. **Sampling density** (`EnvelopeWidthFitter.DefaultMinimumSamplesPerSigma = 3`). A window must be
   sampled at least 3 times per σ to count as a peak shape. For k evenly spread centroids the
   implied σ is ≈ gap·√((k²−1)/12), so samples-per-σ is ~1.2 at k = 3 and falls from there — the two
   populations separate cleanly at 3. On rep2 this rejects 4084 of 4097 windows.
   A threshold set too low is worse than none: the test scales with σ, so it admits the *wide*
   windows preferentially. At a threshold of 2 the surviving subset fitted k = 1.26e-6 — σ = 0.018 at
   m/z 600, roughly 4× an Orbitrap — because the selection was biased, not merely noisy.
2. **Lever arm** (`PeakWidthModelFitOptions.MinimumMzRatio = 1.5`) and **evidence**
   (`MinimumMeasurements = 30`). The 13 windows surviving guard 1 spanned m/z 601–725, a ratio of
   1.2. A power law fitted over that band is an extrapolation dressed as a fit.
3. **Slope agreement.** `PeakWidthModelFit.ContradictsWidthLaw` requires the slope to be off by a
   materially wrong amount (>0.25, which is ~17% error in σ at m/z 1500) **and** for that to exceed
   3 standard errors. The absolute term is not redundant: on data lying exactly on the law the
   residuals are rounding error, the standard error collapses toward zero, and a
   standard-error test alone rejects the best possible fit. With the guard at 1 disabled, rep2 trips
   this too — 0.942 ± 0.081 is both 0.56 off and 6.9 standard errors out.

An earlier revision of the standard error divided by the sum of the weights rather than the number
of measurements, treating the profile samples inside one window as independent replicates. That
understated it by ≈√(samples per window) and made the accept/reject decision a function of the
detector's sampling rate. It is fixed, and pinned by
`StandardErrorDoesNotDependOnHowFinelyEachWindowWasSampled`. The figure quoted above is the
corrected one; the uncorrected formula reported 10.6 standard errors.

With all three in place the pipeline refuses on rep2, says why, falls back to a constant width, and
points at `MZLIB_TOPDOWN_SIM_PEAK_WIDTH_K` for supplying k directly. **The consequence for the
benchmark is that this data set cannot supply its own width law.** Fitting one needs a profile-mode
acquisition; otherwise k has to come from the instrument's known resolving power (R = 60,000 at
m/z 400 gives k ≈ 3.5e-7).

### Numbers from the validation runs

Supplying k = 3.5e-7 reproduces the physics table in Issue 1 essentially exactly — σ = 0.00514 /
0.01107 / 0.02033 at m/z 600 / 1000 / 1500 against the predicted 0.0051 / 0.0111 / 0.0203 — which is
the end-to-end confirmation that the width model is threaded through correctly.

Full-fidelity rep2 (200 models, 2906 scans):

| width model | peaks | features |
| --- | --- | --- |
| Constant σ = 0.012 | 1,980,950 | 1,848 |
| Orbitrap k = 3.5e-7 | 1,901,912 | 1,861 |

The ~4% drop in peak count is the predicted effect of wider high-m/z peaks merging centroids on the
shared axis.

### Issue 3 note

`ResidualFractionFallsMonotonicallyAndConvergesQuickly` still passes unmodified after the metric was
deduplicated, so the loosening the plan allowed for was not needed.

`OverpredictedFraction` is **not** deduplicated, and its documentation now says so. An unmatched
sample has no peak index to key on, so two proteoforms probing the same empty neighbourhood each
contribute. It remains a within-run diagnostic and is not comparable across runs of differing
proteoform density — the same caveat Issue 3 removed from the unexplained fraction. Fixing it needs
a different key (e.g. scan plus quantized m/z) and is not done.

### Known limitations carried forward

- **`AbundanceFitter` is not numerically identical on the `ConstantPeakWidth` path.** It now
  predicts with `kernel.Evaluate` (which carries neighbouring isotopologue tails) rather than
  `Intensity(i)/(σ√2π)`. This is deliberate — it is what makes an abundance fitted under a width
  model reproduce the height rendered under it, and it matches what `GlobalAbundanceRefitter`
  already did — but `Abundance`, `Residual` and `SamplesUsed` all move for pre-existing constant-σ
  fixtures. The bit-identity claim applies to `IsotopeEnvelopeKernel.Evaluate` and
  `ForwardModel.Rasterize`, not to the fitter.
- **`IsotopeEnvelopeKernel` is not thread-safe** and is now documented as such. No current call site
  shares one across threads, but `ForwardModel` holds kernels as instance state, so parallelising
  `Rasterize` would corrupt the lazy caches.
- **The dedup key carries no extractor identity.** Truths extracted from two different
  `GroundTruthExtractor` instances would alias if passed to one `Refit` call. No caller does this.

### Bugs found in review

Two agent reviews over this change set and the two commits beneath it turned up four defects that
were fixed, all of which the 107 tests passing at the time missed:

1. **`FeatureGroundTruth` read the monoisotopic m/z off `CentroidMzs(z)[0]`.** Averagine drops
   isotopologues below an absolute probability of 1e-8, and measurement shows the monoisotopic peak
   is retained through 31 kDa and gone from 32 kDa — the offset reaches 1.0 Da at 32 kDa, 2.0 at
   40 kDa and 5.0 at 50 kDa. Above that threshold the column reported a heavier isotopologue's m/z
   under a header saying monoisotopic, for exactly the large proteoforms top-down cares about. Now
   computed from the mass. Pinned by
   `MonoisotopicMzStaysCorrectAboveTheMassWhereAveragineDropsTheMonoisotopicPeak`.
2. **The slope standard error treated within-window samples as replicates** (see above).
3. **The outlier pass readmitted the measurements it had just rejected** when too few survived,
   reporting "n/n used" — indistinguishable from finding no outliers. It now refuses the fit.
4. **The m/z span guard ran before outlier rejection** and reported the pre-rejection span, so
   far-out measurements could satisfy it and then be discarded, leaving the narrow-band
   extrapolation the guard exists to refuse. Both now operate on the surviving measurements.

Also addressed: `CentroidMzs` returned its internal cache by reference while `Evaluate`
binary-searches it (a caller sorting the array in place would silently corrupt evaluation — the
clone that used to absorb this was removed in `c461b6eb`); `SourcePeakMzs` became a `float` delta
array, halving that tensor's footprint in a class that otherwise avoids dense tensors; and
`BuildMzGrid` no longer drops the top of its padding.
