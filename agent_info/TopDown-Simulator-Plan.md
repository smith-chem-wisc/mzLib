# Top-Down MS1 Simulator — Project Plan

**Status:** planning complete, Phase 1 starting
**Location:** `mzLib/TopDownSimulator/` (shipped via the mzLib NuGet package; a thin CLI tool may wrap it later)
**Goal:** A parametric forward model for top-down MS1 data that can (a) reproduce observed raw data when fit to a known proteoform ID, and (b) generate realistic synthetic MS1 for unseen proteoforms using population priors. Output deliverables are a library API, synthetic `.mzML` writers, a short technical note, and a published tool.

---

## 1. Forward model

For scan time `t_s` and m/z bin `b`, the observed MS1 intensity is a sum over proteoforms `p` and charges `z`:

```
I(s, b) = Σ_p A_p · g_p(t_s) · Σ_z f_p(z) · φ(b; M_p, z, σ_m)  +  ε(s, b)
```

| Factor | Symbol | Parametric form | Fit source |
|---|---|---|---|
| Abundance | `A_p` | scalar (per-proteoform) | NNLS on residual once shape params fixed |
| RT profile | `g_p(t)` | EMG(μ_rt, σ_rt, τ) | summed-charge XIC in the ID's RT neighborhood |
| Charge distribution | `f_p(z)` | Gaussian on `z` (μ_z, σ_z); upgrade to skew-normal if asymmetric residuals | per-charge apex intensities |
| Isotope envelope | `φ(b; M, z, σ_m)` | averagine via `Chemistry.IsotopicDistribution`, convolved with Gaussian of width σ_m | theoretical + σ_m from observed FWHM |
| Noise | `ε` | baseline + Gaussian + Poisson shot | empty-region fit around the envelope |

- `M_p` is fixed from the MetaMorpheus ID.
- `σ_m` is instrument-shared — Phase 1 uses a constant; Phase 2 upgrades to the Orbitrap-style `σ_m ∝ m/z^1.5 / R`.
- All shape params fit sequentially (`σ_m` → `g_p` → `f_p(z)` → `A_p`). Start linear-space for abundance, log-space for shape fits.

---

## 2. Project layout

New project at `mzLib/TopDownSimulator/`, sibling to `TopDownEngine`. Ships in the mzLib NuGet.

```
mzLib/TopDownSimulator/
├── TopDownSimulator.csproj       # net8.0, refs Chemistry, MassSpectrometry,
│                                 # Readers, UsefulProteomicsDatabases, MzLibUtil,
│                                 # SpectralAveraging
├── Model/
│   ├── ProteoformModel.cs        # θ_p = {A, μ_z, σ_z, μ_rt, σ_rt, τ, M}
│   ├── ForwardModel.cs           # evaluate I(s, b) on a scan × m/z grid
│   ├── ChargeStateDistribution.cs# interface + Gaussian/SkewNormal impls
│   ├── EmgProfile.cs             # exponentially-modified Gaussian
│   ├── IsotopeEnvelopeKernel.cs  # averagine → φ, convolved with σ_m
│   └── NoiseModel.cs             # baseline + Gaussian + Poisson
├── Extraction/
│   ├── MmResultLoader.cs         # parse MM PSM table (reuse Readers parser)
│   ├── CoeluterFinder.cs         # MM IDs within an (RT, M) neighborhood
│   └── GroundTruthExtractor.cs   # pull scan × m/z tensor via PeakIndexingEngine
├── Fitting/
│   ├── ParameterFitter.cs        # top-level: (raw, MM id) → θ_p
│   ├── RtProfileFitter.cs        # EMG fit on the summed-charge XIC
│   ├── ChargeDistributionFitter.cs
│   ├── EnvelopeWidthFitter.cs    # σ_m calibration
│   └── AbundanceFitter.cs        # NNLS for A_p (per-proteoform or joint)
├── Simulation/
│   ├── Simulator.cs              # θ_p list → synthetic MsDataFile
│   ├── GridRasterizer.cs         # evaluate forward model on output grid
│   └── ScanBuilder.cs            # build MsDataScan objects
└── Comparison/
    ├── SpectralAngle.cs          # per-scan cosine similarity
    ├── XicCorrelation.cs         # per-charge Pearson
    ├── ResidualAnalyzer.cs       # model-unexplained energy
    └── ComparisonReport.cs       # TSV/JSON + plot-ready arrays
```

Tests at `mzLib/Test/TopDownSimulator/` mirroring the subfolders.

**Mzlib integration notes:**
- Public surface goes into the NuGet. Treat every public API as a potential MetaMorpheus breaking change (integration CI).
- Keep fitting helpers `internal` unless they need to be called from outside.
- Do not add to `mzLib.nuspec` until the API has settled enough that a cross-repo consumer would be reasonable.

---

## 3. Phased plan

### Phase 1 — Single-proteoform replay (MVP)
Pick one MetaMorpheus ID, pull the ground-truth scan × m/z tensor, fit θ_p sequentially (`σ_m` → `g_p` → `f_p(z)` → `A_p`), regenerate synthetic MS1 on the same scan grid, and report per-scan cosine + per-charge XIC correlation. Validates the isotope kernel, EMG, and peak-width model end-to-end.

### Phase 2 — Charge distribution calibration
Start Gaussian-on-`z`. If residuals are systematically asymmetric across many proteoforms, upgrade to skew-normal. Add the resolution-scaled `σ_m(m/z)` model here.

### Phase 3 — Co-eluting proteoforms
`CoeluterFinder` pulls MM IDs within an RT and `M` neighborhood of the anchor. Fit each independently, then joint-NNLS on abundances against the combined raw tensor. Compare the multi-proteoform forward model to the raw window.

### Phase 4 — Noise model + population priors
Fit the noise model in empty neighborhoods around envelopes. Run the fitter across *all* MM results in the file, characterize distributions over `(μ_z, σ_z, σ_rt, τ, σ_m, noise params)`. **These priors are the generative payload of the tech note** — the model you can sample from for unseen proteoforms.

### Phase 5 — Synthetic mzML + harness integration
`Simulator → MsDataFile → MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra` plus a ground-truth sidecar (TSV with θ_p per proteoform). Feeds directly into the TopDownEngine test harness (task T20.1).

---

## 4. Metrics for the technical note

- Per-scan cosine similarity (simulated vs. real MS1) — histogram across all fitted proteoforms
- Per-charge XIC Pearson correlation
- Residual-energy fraction `||y − ŷ||² / ||y||²`
- Bland–Altman plot on per-charge intensities
- Prior distributions from Phase 4 as the generative deliverable

---

## 5. Decisions (locked)

1. **Parametric forward model**, not empirical re-synthesis.
2. **Sequential per-factor fitting** (`σ_m` → `g_p` → `f_p(z)` → `A_p`) — more debuggable than joint nonlinear least squares. Can upgrade to joint LM later if needed.
3. **Phase 1 uses a constant `σ_m`**; Phase 2 upgrades to `σ_m(m/z)` for Orbitrap.
4. **Log-intensity for shape fits, linear for abundance** — log handles dynamic range for shapes, linear keeps NNLS clean.
5. **Start Gaussian for `f_p(z)`**; revisit skew-normal only if residuals warrant.
6. **Project lives inside mzLib, ships via the mzLib NuGet.** A thin CLI wrapper may follow if users want it standalone.

## 6. Open questions deferred to each phase

- Exact co-eluter window width (Phase 3) — empirical from data, not upfront.
- Whether to fit noise parameters per-scan or per-file-region (Phase 4).
- Whether the simulator needs to synthesize MS2 (later; current focus is MS1 only).
