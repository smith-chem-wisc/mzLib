# Top-Down Search Engine: Algorithm Design

## Overview

A novel top-down proteomics search algorithm that applies the same identification-free, cross-run, binning-based philosophy developed for the DIA engine (see `Start-Up/DIA-Engine-Algorithm.md`) to intact-protein data. The central bet: **deconvolution is a consensus problem, not a per-scan problem. Bin and cluster first in raw m/z space, then deconvolve on the cross-run consensus.**

Top-down is hard because the MS1 signal of a single proteoform is convolved with (a) its isotopologue distribution and (b) its charge-state envelope. State-of-the-art tools (ProMex, pTop, TopPIC, FLASHDeconv) attempt to deconvolve scan-by-scan before feature detection and search. This algorithm inverts that order: bin the raw m/z signal, cluster features across runs in m/z space, then recover each envelope independently per charge and reconcile via parsimony + NNLS on a locked basis.

### Pipeline at a glance

```
Step 1  Per-file dual-resolution m/z indexing       (fine 0.01 + thick 0.1 Th)
Step 2  Identification-free RT alignment            (thick index)
Step 3  Cross-run feature detection + grouping      (fine index)
        3a  Intensity-ranked 2D flood-fill → FeatureBoxes
        3b  Binomial cross-run box matching → FeatureGroups
        3c  FLASHDeconv log-mass envelope grouping → ChargeStateEnvelopes
Step 4  Per-charge deconvolution + proteoform recovery
        4a  Per-charge averaged MS1 spectra         (SpectralAveraging)
        4b  Per-charge IsoDec deconvolution
        4c  Cross-charge parsimony
        4d  Raw-peak NNLS explained-fraction fit    (global basis)
        4e  Per-run NNLS quantification             (locked basis, per-run y)
        4f  Proteoform database matching
Step 5  MS2 aggregation → synthetic mzML + sidecar
Step 6  Two-stage FDR (PIP-ECHO-inspired)
```

Quantification is not a separate step — it falls out of 4e.

---

## Step 1: Per-File Dual-Resolution m/z Peak Indexing

### Goal
For each input LC-MS file, build **two** bin-indexed views of the raw MS1 centroids without any deconvolution or charge assumption. Both views feed later steps: the thick view drives alignment, the fine view drives feature detection, deconvolution, and quantification.

### Dual-Resolution Indexing
- **Fine index — 0.01 m/z bin width.** Built directly from the scan peaks in a single pass using `PeakIndexingEngine` at its default `BinsPerDalton = 100`. Fine enough that individual isotopologue peaks within a single charge state fall into separate bins, preserving the envelope structure that deconvolution will need in Step 4.
- **Thick index — 0.1 m/z bin width.** Derived from the fine index by combining every 10 consecutive bins, exposed as a `ThickIndexView` wrapper — no second pass over raw peaks. Coarse enough that a proteoform's lit bins are stable across runs under small mass-calibration drift.

### Why m/z, not mass
Binning on m/z avoids committing to a charge-state interpretation before we have cross-run evidence. **The DIA engine's 1-Th forbidden-region compression does not apply to top-down.** With realistic charge states z = 5 to ≈ 30+, isotope peaks from different `(proteoform, charge)` combinations can land essentially anywhere on the m/z axis; there are no forbidden regions to discard. So binning is not a compression step — it's a uniform discretization on which downstream feature grouping can operate. The cross-run win comes from Step 3's reproducibility matching, not from Step 1's sparsity.

### Output
- Per-file `PeakIndexingEngine` at 0.01 m/z resolution (fine index)
- Per-file `ThickIndexView` wrapper at 0.1 m/z resolution (thick index)

---

## Step 2: Identification-Free RT Alignment

### Goal
Align n top-down LC-MS files to a common retention time axis before any deconvolution or proteoform identification.

### Signal Selection
- For DDA top-down: use the full MS1 chromatogram
- For DIA top-down: pick the highest-TIC isolation window

### Alignment Procedure
Operates on the **thick (0.1 m/z) indices only**; the fine indices are untouched at this stage.

1. Select high-signal thick-index bins as RT alignment anchors across all n files.
2. Build an XIC per anchor in each file using `IndexingEngine.GetXic` + `PeakSpline` smoothing.
3. Detect XIC extrema and match across files (adapted from `FlashLFQ/IsoTracker/XICGroups`, with the identification dependency stripped).
4. Fit an RT transformation (piecewise linear or monotone spline, TBD).

### Output
A per-file RT warp mapping local RT → common RT axis. Applied to RT ranges at query time; the fine index itself is not re-timestamped.

---

## Step 3: Cross-Run Feature Detection & Grouping

### Conceptual Model
- **Rows** = individual top-down runs
- **Columns** = (RT × m/z) cells in the fine index, interpreted on the common aligned RT axis
- **Values** = binned m/z signal intensity

### 3a. Intensity-Ranked 2D Flood-Fill

**Goal.** Draw rough `(RT, m/z)` boxes around candidate features in each run (or in the cross-run consensus). No assumption about charge state, isotope spacing, or reproducibility at this stage — just a candidate list of bounding boxes.

**Procedure.**
1. **Find the global maximum.** Locate the most intense cell in the fine index.
2. **Trace outward in RT.** From the seed, walk the scan index in both directions using `IndexingEngine.GetXicByScanIndex` until `maxMissedScanAllowed` / `maxRTRange` / local-minimum (via `ExtractedIonChromatogram.FindPeakBoundaries`) ends the walk. This is the 1D flood-fill that `IndexingEngine.GetAllXics` already implements.
3. **Grow in m/z.** Walk ±1 fine bin (±0.01 Th) from the seed, gluing adjacent-bin XICs in whose RT range overlaps the seed with intensity monotonically decreasing across the frontier. Stop when the adjacent-bin XIC fails the overlap-and-decrease test (the m/z local minimum). This is the only genuinely new routine.
4. **Draw a bounding box.** Emit a `FeatureBox { MzRange, RtRange, SeedIntensity, TotalIntensity, SourceFile }`.
5. **Black out and repeat.** Mark consumed peaks in the `matchedPeaks` dictionary (already the mechanism in `GetAllXics`) so they cannot seed again. Return to step 1 on the remaining signal.
6. **Stop** when the next seed's intensity drops below a noise floor or a configured feature-count cap is reached.

**Variants to benchmark.** Per-file then intersect vs. flood-fill on the cross-run consensus (after summing/median of aligned fine indices).

### 3b. Cross-Run Binomial Feature Matching

**Goal.** For each donor `FeatureBox`, decide whether the same feature exists in each other acceptor run. Unit of matching is the whole box, not a single isotopologue.

**Key design point.** Top-down MS1 boxes contain tens to hundreds of peaks across their envelope and scan range — an order of magnitude more than a bottom-up peptide XIC's 3–5 isotopologues. That peak-count abundance makes **binomial scoring sharp** at the box level, which it isn't in bottom-up — and is why we cannot just copy FlashLFQ's feature-vector MBR scorer.

**Procedure.** For each donor `FeatureBox` and each acceptor run:
1. **Predicted-position probe.** Warp the donor box's RT range onto the acceptor's RT axis via Step 2's alignment (m/z is unchanged). Query the acceptor's fine index inside the predicted `(RT, m/z)` window; count matched peaks → `k_true`.
2. **Shifted null probes.** Query the same acceptor run at `(RT, m/z)` positions shifted outside the tolerances. Count matched peaks at each shifted probe → `k_null_1, …, k_null_j`.
3. **Binomial score.** Model matched-peak count as `Binomial(n = peaks_in_donor_box, p = mean(k_null_i) / n)`. Score = `1 − BinomialCDF(k_true − 1; n, p)`.
4. **Accept or reject.** Box is matched in the acceptor if the binomial p-value clears the threshold (calibrated in Step 6 Stage 1).

**Output.** A `FeatureGroup`: a `FeatureBox` plus the set of acceptor runs in which it was matched, each with its own `(k_true, score, intensity)`. A single proteoform produces multiple simultaneous `FeatureGroup`s — one per charge state in its envelope — all at the same RT. These co-eluting feature groups are the raw material for 3c.

### 3c. Charge-State Envelope Grouping (FLASHDeconv Log-Mass Trick)

**Goal.** Cluster co-eluting `FeatureGroup`s that belong to the same proteoform at different charge states. Output one candidate `(M, {z_i}, {FeatureGroups})` per RT cluster. See `FLASHDeconv-LogMass-Summary.md` for the derivation and source-code references.

**Key equation.**
```
log(m/z − m_proton) = log M − log z
```
`log z` depends only on charge, not on mass, so the set of log-m/z values produced by a single proteoform's charge envelope is a mass-independent template `{−log 1, −log 2, −log 3, …}` shifted by `log M`. The template shape is identical for every proteoform; only the shift changes. Inter-charge spacing `log(z+1) − log(z) ≈ 1/z` compresses at high charge — the template is non-uniform but fixed.

**Procedure.**
1. **RT-cluster the feature groups.** Group `FeatureGroup`s from 3b into RT clusters by RT proximity.
2. **Log-transform centroids.** `x_i = log(group[i].MzCenter − m_proton)` for each group in the cluster.
3. **Precompute template.** `template = { −log(z) : z ∈ 1..z_max }`, default `z_max = 60`.
4. **Binned sparse sweep.** Bin log-m/z values at width `2.5e-6 × tol_ppm`. For each candidate anchor (one per observed `x_i`, assumed to be at some charge `z_anchor ∈ 1..z_max`), compute implied `log M = x_i + log(z_anchor)`, shift the template to that anchor, count how many template positions hit non-empty bins → `k_hits`.
5. **Harmonic suppression.** For each candidate with `k_hits ≥ min_charges` (default 3), check harmonic aliases at orders `{2, 3, 5, 7, 11}`; suppress the weaker of any colliding pair.
6. **Greedy envelope assignment.** Pick the highest-scoring candidate, consume its matched `FeatureGroup`s, repeat on the leftovers.
7. **Emit envelopes.** `ChargeStateEnvelope { MonoisotopicMass, Charges[], MatchedFeatureGroups[], TemplateScore }` per cluster. Leftovers go forward as singletons for Step 4 to adjudicate.

**What 3c does not do.** No averagine scoring, no final mass refinement, no sequence matching. The emitted `M` is bin-accurate; Step 4d refines it against the full consensus envelope.

---

## Step 4: Per-Charge Deconvolution + Proteoform Recovery

Deconvolution finally happens — but never on a single raw scan, and never on the cross-run consensus as a single lump. Each `ChargeStateEnvelope` is deconvoluted **charge by charge** on per-charge averaged spectra, reconciled via parsimony, and validated + quantified via a raw-peak forward model.

### 4a. Per-Charge Averaged MS1 Spectra
For each `ChargeStateEnvelope` from Step 3c and each charge state `z_i` it contains:
1. Take the `FeatureBox` at `z_i` (its RT extent, m/z extent, member runs).
2. Select the top-N highest-TIC MS1 scans within that box's RT range across all contributing runs (~1 minute of RT is the typical window).
3. Average those scans via `SpectralAveraging` to produce one high-SNR averaged `MzSpectrum` **for that charge state only**.

Per-charge averaging preserves the isotopologue grid at each z. A single combined average over all charges smears it because different charges have different m/z spacing.

### 4b. Per-Charge Deconvolution (IsoDec Default)
Run `IsoDecAlgorithm` on each per-charge averaged spectrum. Output is a list of candidate monoisotopic masses for that charge state.

**IsoDec is the default** because its ML-trained envelope scorer is more robust to partial envelopes than `ClassicDeconvolutionAlgorithm`, and partial envelopes are exactly what we expect when averaging only the top-N scans inside a single charge's box. Classic remains available as a fallback.

### 4c. Cross-Charge Parsimony
Reconcile the per-charge mass lists from 4b:
- Group masses across charges by ppm tolerance (default 10 ppm).
- Rank candidates by **charge support count** — how many distinct charges independently recovered this mass.
- Emit `ParsimonyCandidateSet { MonoisotopicMass, SupportingCharges[], PerChargeIntensities[], SupportCount }` per envelope, sorted by support count.

A high support count is the cheap, model-free confidence signal: a real proteoform should be recoverable from most of its charge states independently; a noise-driven mass at one charge is unlikely to recur at others.

### 4d. Raw-Peak NNLS Explained-Fraction Fit
Validate the parsimony set against the averaged spectra via a non-negative least squares forward fit.

1. For each parsimony candidate `(M, z)`, generate the predicted isotopologue pattern at m/z via `Averagine` → `IsotopicDistribution`.
2. Assemble basis matrix `B` whose columns are the predicted patterns sampled onto the observed averaged spectra's m/z grid. `y` is the stacked per-charge averaged spectra from 4a.
3. Solve `min ||B·x − y||²  s.t. x ≥ 0`.
4. Compute `explained_fraction = ||B·x|| / ||y||`.
5. If `explained_fraction < threshold` (default 0.7), extract residual peaks `y − B·x`, return to 4c for a second parsimony pass, refit.

The fit is at the **raw peak level**, not the deconvoluted-mass level. Throwing away the per-isotopologue constraint at the scoring step wastes the whole reason per-charge averaging was done in 4a.

### 4e. Per-Run NNLS Quantification
Step 4d answered "are these species real?" but collapsed per-run information. To get quantification, **re-solve NNLS per run against the same locked basis `B`**:

1. For each run `r` in the envelope's supporting runs:
   - Average the top-N scans inside the envelope's RT window **from run r only**, per charge → per-run `y_r`.
   - Solve `min ||B·x_r − y_r||²  s.t. x_r ≥ 0` with the **same B** as 4d.
   - Compute `explained_fraction_r = ||B·x_r|| / ||y_r||` as per-run quant confidence.
2. Emit a `PerRunQuant { Species: (M, z_set), Run → (intensity, explained_fraction, residual) }` table.

**Properties.**
- **Shared basis across runs** — species identity is locked from the global fit.
- **Co-eluting proteoform deconvolution is correct by construction** — NNLS distributes intensity per-run according to observed peak ratios.
- **Sparse runs degrade gracefully** — the per-run explained fraction tells downstream stats how much to trust each value.
- **Replaces naive envelope summing** — single-column basis is the degenerate special case.

This is the engine's quantification output; there is no separate Step 6 for quant.

### 4f. Proteoform Database Matching
For each parsimony candidate that clears the explained-fraction threshold, match against a proteoform database (UniProt + GPTMD-style PTM expansion + terminal truncations) within mass tolerance. This is identification-level output, not a quant gate.

### Why this structure (vs. "deconvolute once on the consensus")
- **Per-charge averaging preserves the isotopologue grid** at each z.
- **Parsimony is cheap and model-free** — catches the low-hanging fruit before any forward-model cost.
- **Raw-peak NNLS** validates the parsimony set and surfaces missed components via residuals, with a basis kept small by the parsimony filter.
- **Restartable chain** — a poor explained fraction sends you back to 4c with residual peaks, no re-averaging or re-deconvoluting.
- **Quant and identification share infrastructure** — same basis, same solver, same averaging primitive.

---

## Step 5: MS2 Aggregation Per Envelope (Synthetic mzML + Sidecar Output)

The engine is identification-agnostic through Step 4f, but DDA runs will have triggered MS2 scans on many of these envelopes. Aggregating those MS2 scans across runs gives a high-SNR fragment spectrum that can be handed off to a downstream top-down search tool as a synthetic `.mzML`.

### 5a. MS2 Scan Selection
For each `ChargeStateEnvelope`, across every contributing run, collect MS2 scans satisfying:
- **RT filter:** scan's RT falls inside the envelope's union-of-per-charge-FeatureBox RT range.
- **Isolation filter:** scan's precursor isolation window intersects the m/z extent of at least one charge state in the envelope.

### 5b. Averaging Strategies (All Worth Benchmarking)
1. **Per charge state.** One averaged MS2 per `(envelope, charge)`. Max specificity, no cross-charge SNR gain.
2. **Per envelope.** Pool all MS2 across charges; one averaged spectrum per envelope. Max SNR, mixes precursor charges (OK for HCD, messier for ETD/UVPD).
3. **Per envelope region (low / mid / high charge bucket).** Intermediate: 2–3 groups along the envelope's charge range.

Implement all three; benchmark and pick a default.

### 5c. Synthetic mzML + Sidecar Output
- **Synthetic `.mzML`** — a new `MsDataFile` of averaged MS2 scans with precursor metadata, monotonic synthetic scan numbers, written via `MzmlMethods` / `AveragedSpectraWriter`. Consumable directly by downstream top-down identification tools. **No MGF writer needed** — mzML is a superset.
- **Envelope sidecar `.tsv`** — one row per `ChargeStateEnvelope` with `M`, charge set, RT range, per-run quant from 4e, explained fraction, and the synthetic scan numbers produced under each averaging strategy. The sidecar's scan numbers are the key that lets a downstream tool cross-reference fragment matches back to the engine's quant output.

---

## Step 6: FDR Control (Two-Stage, PIP-ECHO-Inspired)

### Stage 1: Per-Run Membership Scoring
For each run within a `ChargeStateEnvelope`, score the probability that the run's signal is a true member. Features:
- RT deviation from cluster centroid (post-alignment)
- **Envelope completeness in this run** — fraction of expected `(M, z)` peaks present
- **Intra-envelope isotopologue consistency** — does the observed pattern match `Averagine` at the assigned charge?
- Intensity consistency relative to other runs in the envelope
- Per-run NNLS explained fraction from 4e

Calibration: random-RT or random-m/z null model, as in PIP-ECHO.

### Stage 2: Cluster-Level Target-Decoy FDR
- Propagate target and decoy proteoform identities into each envelope
- Compete target vs. decoy at the `(deconvolved M, matched DB entry)` level
- **PTM-aware decoy generation** — reversed sequence with mirrored PTM positions rather than naive reversal (open question)

### Why Two Stages
- Stage 1: *"Is this run's signal really part of this envelope?"* (within-envelope error)
- Stage 2: *"Is this envelope really the proteoform we matched?"* (identification error)

Orthogonal error sources, same decomposition as the DIA engine.

---

## Relationship to the DIA Engine

| | DIA Engine | Top-Down Engine |
|---|---|---|
| Axis of binning | m/z (MS1 + MS2) | m/z (MS1 only, dual-resolution) |
| Bin width rationale | coarse (~1 Th), forbidden-region compression | fine (0.01 Th), no compression — m/z space is fully populated |
| What gets deferred | peptide identification | charge-state deconvolution *and* identification |
| Hardest convolution | chimeric MS2 | charge-state envelope + isotopologues |
| Deconvolution step | none (search directly) | Step 4, per-charge on averaged spectra |
| Quant source | MS2 fragment integration | Step 4e per-run NNLS on MS1 envelope |
| FDR framework | two-stage PIP-ECHO | two-stage PIP-ECHO (PTM-aware decoys) |

Both engines share the same three-part core: **identification-free alignment → cross-run feature grouping → match against consensus**. The top-down version adds deconvolution, but pushes it to after cross-run consensus is built, and per-charge rather than per-scan.

---

## Prior Art / Inspiration

- **DIA Engine Algorithm** (`Start-Up/DIA-Engine-Algorithm.md`) — direct parent; same alignment/grouping/FDR philosophy
- **PIP-ECHO** (Smith Lab) — two-stage FDR framework
- **FLASHDeconv** (Jeong / Kohlbacher) — log-mass template match for charge-envelope detection; Step 3c is adapted from it
- **IsoDec** (mzLib) — ML-trained per-charge deconvolution scorer used in Step 4b
- **ProMex / TopPIC / pTop** — existing top-down feature finders to benchmark against

---

## Open Questions

### Step 1 — Indexing
- Do 0.1 / 0.01 m/z fixed widths scale correctly to high m/z, or should they be ppm-based?
- `ThickIndexView` lazy vs. materialized — current plan is lazy.

### Step 2 — Alignment
- Anchor selection strategy: top-N by intensity, or top-N by cross-run reproducibility?

### Step 3 — Feature Detection & Grouping
- Flood-fill noise-floor stopping: fixed intensity, percentile, or derived per-file estimate?
- m/z-axis growth: concentric ring expansion vs. seeded bin walk — benchmark both.
- Per-file-then-intersect vs. on-the-consensus for the flood-fill — benchmark both.
- Binomial null probes: shift in RT only, m/z only, or both? Different error modes.
- Number of shifted null probes for robust `p_null` estimation.
- RT-cluster formation in 3c: pure RT proximity vs. co-intensity-correlation across runs.

### Step 4 — Deconvolution & Quant
- Top-N scan selection for 4a / 4e: fixed N, or dynamic by TIC fraction?
- Explained-fraction threshold (default 0.7) needs empirical tuning.
- NNLS solver: check `SharpLearning.Optimization` for a built-in, else ~80-line Lawson-Hanson.
- Residual-driven re-parsimony loop — how many iterations max?

### Step 5 — MS2 Aggregation
- Which of the three averaging strategies is the default? Benchmark-dependent.
- How often is the same proteoform fragmented in multiple runs in DDA top-down? Enough for cross-run MS2 averaging to pay off?

### Step 6 — FDR
- PTM-aware decoy generation: reversed sequence with mirrored PTM positions, shuffle with preserved PTM count, or something else?
- Integration with Consortium for Top-Down Proteomics proteoform-level FDR definitions.
- How does aggregation degrade when only a small subset of runs contains a given proteoform? Minimum `n`?
