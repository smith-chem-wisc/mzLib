# Top-Down Search Engine: Algorithm Design

## Overview

A novel top-down proteomics search algorithm that applies the same identification-free, cross-run, binning-based philosophy developed for the DIA engine (see `Start-Up/DIA-Engine-Algorithm.md`) to intact-protein data. The central bet: **deconvolution is a consensus problem, not a per-scan problem. Bin and cluster first in raw m/z space, then deconvolve on the cross-run consensus.**

Top-down is hard because the MS1 signal of a single proteoform is convolved with (a) its isotopologue distribution and (b) its charge-state envelope. State-of-the-art tools (ProMex, pTop, TopPIC, FLASHDeconv) attempt to deconvolve scan-by-scan before feature detection and search. This algorithm inverts that order: bin the raw m/z signal, cluster features across runs in m/z space, and only then attempt to resolve the charge-state envelope — now against a consensus signal built from n runs instead of a single noisy scan.

---

## Step 1: Retention Time Alignment (Identification-Free)

### Goal
Align n top-down LC-MS files to a common retention time axis before any deconvolution or proteoform identification.

### Signal Selection
- Compute integrated TIC across MS1 (top-down identifications are dominated by MS1 signal)
- For DDA top-down: use the full MS1 chromatogram
- For DIA top-down: pick the highest-TIC isolation window as in the DIA engine

### Binning Strategy (m/z space, no deconvolution, dual-resolution)

Binning is done directly on the **raw m/z axis**. Nothing is projected into mass space at this stage — the goal is precisely to avoid committing to a charge-state interpretation before we have cross-run evidence.

**Important departure from the DIA engine:** the DIA binning trick compresses signal by exploiting the quasi-integer mass defect of singly-charged tryptic fragments — 1 Th bins work because the space between integer-Th clusters is chemically forbidden and can be thrown away. **That compression does not apply to top-down.** With a realistic charge-state range of z = 5 to z ≈ 30+, isotope peaks from different (proteoform, charge) combinations can land essentially anywhere on the m/z axis. There are no forbidden regions to discard. So binning here is not a compression step — it is just a uniform discretization of the m/z axis to give feature grouping a tractable matrix to operate on.

**Dual-resolution indexing.** For each input file, build two indices in a single pass:

- **Thick index — 0.1 m/z bin width.** Used for the first-pass RT alignment (below). Coarse enough that a proteoform's lit bins are stable across runs under small mass-calibration drift, and sparse enough that anchor selection is cheap.
- **Fine index — 0.01 m/z bin width.** Reserved for Step 2 cross-run feature grouping and Step 3 consensus deconvolution. Fine enough that individual isotopologue peaks within a single charge state fall into separate bins, preserving the envelope structure that deconvolution will need.

The **fine index is built first** (one pass over the scan peaks at 0.01 m/z resolution), and the **thick index is derived from it by combining every 10 consecutive bins**. Peak-to-bin assignment happens exactly once; the thick index is a cheap downstream array transform, not a second indexing pass.

At the fine resolution, a single proteoform at charge z contributes one bin per isotopologue (isotopologues are 1/z Th apart; at z = 30 that is ~0.033 Th, well separated at 0.01 m/z bins), and the full envelope is a dense constellation of lit bins. At the thick resolution, adjacent isotopologues at high charge collapse into the same bin, which is exactly what makes the thick index a good *reproducibility* signal for alignment — it is robust to minor m/z drift.

**No charge-agnostic projection, no voting into mass bins.** The raw m/z binned signal is the only primitive. A proteoform with charge states z₁…zₖ will produce a structured pattern of lit bins — one dense isotopologue series per charge state, each series spaced 1/z_i Th apart. Recovering `M` and `{z_i}` from that pattern is the job of Step 3, not Step 1.

**No first-pass charge window restriction.** Because we are not deconvolving, there is nothing to restrict. Every m/z peak in the instrument's acquisition range contributes to a bin regardless of its (unknown) charge.

**Where the cross-run savings actually come from.** Without the DIA-style sparsity compression, the win has to come entirely from **cross-run reproducibility**: a real proteoform's lit-bin pattern repeats at the same RT and same m/z values across many runs, while noise does not. Step 2 is what does the work; Step 1 just provides a discretization on which Step 2 can operate.

### Alignment
- Use the **thick (0.1 m/z) index** for alignment; the fine index is untouched here
- Select high-signal thick-index bins as RT alignment anchors across all n files
- Fit an RT transformation (piecewise linear or spline, TBD) — same as DIA engine
- Output: a common RT axis; the fine index is then re-timestamped onto this axis for Step 2

Alignment happens on **m/z bins**, same axis as the instrument measured. Two runs whose proteoform charge-state distributions differ will still align, because the most intense thick bins — wherever they sit in the envelope — will tend to repeat across runs.

---

## Step 2: Cross-Run Feature Grouping (Identification-Free, m/z Space)

### Conceptual Model
- **Rows** = individual top-down runs
- **Columns** = (RT bin × m/z bin) cells on the common aligned axis
- **Values** = binned m/z signal intensity

### Rough Feature Detection — First Approach: Intensity-Ranked Flood-Fill

Before worrying about cross-run consistency, we need a mechanism for drawing boundaries around candidate features in the (RT, m/z) plane. First attempt, deliberately simple:

1. **Find the global maximum.** Locate the single most intense (RT bin × m/z bin) cell in the (aligned) signal — whether that signal is per-file or summed across runs is itself a knob to try (see below).
2. **Grow outward until local minima.** Starting from that seed cell, iterate outward in both the m/z direction and the RT direction. Keep extending a neighborhood as long as intensity is monotonically decreasing; stop at the first local minimum on each expansion frontier.
3. **Draw a bounding box.** The rectangle (or, optionally, the actual reached set) defined by the m/z and RT limits of the reached minima becomes the feature's extent.
4. **Black out and repeat.** Zero (or flag) the boxed-out region in the working signal so it cannot seed future features. Re-find the next most intense cell in what remains and repeat from step 1.
5. **Stop** when the next seed's intensity drops below a noise threshold or a configured feature count is reached.

This is deliberately a **rough** detector — it makes no assumption about charge-state envelope shape, isotope spacing, or reproducibility across runs. Its only job is to produce a candidate list of (RT, m/z) boxes that Step 2's grouping logic and Step 3's deconvolution can then refine.

### Where to Apply the Flood-Fill

Two variants to try, probably both:
- **Per-file then intersect.** Run flood-fill on each run's own fine index, then cross-reference the resulting box lists across runs to find reproducible boxes. Pro: no coordinate confusion. Con: noise in individual runs muddies the intensity ranking.
- **On the cross-run consensus.** Sum (or median) the aligned fine-index intensities across runs first, then flood-fill the consensus signal. Pro: noise averages down, faint features become visible. Con: a proteoform present in only one run gets washed out.

### Cross-Run Feature Matching (Box Propagation)

Once a `FeatureBox` exists in one run, we need to decide whether the same feature is present in each other replicate. This is conceptually match-between-runs / peptide-identity-propagation, but with two important differences from the bottom-up case:

1. **Top-down MS1 features contain many more peaks per box** than a bottom-up peptide isotope trace. A bottom-up feature is typically 3–5 isotopologues at a single charge state; a top-down proteoform box contains tens to hundreds of m/z peaks spread across isotopologues and scans. This peak-count abundance makes **binomial scoring statistically sharp** in a way it is not for bottom-up MBR.
2. **We are matching before identification**, so there is no peptide sequence, no predicted fragment pattern, and no prior RT model beyond the empirical Step 1 alignment.

### Matching Procedure

For each `FeatureBox` discovered in a "donor" run, and each "acceptor" run:

1. **Predicted-position probe.** Using the Step 1 RT warp, transform the donor box's RT range into the acceptor run's RT axis. The m/z range is unchanged (m/z is not warped). Query the acceptor's fine index inside this predicted `(RT, m/z)` window and count how many peaks match within `(ppm_tol, rt_tol)` — call this `k_true`.
2. **Shifted-position null probes.** Query the same acceptor run at **shifted** `(RT, m/z)` positions — RT shifts well outside the RT tolerance, m/z shifts well outside the ppm tolerance — and count matched peaks at each shifted position: `k_null_1, k_null_2, …`. These are the competitor null draws, analogous to PIP-ECHO's random-RT competition.
3. **Binomial score.** Model matched-peak count as a binomial with `n = (peaks in donor box)` and background rate `p_null = mean(k_null_i) / n` estimated from the shifted probes. The score of the true-position match is the binomial-tail probability of observing `≥ k_true` hits under `p_null`.
4. **Accept or reject.** A box is accepted as matched in the acceptor run if the binomial p-value clears the match threshold (threshold calibrated later in Step 4 Stage 1).

### Output
An **m/z feature group**: a `FeatureBox` plus the set of acceptor runs in which it was matched, each with its own `k_true`, null distribution, binomial score, and local-intensity readout. A box with matches in many runs is a reproducible feature group; a box with matches in few runs either is a low-prevalence proteoform or a noise artifact, and is demoted in Step 4.

Crucially, a single proteoform produces **multiple simultaneous feature groups** — one box per charge state in its envelope, all at the same RT. These co-eluting boxes are the raw material for Step 2c (charge-state envelope grouping), which is the bridge to Step 3.

---

## Step 2c: Charge-State Envelope Grouping (FLASHDeconv Log-Mass Trick)

### Goal
Take the cross-run-matched m/z feature groups from Step 2, cluster the ones that co-elute in RT, and identify **which subsets of them are the same proteoform observed at different charge states**. Output one candidate `(monoisotopic mass M, charge-state set {z_i})` per RT cluster — ready for Step 3 deconvolution and DB matching.

This is a pure grouping step. No per-scan deconvolution, no averagine scoring, no identification. It operates on the sparse list of feature-group m/z centroids within each RT co-elution cluster.

### The log-mass trick (adapted from FLASHDeconv)

See `FLASHDeconv-LogMass-Notes.md` and `FLASHDeconv-LogMass-Summary.md` for the full derivation and source-code references. The key equation is:

```
log(m/z − m_proton) = log(M / z) = log M − log z
```

Because `log z` is a function of charge only and does not depend on the proteoform mass, **the set of log-m/z values produced by a single proteoform's charge-state envelope is a mass-independent template `{−log 1, −log 2, −log 3, …}` shifted by a mass-dependent offset `log M`**. The template shape is identical for every proteoform; only the shift changes.

Note: inter-charge spacing is `log(z+1) − log(z) ≈ 1/z`, which **compresses at high charge** — the template is not uniformly spaced, it is a fixed non-uniform pattern.

### Procedure

1. **Take an RT cluster** — the set of co-eluting m/z feature groups output by Step 2's cross-run matching, grouped by RT proximity.
2. **Log-transform the feature m/z centroids**: `x_i = log(m/z_i − m_proton)` for each feature group `i` in the cluster.
3. **Build the shifted template.** Precompute `template = {−log 1, −log 2, …, −log z_max}` once. For each candidate monoisotopic mass `M`, the expected log-m/z positions for its charge-state envelope are `log M + template`.
4. **Sweep candidate masses via binned sparse template matching.** Bin the log-m/z axis at a ppm-tolerant width (FLASHDeconv uses `bin_mul_factor = 2.5 / (tol_ppm × 1e-6)`). For each candidate mass bin, count how many template positions land on a bin occupied by a feature group in this RT cluster. High hit counts identify candidate envelopes.
5. **Harmonic suppression.** For each candidate envelope, check against harmonic aliases at orders `{2, 3, 5, 7, 11}` (an envelope detected at `M` may be a z=1 ghost of an actual envelope at `2M`). Suppress the weaker of the two.
6. **Emit envelope assignments.** For each surviving candidate, output `(M, {z_i}, {matched FeatureGroups})` — the monoisotopic mass, the contributing charge states, and the cross-run feature groups that make up the envelope.
7. **Leftovers.** Feature groups that are not assigned to any envelope at this stage are carried forward as singletons. They might be real low-charge-count proteoforms, or they might be noise; Step 3 and Step 4 decide.

### What is NOT done in Step 2c

- **No averagine / isotope-pattern scoring.** Isotopologue shape is a Step 2 concern (inside the flood-fill) or a Step 3 concern (inside consensus deconvolution). Keeping Step 2c purely a charge-grouping step avoids committing to a particular isotope model at the grouping stage.
- **No final mass refinement.** The `M` emitted by the template match is bin-accurate, not centroid-accurate. Step 3 refines `M` against the full cross-run consensus envelope.
- **No sequence matching.** Step 2c is still identification-free.

### Why this works on cross-run feature lists (not just per-scan spectra)

FLASHDeconv runs the same primitive on a per-scan centroided peak list. Our input is different — a sparse list of cross-run `FeatureBox` m/z centers within an RT cluster — but the primitive is position-agnostic: the log-transform and template match care about m/z values, not about where they came from. Running it on cross-run feature groups rather than raw scans inherits the cross-run √n noise improvement from Step 2 before charge-state assignment is attempted, which is exactly the point of the whole engine's design.

---

### Output
A set of **m/z feature groups**, each with:
- A retention time position (common axis)
- An m/z bin position (centroid, refined sub-bin by weighted average of contributing peaks)
- A per-run intensity vector
- Adjacency metadata: which other feature groups co-elute (same RT cluster) — these are the candidate charge-state siblings

Nothing in this step commits to a mass, a charge state, or an identity. The output is a purely data-driven decomposition of the binned LC-MS signal into reproducible (RT, m/z) features.

---

## Step 3: Per-Charge Deconvolution + Parsimony + Envelope Explained-Fraction

Deconvolution finally happens here — but never on a single raw scan, and never on the cross-run consensus as a single lump. Instead, each charge-state envelope from Step 2c is deconvoluted **charge by charge**, and the per-charge mass lists are reconciled via parsimony and a raw-peak forward-model fit.

### 3a. Per-Charge Averaged Spectra
For each `ChargeStateEnvelope` from Step 2c, iterate its member charges. For each charge state `z_i`:
1. Take the `FeatureBox` at that charge state (RT extent, m/z extent, member runs).
2. Select the top-N highest-TIC MS1 scans within the box's RT range across all runs where the box is present (~1 minute of retention time is the typical window).
3. Average those scans via `SpectralAveraging` to produce a single high-SNR averaged MS1 spectrum **for that charge state only**.

Averaging per-charge rather than once over the whole envelope keeps isotope spacing crisp at each z — a single combined average would smear the isotopologue grid because different charges have different m/z spacing.

### 3b. Per-Charge Deconvolution (IsoDec)
Run **IsoDec** (`IsoDecAlgorithm`, default) on each per-charge averaged spectrum. Each run returns a candidate monoisotopic mass list for that charge state. These are *candidates*, not commitments — the parsimony step in 3c is what decides which ones are real.

`IsoDec` is the default because its ML-trained envelope scorer is more robust to partial envelopes than `ClassicDeconvolutionAlgorithm`, and partial envelopes are exactly what we expect when averaging only the top-N scans inside a single charge's box. `ClassicDeconvolutionAlgorithm` remains available as a fallback.

### 3c. Cross-Charge Parsimony
For each envelope, collect the per-charge deconvoluted mass lists from 3b and find recurring masses:
- Group masses across charge states by ppm tolerance (e.g., 10 ppm).
- Rank candidates by **charge support count** — how many distinct charge states in the envelope yielded this mass.
- A high charge-support count is the cheap, model-free confidence signal: a real proteoform should be recoverable from most of its charge states independently, whereas a noise-driven mass at one charge state is unlikely to recur at others.
- Emit a `ParsimonyCandidateSet { MonoisotopicMass, SupportingCharges[], PerChargeIntensities[] }` per envelope, sorted by support count.

### 3d. Explained-Fraction via Raw-Peak NNLS Fit
Given the parsimony candidate set, ask: **how much of the observed averaged-spectrum signal is explained by these candidates?**

1. For each candidate `(M, z)` pair, generate the predicted isotopologue pattern at that m/z via `Averagine` → `IsotopicDistribution` → m/z peaks at charge z.
2. Assemble a basis matrix `B` whose columns are the predicted isotopologue patterns (one column per `(M, z)` basis element), sampled onto the observed averaged spectrum's m/z grid.
3. Solve the non-negative least squares problem `min ||B·x − y||²   s.t. x ≥ 0`, where `y` is the stacked per-charge averaged spectra and `x` is the vector of per-species intensity coefficients.
4. Compute the **explained fraction** `||B·x|| / ||y||`. A high fraction (e.g., > 0.7) means the parsimony set captures most of the envelope; a low fraction flags an incomplete model and seeds a second pass (pull unexplained residual peaks, re-parsimony, re-fit).

The NNLS fit is done at the **raw peak level**, not the deconvoluted-mass level, because the per-isotopologue constraint is exactly the signal that distinguishes a real proteoform from a random mass hit — throwing it away at the scoring step wastes the whole reason we averaged per charge in 3a.

### 3e. Per-Run NNLS Quantification (the Quant Output)
Step 3d gave us a validated basis `B` of predicted isotopologue patterns for the `(M, z)` species that survived parsimony on the globally averaged spectrum. That fit answered "are these species real?" but collapsed per-run information. To get quantification, re-solve NNLS **per run** against the same basis:

1. For each run `r` in the envelope's supporting runs:
   - Select the top-N scans inside the envelope's RT window **from run r only**.
   - Average them per charge via `SpectralAveraging` → per-run, per-charge `MzSpectrum` vector `y_r`.
   - Solve `min ||B · x_r − y_r||²  s.t. x_r ≥ 0` with the **same basis B** as Step 3d.
   - `x_r` is the per-species intensity vector for run r.
   - Compute `explained_fraction_r = ||B · x_r|| / ||y_r||` as the per-run quant confidence.
2. Emit a `PerRunQuant { Species: (M, z_set), Run → (intensity, explained_fraction, residual) }` table.

Key properties:
- **Shared basis across runs** — species identity is locked from the global fit, so the per-run solves are pure intensity inference. No risk of the same physical species getting assigned to different basis elements in different runs.
- **Co-eluting proteoform deconvolution is correct by construction**. If two near-isobaric species are both in `B`, NNLS distributes intensity between them per-run according to the observed peak ratios — rather than double-counting them, which a naive envelope-sum would do.
- **Sparse runs degrade gracefully**. A run with few contributing scans still fits against the full global basis; its per-run explained fraction tells downstream stats how much to trust that quant value.
- **Replaces naive envelope summing** in the old Step 5. A one-column basis is the special case.

### 3f. Database / Library Matching
Parsimony-selected masses that clear the explained-fraction threshold are matched against a proteoform database (UniProt + GPTMD-style PTM expansion + terminal truncations) within mass tolerance.

Optional MS2 fragment match: for RT clusters where any contributing run triggered an MS2 scan, aggregate the MS2 fragments across runs and match against predicted b/y (HCD) or c/z (ETD/UVPD) ions from the candidate proteoform. MS2 consensus gets the same √n benefit as MS1 averaging.

### Why this structure (vs. "deconvolute once on the consensus")
- **Per-charge averaging preserves the isotopologue grid** at each z. A single combined average of all charges smears it.
- **Parsimony is cheap and model-free** — it needs only ppm-tolerance matching and a count. It catches the low-hanging fruit before any forward-model cost is paid.
- **Raw-peak NNLS is the forward model** that validates the parsimony set and surfaces missed components via residuals. Doing it after parsimony keeps the basis small (one column per parsimony candidate × supporting charge) instead of one per `(every possible mass × every possible charge)`.
- **The chain is restartable**: a poor explained fraction sends you back to 3c with the residual peaks, without re-averaging or re-deconvoluting.

---

## Step 3g: MS2 Aggregation Per Envelope (mzML + Sidecar Output)

The engine is identification-agnostic through Step 3f, but DDA runs will have triggered MS2 scans on many of these envelopes. Aggregating those MS2 scans across runs gives a high-SNR fragment spectrum that can be handed off to a downstream top-down search tool (e.g., MetaMorpheus, TopPIC, ProSight) as a synthetic `.mzML` file, without the engine itself having to know anything about fragment matching.

### 3g.1 MS2 Scan Selection
For each `ChargeStateEnvelope` from Step 2c, scan every contributing run for MS2 scans satisfying:
- **RT filter:** scan's RT falls inside the envelope's elution window (use the union of per-charge `FeatureBox` RT ranges).
- **Isolation filter:** scan's precursor isolation window intersects the m/z extent of at least one charge state in the envelope.

Both filters are cheap: RT is a `ScanInfo` lookup, isolation is `MsDataScan.IsolationRange` vs. per-charge m/z bounds.

### 3g.2 Averaging Strategies (All Worth Benchmarking)
Once the MS2 pool for an envelope is assembled, there are three ways to group and average:

1. **Per charge state.** One averaged MS2 per (envelope × charge). Maximum resolution on precursor-specific fragmentation, but no SNR benefit from pooling across charges. Most MS2 scans per group, fewest groups per envelope.
2. **Per envelope.** Pool all MS2 scans for the envelope regardless of precursor charge; one averaged MS2 per envelope. Maximum SNR, but mixes fragment spectra from different precursor charges — OK for HCD (charge-reduction products are predictable) but potentially messier for ETD/UVPD.
3. **Per envelope region (low / mid / high charge).** Intermediate: bucket charge states into 2–3 groups by their position along the envelope (low-charge → higher-m/z region vs. high-charge → lower-m/z region), average within each bucket. Compromise between the two extremes.

Initial plan: implement all three, benchmark on a dataset with known proteoforms, pick a default.

### 3g.3 mzML + Sidecar Output
Emit two files per project (or per run, configurable):

1. **Synthetic `.mzML` file** — a new `MsDataFile` whose scans are the averaged MS2 groups. Each synthetic MS2 scan carries the averaged peak list, the precursor m/z and charge (or charge set, for pooled strategies), the RT of the averaging window's midpoint, and the usual mzML selected-ion metadata. Fake scan numbers are assigned monotonically. This `.mzML` goes straight into any downstream top-down identification tool without format conversion. mzLib already has mzML writing via `Readers/MzML/MzmlMethods.cs` and `SpectralAveraging/AveragedSpectraWriter.cs`.
2. **Envelope sidecar file** (.tsv) — one row per `ChargeStateEnvelope` with: monoisotopic mass, charge-state set, elution RT range, per-run `PerRunQuant` from Step 3e, explained fraction, and the synthetic mzML scan numbers that came from averaging this envelope's MS2 pool under each strategy. The linkage between the synthetic mzML's scan numbers and the sidecar's envelope rows is what lets a downstream identification tool cross-reference a fragment-match result back to the engine's quant output.

### Reuse vs. Build
- **Reuse:** `MsDataScan.IsolationRange`, `MsDataScan.MsnOrder == 2` filter, `ScanInfo` for RT lookup, `SpectralAveraging` for the averaging primitive, **`MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra` (or similar)** + `AveragedSpectraWriter` for mzML output. No MGF writer needed.
- **Build:** MS2 pool construction per envelope, the three grouping strategies as pluggable functions, the sidecar .tsv writer, a thin wrapper that constructs `MsDataScan` objects from averaged peak lists with correct precursor metadata, and the logic that tracks the synthetic-scan → envelope mapping for the sidecar.

---

## Step 4: FDR Control (Two-Stage, PIP-ECHO-Inspired)

Same two-stage decomposition as the DIA engine.

### Stage 1: Per-Run Membership Scoring
For each run within an m/z feature group (or within an RT cluster of co-eluting m/z groups), score the probability that the run's signal is a true member.

Features for scoring:
- RT deviation from the cluster centroid (post-alignment)
- **Envelope completeness in this run** — of the co-eluting m/z bins that define the RT cluster's proteoform, how many are actually present with consistent intensity in this specific run? A real proteoform in a real run shows most of the envelope; a spurious hit shows one bin
- **Intra-bin isotopologue pattern consistency** — at higher MS1 resolution, the raw peaks inside a bin should follow the averagine envelope at the assigned charge; check this once a charge state is hypothesized in Step 3
- Intensity consistency relative to other runs in the cluster
- (Additional features TBD)

Calibration: random-RT or random-m/z null model, as in PIP-ECHO.

### Stage 2: Cluster-Level Target-Decoy FDR
- Propagate target and decoy proteoform identities into each RT cluster
- Compete target vs. decoy matches at the deconvolved (intact mass + fragments) level
- **PTM-aware decoy generation** — a reversed sequence with mirrored PTM positions is a more honest null than naive sequence reversal for proteoform search (open question below)

Combined FDR formula inherits PIP-ECHO logic, adapted.

### Why Two Stages
- Stage 1: *"Is this run's signal really part of this RT cluster?"* (within-cluster error)
- Stage 2: *"Is this cluster really the proteoform we deconvolved?"* (identification error, now including the deconvolution step)

These are orthogonal error sources — same decomposition as the DIA engine.

---

## Step 5: Quantification (MS1 Envelope Intensity)

Quantification returns to MS1 (unlike the DIA engine, which quants on MS2):
- Top-down MS1 is not chimeric the way DIA MS1 is — co-eluting proteoforms are generally separable in m/z and certainly after deconvolution
- Top-down MS2 is sparse (few precursors fragmented per cycle) and a poor quantification substrate
- The charge-state envelope is the native quantitative signal

### Quantification Strategy
For each proteoform passing FDR:
- Sum the raw m/z bin intensities belonging to its envelope (all charge states, all isotopologues) across the elution profile
- Weight per-run contributions by the Stage 1 membership score — runs with weak or incomplete envelopes contribute less
- Output: a proteoform-level intensity matrix (runs × proteoforms), ready for downstream DEA

---

## Full Pipeline Summary

```
n top-down files (.raw / .mzML)
    ↓
[Step 1] RT Alignment (dual-resolution raw m/z binning, no deconvolution)
    - Build per file: fine (0.01 m/z) index from raw peaks; thick (0.1 m/z) derived by 10-bin coarsening
    - No charge assumption, no mass projection, no forbidden-region compression
    - Thick index → anchor selection → RT warp across all n files
    - Fine index reserved for Step 2 and Step 3
    ↓
[Step 2] Rough Feature Detection + Cross-Run Grouping (still in m/z space)
    - Intensity-ranked flood-fill: seed at global max, expand to local minima, box out, repeat
    - Try per-file-then-intersect and on-the-consensus variants
    - Propagate each box to other runs at predicted (RT,m/z), score peak-count
      match via binomial tail vs. shifted-position null probes
    - Boxes that match in many runs → m/z feature groups
    - Co-eluting feature groups at the same RT are candidate charge-state siblings
    ↓
[Step 3] Deconvolution + Proteoform Matching (on consensus)
    - Within each RT cluster, solve (M, {z_i}) from the co-eluting bin hits
    - Match deconvolved mass against proteoform DB
    - Optional MS2 consensus fragment matching
    ↓
[Step 4] Two-Stage FDR Control
    - Stage 1: per-run membership scoring (envelope completeness, isotope consistency)
    - Stage 2: cluster-level target-decoy FDR (PTM-aware decoys)
    ↓
[Step 5] Quantification
    - MS1 envelope summation across elution peak (all charges, all isotopologues)
    - Membership-score-weighted per-run intensities
    - Output: proteoform intensity matrix → DEA Tool
```

---

## Relationship to the DIA Engine

| | DIA Engine | Top-Down Engine |
|---|---|---|
| Axis of binning | m/z (MS1 + MS2) | m/z (MS1 only) |
| Bin width rationale | coarse (~1 Th), forbidden-region compression | fine (~ppm), no compression — m/z space is fully populated |
| What gets deferred | peptide identification | charge-state deconvolution *and* identification |
| Hardest convolution | chimeric MS2 | charge-state envelope + isotopologues |
| Feature group = | reproducible (RT, m/z) fragment pattern | reproducible (RT, m/z) bin — multiple per proteoform |
| Deconvolution step | none (search directly) | Step 3, on the cross-run consensus |
| Quant source | MS2 fragment integration | MS1 envelope summation |
| FDR framework | two-stage PIP-ECHO | two-stage PIP-ECHO (PTM-aware decoys) |

Both engines share the same three-part core: **identification-free alignment → cross-run feature grouping → match against consensus**. The top-down version adds one extra step — deconvolution — but pushes it to after cross-run consensus is built, instead of doing it per scan up front. This is the inversion relative to existing top-down tools.

---

## Prior Art / Inspiration

- **DIA Engine Algorithm** (`Start-Up/DIA-Engine-Algorithm.md`) — direct parent; same alignment/grouping/FDR philosophy
- **PIP-ECHO** (Smith Lab) — the two-stage FDR framework
- **FLASHDeconv** (Kyowon Jeong / Oliver Kohlbacher) — fast top-down deconvolution, useful as a Step 3 benchmark
- **ProMex / TopPIC / pTop** — existing top-down feature finders to benchmark against

---

## Open Questions

- Are 0.1 / 0.01 m/z the right widths for thick / fine, or do they need to scale with m/z (e.g. ppm-based) at high m/z?
- How are adjacent m/z bins merged when a real peak straddles a bin boundary? Sliding window or overlap?
- Should the thick index be materialized as its own coarsened array, or exposed lazily as a query-time view over the fine index? (Current plan: lazy view wrapper over the fine index, no separate materialization.)
- In Step 2, how are co-eluting m/z feature groups grouped into "RT clusters" before deconvolution? Pure RT proximity, or co-intensity-correlation across runs?
- For the flood-fill feature detector: what defines a "local minimum" in 2D — strict monotone decrease along the expansion frontier, or a smoothed/threshold-based criterion? How are noisy plateaus handled?
- Does the flood-fill run on per-file signal first (then cross-referenced) or on the cross-run consensus (then decomposed)? Benchmark both.
- What noise-floor stopping criterion ends the flood-fill loop? Fixed intensity, percentile, or a derived per-file estimate?
- For the binomial cross-run matching: what is the right `n` — raw peak count inside the donor box, or a downsampled "resolution-independent" count? And what is the right number of shifted null probes to estimate `p_null` robustly?
- Should the shifted probes shift in RT only, m/z only, or both? Shifting only in RT tests against same-mass-different-feature contamination; shifting only in m/z tests against same-RT-different-mass contamination; these are different error modes and could both be scored.
- In Step 3, the deconvolution on consensus is the novel step. What's the scoring function? Candidate-by-candidate likelihood against averagine, or joint optimization over the whole RT cluster?
- **PTM-aware decoy generation**: reverse sequence with mirrored PTM positions? Shuffle with preserved PTM count? Still open in the field.
- How does cross-run aggregation degrade when only a small subset of runs contains a given proteoform? What is the minimum n?
- Native / intact MS (higher masses, wider charge distributions) — does the m/z bin count need to scale?
- MS2 consensus only works if the same proteoform is fragmented in multiple runs. How often does that actually happen in DDA top-down? Does this favor a targeted/inclusion-list acquisition strategy?
- Integration with proteoform-level FDR as defined by the Consortium for Top-Down Proteomics — does the two-stage FDR map onto their definitions?
