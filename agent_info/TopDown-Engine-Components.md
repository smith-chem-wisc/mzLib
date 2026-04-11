# Top-Down Search Engine: Components and Code Provenance

Companion to `TopDown-Engine-Algorithm.md`. This document decomposes the algorithm into concrete software components and, for each one, identifies:

- **Reuse** — existing mzLib code (with file paths) that can be used as-is or lightly adapted
- **Adapt** — existing mzLib code that needs meaningful modification or a new subclass to fit the top-down cross-run, consensus-first design
- **Build** — new components with no existing mzLib analog

All mzLib paths are relative to `mzLib/mzLib/` (the solution directory).

---

## 1. File I/O and Scan Access

**What it does.** Read Thermo `.raw` and `.mzML` (plus Bruker `.d` eventually), expose MS1 scans with peak arrays and retention times.

| Need | Status | Where |
|---|---|---|
| Raw/mzML readers | **Reuse** | `Readers/` — `Mzml`, `ThermoRawFileReader`, the `MsDataFile` abstraction |
| Per-scan peak arrays & RT | **Reuse** | `MassSpectrometry/MsDataScan.cs`, `MassSpectrometry/MzSpectra/MzSpectrum.cs` |
| Scan metadata index | **Reuse** | `MassSpectrometry/PeakIndexing/ScanInfo.cs` |

Nothing new here. This is the layer the entire algorithm sits on top of.

---

## 2. Step 1 — Raw m/z Peak Indexing (Dual-Resolution, Per-File)

**What it does.** For each of the n input files, build **two** RT × m/z indices of centroided MS1 peaks at different bin widths, with no deconvolution or charge assumption. The two indices serve different downstream stages:

- **Thick index — 0.1 m/z bin width** — feeds RT alignment (Step 1c). Coarse enough that a proteoform's lit bins are stable across runs under mass-calibration drift; sparse enough that anchor selection is cheap.
- **Fine index — 0.01 m/z bin width** — feeds cross-run feature grouping (Step 2) and consensus deconvolution (Step 3). Fine enough that individual isotopologue peaks within a charge state fall into separate bins, preserving the envelope structure needed for deconvolution.

**The fine index is built first, directly from the scan peaks using the existing `PeakIndexingEngine` at its default `BinsPerDalton = 100` (0.01 m/z bins). The thick index is then derived from the fine index by combining every 10 consecutive bins into one 0.1 m/z bin.** This avoids a second pass over the raw peaks: peak-to-bin assignment is done exactly once. The thick-derivation step is just an array coarsening over the bin dimension — each thick bin's peak list is the concatenation of its 10 source fine bins (or, if memory is a concern, an intensity-summed summary keyed by scan index).

| Need | Status | Where / Notes |
|---|---|---|
| Generic indexing base class (bin-index → list of peaks, with binary search on scan index) | **Reuse** | `MassSpectrometry/PeakIndexing/IndexingEngine.cs` — `abstract IndexingEngine<T>` provides the bin-array storage, `GetIndexedPeak`, `GetXic`, `GetXicByScanIndex`, and the binary search primitives we need |
| Per-file fine-resolution engine | **Reuse as-is** | `FlashLFQ/PeakIndexingEngine/PeakIndexingEngine.cs` — concrete `PeakIndexingEngine : IndexingEngine<IndexedMassSpectralPeak>`. Default `BinsPerDalton = 100` is already 0.01 m/z. Call `PeakIndexingEngine.InitializeIndexingEngine(scans)` unchanged. |
| Thick-index derivation (10× bin coarsening) | **Build** | Small new helper: take a built `PeakIndexingEngine`'s `IndexedPeaks` array and produce a coarsened view where bin `k` of the thick index aggregates fine bins `[10k, 10k+10)`. Preserves RT/scan indexing because the underlying `IIndexedPeak` objects still carry their scan metadata. |
| Indexed peak record | **Reuse** | `MassSpectrometry/PeakIndexing/IndexedMassSpectralPeak.cs` — already has `mz`, `intensity`, `retentionTime`, `zeroBasedScanIndex` |
| `IIndexedPeak` interface | **Reuse** | `MassSpectrometry/PeakIndexing/Interfaces/IIndexedPeak.cs` |
| Mass-based variant `MassIndexingEngine` | **Do not use** | `MassSpectrometry/PeakIndexing/MassIndexingEngine.cs` — indexes *deconvolved* envelopes by mass bin. Exactly the thing we are trying to avoid in the first pass. |

**Note on `BinsPerDalton`.** No modification to `PeakIndexingEngine` or to `IndexingEngine<T>` is needed. The default `BinsPerDalton = 100` gives us the fine (0.01 m/z) index for free, and the thick index is a pure downstream transform on the bin array.

Deliverable: a small `ThickIndexView` (or similar) that wraps a built `PeakIndexingEngine` and exposes the same `GetIndexedPeak` / `GetXic` / `GetXicByScanIndex` surface but against 0.1 m/z bins. Implementation is a 10-bin stride over the fine `IndexedPeaks[]` array; binary searches on `zeroBasedScanIndex` within each coarsened bin reuse the existing `BinarySearchForIndexedPeak` primitive unchanged. The top-down engine then owns one `PeakIndexingEngine` plus a derived `ThickIndexView` per input file.

---

## 3. Step 1 — Cross-Run Retention Time Alignment

**What it does.** Given the n per-file **thick (0.1 m/z)** indices, pick high-signal m/z bins as anchors and fit an RT warp from each file onto a common axis. The fine index is not touched during alignment — it is reserved for feature grouping and deconvolution in Steps 2–3 once the warp exists.

| Need | Status | Where / Notes |
|---|---|---|
| RT calibration data point + sort/score | **Adapt** | `FlashLFQ/MBR/RetentionTimeCalibDataPoint.cs` — currently keyed on matched identifications (`DonorFilePeak.Apex.IndexedPeak.RetentionTime`). We need an identification-free variant keyed on a matched *m/z bin anchor* instead of a peptide ID. |
| RT window / info | **Reuse** | `FlashLFQ/MBR/RtInfo.cs` |
| XIC-based cross-run alignment via shared extrema | **Adapt** | `FlashLFQ/IsoTracker/XICGroups.cs` + `Extremum.cs` + `PeakRegion.cs` already align XICs across runs via shared extrema, but require an `IdList` input. Strip the identification dependency and drive alignment off the top-N reproducible m/z bins instead. |
| XIC construction from the index | **Reuse** | `IndexingEngine<T>.GetXic(...)` and `GetXicByScanIndex(...)` in `IndexingEngine.cs` |
| Peak splines for smoothing XICs before extremum detection | **Reuse** | `MassSpectrometry/PeakIndexing/PeakSpline/` — `xICCubicSpline.cs`, `XicLinearSpline.cs`, `Bspline.cs` |
| Existing peak predictors from `Chromatography/` | **Do not use** | All of `Chromatography/RetentionTimePrediction/` (SSRCalc, Chronologer, CZE) does sequence-based *prediction* — we need cross-run empirical *alignment*. Wrong tool. |

Deliverable: a new `IdentificationFreeRtAligner` that takes N `TopDownRawMzIndexingEngine` instances, picks anchor bins, constructs XICs for each anchor in each file, detects extrema, and fits an RT warp (piecewise linear or monotone spline, TBD).

---

## 4. Step 2 — Rough Feature Detection + Cross-Run Grouping

**What it does.** Draw rough (RT, m/z) boxes around candidate features via an intensity-ranked flood-fill, then identify boxes that repeat across runs as m/z feature groups. Co-eluting feature groups at the same RT are the candidate charge-state siblings for Step 3.

### First-approach algorithm (see algorithm doc Step 2 for the full description)
1. Find global max cell in the signal
2. Expand outward in m/z and RT until local minima are reached
3. Draw bounding box, black out the region
4. Repeat until a noise-floor stopping criterion is hit

Two variants to benchmark: **per-file flood-fill then intersect** vs. **flood-fill on the cross-run consensus**.

### Major reuse: `IndexingEngine<T>.GetAllXics` is already a 1D flood-fill

**This is the biggest single reuse in the whole project.** `IndexingEngine<T>.GetAllXics` at `mzLib/MassSpectrometry/PeakIndexing/IndexingEngine.cs:198` already implements the intensity-ranked flood-fill skeleton — it is missing only the m/z-axis expansion:

| Flood-fill step | Existing mzLib code | Notes |
|---|---|---|
| Find global max, iterate intensity-ranked | `IndexingEngine.cs:202` — `IndexedPeaks... .OrderByDescending(p => p.Intensity).ToList()` | Exact pattern we want |
| "Already claimed" masking | `IndexingEngine.cs:201, 205, 213–216` — `matchedPeaks` dictionary prevents re-seeding | This is the "box out" mechanism |
| Trace outward in RT from a seed | `IndexingEngine.cs:208` — `GetXicByScanIndex(peak.M, ..., maxMissedScanAllowed, maxRTRange, ..., matchedPeaks)` | Walks scan index in both directions, stops on missed-scan count or RT range |
| Local-minimum peak boundary detection | `ExtractedIonChromatogram.FindPeakBoundaries` at `ExtractedIonChromatogram.cs:61` | Walks both directions from apex, tracks the running valley, cuts when current intensity exceeds valley by a discrimination factor. **Exactly the "iterate out to local minima" primitive.** |
| Cut peak to boundary region | `ExtractedIonChromatogram.CutPeak` at `ExtractedIonChromatogram.cs:136` → `RemovePeaks` at `:115` | Trims the reached set down to the bounded region |
| Noise-floor stopping | `numPeakThreshold` parameter in `GetAllXics` | Rejects XICs with too few peaks; can be extended with an intensity floor |

The 1D (RT-only) flood-fill is **already shipping**. The missing piece is **m/z-axis growth** — today `GetXicByScanIndex` stays inside `GetBinsInRange(mz, ppmTolerance)` (see `IndexingEngine.cs:229`), which is a fixed tolerance window around the seed m/z, so the RT walk never expands into adjacent m/z bins when signal extends there. Two ways to add m/z growth:

- **(a) Concentric ring expansion.** After each successful `GetAllXics` trace, re-run `GetBinsInRange` with an expanded m/z tolerance and see whether intensity is still decreasing along the m/z frontier. Stop when a local minimum is hit in m/z. Natural extension of the existing API.
- **(b) Seeded m/z walk.** Keep `GetAllXics`'s per-bin trace as the atomic operation, then post-process: for each XIC at seed m/z, walk ±1 fine bin repeatedly and check whether the adjacent-bin XIC's RT range overlaps the seed XIC's RT range with monotonically decreasing intensity. Glue adjacent-bin XICs together into a 2D box until the m/z frontier hits a local minimum.

Both approaches reuse the entire existing pipeline — no new peak-boundary logic, no new masking scheme, no new intensity ordering.

### Component table

| Need | Status | Where / Notes |
|---|---|---|
| Intensity-ranked seed iteration over bins | **Reuse** | `IndexingEngine<T>.GetAllXics` at `IndexingEngine.cs:198` — uses `OrderByDescending(p => p.Intensity)` over `IndexedPeaks` |
| "Claimed" masking to prevent re-seeding | **Reuse** | `matchedPeaks` dict in `GetAllXics` (`IndexingEngine.cs:201, 215`) |
| RT-axis outward trace from a seed | **Reuse** | `IndexingEngine<T>.GetXicByScanIndex` — stops on `maxMissedScanAllowed` + `maxRTRange` |
| Local-minimum boundary detection (RT axis) | **Reuse** | `ExtractedIonChromatogram.FindPeakBoundaries` — discrimination-factor valley detection, bidirectional from apex. This *is* your "iterate out to local minima" routine, already written. |
| Cut to boundary region | **Reuse** | `ExtractedIonChromatogram.CutPeak` / `RemovePeaks` at `ExtractedIonChromatogram.cs:115, 136` |
| Noise-floor stop on seed iteration | **Reuse (extend)** | `numPeakThreshold` in `GetAllXics`; add an intensity-floor parameter as a second stop condition |
| **m/z-axis expansion (the new bit)** | **Build** | Either concentric ring expansion via repeated `GetBinsInRange` calls with growing tolerance, or seeded-bin-walk with adjacent-bin XIC gluing. This is the only genuinely new routine in the 2D flood-fill. |
| Box / region data structure | **Build (small)** | `FeatureBox { MzRange, RtRange, SeedIntensity, TotalIntensity, SourceFile? }` — trivial. |
| Cross-run consensus signal (sum / median of aligned fine indices) | **Build** | Needs the Step 1c RT warp. On-the-fly summation across n files' fine indices at query time, or materialize a synthetic consensus `PeakIndexingEngine` whose `IndexedPeaks[]` holds aggregated pseudo-peaks. |
| Box-to-box cross-run matching (reproducibility) | **Build** | Per-file variant: overlap-based matching between each file's box list. Consensus variant: boxes already come from the aggregate, matching is trivial. |
| Co-elution grouping of adjacent m/z feature boxes into RT clusters | **Build** | Groups boxes whose RT ranges overlap into candidate charge-state siblings for Step 3. |
| Per-box intensity vector across runs | **Reuse** | `IndexingEngine<T>.GetXic` / `GetXicByScanIndex` on each per-file fine index with the box's (m/z, RT) extent. |
| Peak splines (optional smoothing before m/z-frontier decision) | **Reuse** | `MassSpectrometry/PeakIndexing/PeakSpline/` — can smooth an XIC along RT before the discrimination factor is applied |

### Implementation sketch for the first prototype
1. Call `fineIndex.GetAllXics(tolerance, maxMissedScansAllowed, maxRTRange, numPeakThreshold, cutPeakDiscriminationFactor: 0.6)` — this gives you 1D (RT-only) features for free.
2. For each returned `ExtractedIonChromatogram`, attempt m/z growth by querying adjacent fine bins (±1 bin = ±0.01 m/z). If the adjacent-bin XIC overlaps in RT with an intensity lower than the seed XIC and higher than noise, merge it into a `FeatureBox`. Repeat outward until the adjacent-bin XIC fails the overlap-and-monotone-decrease test (the m/z local minimum).
3. Emit a `FeatureBox` per merged cluster.

The 2D flood-fill reduces to "call `GetAllXics`, then extend each XIC outward in m/z using the same discrimination-factor logic that already exists along the RT axis." This is a much smaller prototype than originally scoped — probably on the order of 100 lines of new glue code rather than a from-scratch clustering implementation. Prototype this before anything else downstream.

### 4b. Cross-Run Feature Matching (Box Propagation)

**What it does.** Given a `FeatureBox` found in one donor run, decide whether the same feature exists in each other acceptor run by counting peak matches at the predicted `(RT, m/z)` location and comparing against shifted null probes via a binomial tail test. This is the concrete mechanism for turning a list of per-file boxes into a cross-run feature group. It is conceptually match-between-runs / PIP, but with a top-down-specific scoring primitive.

**Key design point.** Top-down MS1 boxes contain tens to hundreds of peaks across their envelope and scan range — an order of magnitude more than a bottom-up peptide XIC's 3–5 isotopologues. That peak-count abundance makes **binomial scoring sharp** at the box level, which it wouldn't be in bottom-up. This is why we can't just copy FlashLFQ's MBR scorer wholesale, and why the scoring primitive is the new piece worth writing.

| Need | Status | Where / Notes |
|---|---|---|
| RT warp applied to a donor box's RT range in the acceptor's coordinate system | **Adapt** | The per-file RT warp from Step 1c; apply it to `FeatureBox.RtRange`. No new machinery beyond Step 1c's `IdentificationFreeRtAligner`. |
| Peak lookup inside a predicted `(RT, m/z)` window on the acceptor's fine index | **Reuse** | `IndexingEngine<T>.GetXic(mz, retentionTime, ppmTolerance, missedScansAllowed, maxPeakHalfWidth)` already returns all peaks within a predicted `(RT, m/z)` window. Or use `GetBinsInRange` + scan-index filter for a tighter box query. |
| Peak-count match at a given position (the `k` in the binomial) | **Build (small)** | Thin wrapper: call `GetXic` (or a box-query variant), count how many peaks are within `(ppm_tol, rt_tol)` of the donor box's peak positions. Returns an integer. |
| Shifted-position null probes (random-RT / random-m/z competition) | **Build** | Run the same peak-count routine at RT shifts outside `rt_tol` and m/z shifts outside `ppm_tol`. Directly inspired by PIP-ECHO random-RT competition but applied to whole boxes instead of peptide IDs. |
| Binomial tail scoring | **Build** | `1 - BinomialCDF(k_true - 1; n, p_null)` where `p_null = mean(k_null_i) / n`. MathNet.Numerics (already a mzLib dependency via `MathNet.Numerics` in MassSpectrometry.csproj) provides `Binomial.CDF`. |
| Match-acceptance threshold / calibration | **Adapt** | `FlashLFQ/MBR/MbrScorer.cs` + `MbrScorerFactory.cs` — the scaffolding for "here is a score, here is a probability" exists, but the feature set for top-down boxes is different. Reuse the class structure; replace the feature calculation. |
| PEP calibration of the binomial p-values | **Reuse** | `FlashLFQ/MBR/PEP/` — the PIP-ECHO PEP machinery is already lineage-compatible with the null-model approach used here. |
| `RtInfo` (predicted RT + width) | **Reuse** | `FlashLFQ/MBR/RtInfo.cs` — direct reuse for passing predicted RT intervals around. |
| `RetentionTimeCalibDataPoint` | **Adapt** | `FlashLFQ/MBR/RetentionTimeCalibDataPoint.cs` — currently keyed on matched identifications; adapt the data-point idea to be keyed on matched `FeatureBox` anchors instead (same shape as Step 1c, reuse that adaptation). |
| Existing FlashLFQ MBR scorer feature set (6 features for peptide MBR) | **Do not reuse the feature set** | `MbrScorer` uses peptide-specific features (ppm mass error, isotope envelope correlation, intensity distribution, etc.). For top-down boxes the analogous features are (a) the binomial score, (b) retention time deviation from predicted, (c) envelope completeness once we know the charge states, (d) per-run intensity consistency. Different feature set, same scoring framework. |

### Implementation sketch for cross-run matching
1. Given donor `FeatureBox` with `n_peaks` peaks in the donor run and RT/m/z extent `(rt_range, mz_range)`:
2. For each acceptor run:
   - Warp `rt_range` onto the acceptor's RT axis via the Step 1c calibration
   - Call the acceptor's fine index for the peak count inside `(warped_rt_range, mz_range)` → `k_true`
   - For `j` in shifted null probe positions (shifts well outside the tolerances): count peaks at each `(rt_range + Δrt_j, mz_range)` and `(rt_range, mz_range + Δmz_j)` → `k_null_j`
   - `p_null = mean(k_null_j) / n_peaks`
   - `score = BinomialTail(k_true; n_peaks, p_null)`
   - Accept the match if `score < threshold`
3. Emit an m/z feature group containing the donor box plus all accepted acceptor-run matches and their per-run `(k_true, score, intensity)` triples.

### Key differences from FlashLFQ MBR
- **Unit of matching**: `FeatureBox` (many peaks) vs. peptide `ChromatographicPeak` (one isotopologue envelope)
- **Matching primitive**: binomial peak-count tail test vs. feature-set classifier
- **Identification dependency**: none vs. required (FlashLFQ MBR needs an identified donor peptide)
- **Null model source**: shifted-position probes on the same acceptor run vs. decoy peptide IDs
- **Statistical regime**: `n` in the tens-to-hundreds (binomial is sharp) vs. `n` = 3–5 (binomial is noisy; why FlashLFQ uses feature vectors instead)

### 4c. Step 2c — Charge-State Envelope Grouping (FLASHDeconv Log-Mass Trick)

**What it does.** Given an RT cluster of co-eluting `FeatureGroup`s (the output of 4b), group them into charge-state envelopes by performing a binned sparse template match in log-m/z space. The template `{−log 1, −log 2, −log 3, …}` is mass-independent, so it can be precomputed once and swept across candidate monoisotopic masses. Output: `(M, {z_i}, matched FeatureGroups)` per envelope. See `FLASHDeconv-LogMass-Summary.md` for the math and OpenMS source citations.

**Key design point.** This is a pure primitive over a sparse list of m/z centroids — no dense spectrum, no isotopologue fitting, no averagine scoring at this step. It runs on the tens-to-low-hundreds of `FeatureGroup` centroids within a single RT cluster, so the sweep is cheap. Averagine/isotope-envelope scoring is deliberately deferred to Step 3 consensus deconvolution; the flood-fill in Step 2a already captured the isotopologues implicitly inside each `FeatureBox`.

| Need | Status | Where / Notes |
|---|---|---|
| Proton mass constant for the `log(m/z − m_proton)` transform | **Reuse** | `Chemistry/Constants.cs` / `Chemistry/PeriodicTable` (proton mass is already used throughout `ClassExtensions.ToMass`/`ToMz`) |
| m/z ↔ mass conversion baseline | **Reuse** | `Chemistry/ClassExtensions.cs` — `ToMass`, `ToMz` — for sanity-checking the sweep outputs against existing machinery |
| Natural log + basic math | **Reuse** | `System.Math.Log` / `MathNet.Numerics` (already a mzLib dep via `MassSpectrometry.csproj`) |
| Ppm tolerance representation | **Reuse** | `MzLibUtil/Tolerance.cs` — bin width derivation `bin_mul_factor = 2.5 / (tol_ppm × 1e-6)` pulls tolerance from the standard `Tolerance` type |
| Any existing log-space indexing / template match | **None** | `PeakIndexingEngine` bins linearly in m/z, not logarithmically. No existing log-space machinery in mzLib. |
| Averagine isotope model (`Averagine`, `AverageResidue`) | **Do not reuse at this step** | `MassSpectrometry/Deconvolution/AverageResidue/Averagine.cs` exists and is the right tool, but envelope fit belongs in Step 3 consensus deconvolution, not at charge grouping. Using it here would re-introduce the per-scan, per-envelope scoring pattern that the whole design is trying to avoid. |
| Harmonic suppression at orders `{2, 3, 5, 7, 11}` | **Build** | Hardcoded integer multiplier/divisor check against the candidate charge. Reproduce the FLASHDeconv pattern (`harmonic_pattern_matrix_` in `SpectralDeconvolution.cpp`) as a small helper. |
| Envelope record type | **Build (small)** | `ChargeStateEnvelope { MonoisotopicMass, Charges[], MatchedFeatureGroups[], TemplateScore }` — trivial container. |
| Log-mz transform of `FeatureGroup` centroids | **Build (small)** | One-liner: `featureGroup.MzCenter → Math.Log(mz - Constants.ProtonMass)`. Wrap in a helper so the −proton subtraction is never forgotten (see summary doc — the proton subtraction is the whole point of the trick). |
| Template precomputation | **Build (small)** | Compute `{−log 1, −log 2, …, −log z_max}` once per engine instance. `z_max` configurable (default ~60 for top-down). |
| Binned sparse sweep across candidate monoisotopic masses | **Build** | Bin the log-mz values at width `~tol/2.5 ppm`, sweep the template offset across candidate `log M` values, count template hits. This is the only genuinely new indexing primitive in the whole project — but it is small (a dict from log-mz bin index → list of `FeatureGroup` indices + a loop over candidate shifts). |
| Candidate-mass enumeration strategy | **Build** | Either (a) enumerate one candidate per observed feature centroid, assuming it could be charge 1..z_max, or (b) sweep over a dense grid of log-M values. (a) is sparse and matches FLASHDeconv's implementation shape. |
| Leftover handling when a FeatureGroup doesn't fit any envelope | **Build** | Emit as a singleton envelope at implied z=1 or flag as ambiguous for downstream review. |

### Implementation sketch
1. Input: list of `FeatureGroup`s in one RT cluster, each with an m/z centroid and cross-run intensity.
2. Compute `logMz[i] = Math.Log(group[i].MzCenter - Constants.ProtonMass)` for each group.
3. Bin `logMz[i]` at width `2.5e-6 × tol_ppm` (≈ tol/2.5 ppm). Build a dict `bin → [groupIndices]`.
4. Precompute `template = { -Math.Log(z) for z in 1..z_max }`.
5. For each candidate anchor (one per observed logMz, assumed to correspond to some charge z_anchor): for `z_anchor` in `1..z_max`:
   - Implied `logM = logMz[i] + Math.Log(z_anchor)`
   - Shift template to that `logM`: `expectedBins = logM + template` (rebinned)
   - Count how many `expectedBins` hit non-empty bins in the dict → `k_hits`
   - Record `(M, z_set, k_hits)`
6. For each candidate with `k_hits ≥ min_charges` (e.g., 3): apply harmonic suppression — if another candidate at `(M / q)` or `(M × q)` for `q ∈ {2, 3, 5, 7, 11}` has a higher score, discard this one.
7. Greedy assignment: pick the highest-scoring surviving candidate, consume its matched `FeatureGroup`s, repeat on the leftovers until none remain.
8. Emit `ChargeStateEnvelope` records as the output of Step 2.

### Why this stays tiny
- No dense spectrum allocation — the "index" is a dict over tens to low hundreds of bins per RT cluster.
- No averagine math — that's Step 3's job, once we already have `(M, {z_i})` candidates.
- No FFT, no cross-correlation — just template offsets and dict lookups.
- The primitive is reusable: the same log-mz template match could later run on a dense MS1 spectrum (the original FLASHDeconv use case) if we ever want a per-scan fallback.

Prototype cost: probably 150–200 lines of new code, most of it the candidate sweep and harmonic suppression. Everything else (proton mass, log, tolerance, MathNet) is already in mzLib.

---

## 5. Step 3 — Per-Charge Deconvolution + Parsimony + NNLS Explained-Fraction

**What it does.** For each `ChargeStateEnvelope` from Step 2c, independently average and deconvolute each charge state's MS1 scans (IsoDec default), reconcile the per-charge mass lists via ppm-tolerance parsimony, and validate the resulting proteoform set with a raw-peak NNLS forward fit. Output: `ParsimonyCandidateSet` per envelope with an explained-fraction score.

### 5a. Per-Charge Averaged Spectra (Step 3a)

| Need | Status | Where / Notes |
|---|---|---|
| Scan selection inside a `FeatureBox` RT window | **Reuse** | `IndexingEngine<T>.GetXic` / `GetXicByScanIndex` + `ScanInfo` for RT→scan-index mapping. Filter to top-N by TIC inside the window. |
| TIC per scan | **Reuse** | `MsDataScan.TotalIonCurrent` |
| Spectrum averaging across scans | **Reuse** | `SpectralAveraging/` — existing library does exactly this. Call with the top-N selected scans for one charge's `FeatureBox` to get one averaged `MzSpectrum` per charge. |
| Per-charge averaging loop | **Build (small)** | Thin wrapper: `foreach charge in envelope → pick top-N scans in its box → SpectralAveraging.Average → store`. |

### 5b. Per-Charge Deconvolution (Step 3b)

| Need | Status | Where / Notes |
|---|---|---|
| **IsoDec deconvolution (default)** | **Reuse as-is** | `MassSpectrometry/Deconvolution/Algorithms/IsoDecAlgorithm.cs` + `IsoDecDeconvolutionParameters.cs`. Takes an `MzSpectrum`, returns candidate isotopic envelopes with `(mass, charge, score)`. Default here because its ML-trained scorer handles partial envelopes better than Classic. |
| Classic deconvolution (fallback) | **Reuse as-is** | `MassSpectrometry/Deconvolution/Algorithms/ClassicDeconvolutionAlgorithm.cs` + `ClassicDeconvolutionParameters.cs` |
| `Deconvoluter` entry point | **Reuse** | `MassSpectrometry/Deconvolution/Deconvoluter.cs` — standard `Deconvolute(MzSpectrum, DeconvolutionParameters)` call per averaged spectrum. No new subclass needed; we're deconvoluting one-charge-at-a-time on real averaged spectra, not on a synthetic consensus, so the existing API fits. |
| Isotopic envelope record | **Reuse** | `FlashLFQ/IsotopicEnvelope.cs` and the envelope type returned by the deconvolution algorithms |

**Note:** this is the biggest change from the previous plan. The earlier `ConsensusDeconvolutionAlgorithm : DeconvolutionAlgorithm` subclass is no longer needed — we hand real averaged `MzSpectrum` objects to the existing `IsoDecAlgorithm` unchanged. Lower implementation cost, higher reuse.

### 5c. Cross-Charge Parsimony (Step 3c)

| Need | Status | Where / Notes |
|---|---|---|
| Group masses across charge states within ppm tolerance | **Build (small)** | Take the per-charge deconvoluted mass lists, bucket into ppm-width groups, count distinct supporting charges. Dict + loop. |
| Ppm tolerance type | **Reuse** | `MzLibUtil/Tolerance.cs` |
| `ParsimonyCandidateSet` record | **Build (small)** | `{ MonoisotopicMass, SupportingCharges[], PerChargeIntensities[], SupportCount }` — trivial container |
| Ranking / tie-breaking | **Build (small)** | Sort by `SupportCount` desc, then by summed intensity |

### 5d. Raw-Peak NNLS Explained-Fraction Fit (Step 3d)

| Need | Status | Where / Notes |
|---|---|---|
| Predicted isotopologue pattern for a `(M, z)` candidate | **Reuse** | `Averagine` + `IsotopicDistribution` (`Chemistry/IsotopicDistribution.cs`) → m/z peak list at charge z. Already used by `ClassicDeconvolutionAlgorithm.FindIsotopicEnvelope`. |
| Basis matrix assembly | **Build** | For each parsimony candidate × supporting charge pair, generate the predicted peak pattern, sample onto the observed averaged spectrum's m/z grid. Result: `B ∈ ℝ^(n_observed_peaks × n_basis)`. |
| Observed-spectrum vector | **Reuse** | `MzSpectrum.YArray` from the per-charge averaged spectra, stacked into a single `y`. |
| **NNLS solver** | **Investigate then Build or Adapt** | MathNet.Numerics does **not** expose NNLS directly. Two options to investigate: (a) `SharpLearning.Optimization` is already a mzLib dep — check whether it has NNLS or only unconstrained/bounded optimizers. (b) Write a small active-set NNLS (Lawson-Hanson) — ~80 lines, standard reference algorithm. Prefer (a) if available; fall back to (b). |
| Explained fraction `||B·x|| / ||y||` | **Build (trivial)** | One dot product after NNLS solve |
| Residual peak extraction (for second-pass parsimony) | **Build (small)** | `y - B·x`, threshold, emit as a new candidate peak list |
| `EnvelopeFitResult` record | **Build (small)** | `{ ParsimonyCandidateSet, Coefficients, ExplainedFraction, ResidualPeaks }` |

### 5e. Per-Run NNLS Quantification (Step 3e / the quant output)

**What it does.** Lock the basis `B` from 5d, then re-solve NNLS independently per run against per-run averaged spectra. This is the actual quantification output of the engine — Step 5 in the old design is collapsed into this.

| Need | Status | Where / Notes |
|---|---|---|
| Per-run scan selection inside the envelope's RT window | **Reuse** | Same pattern as 5a, but filtered to one run's scans |
| Per-run, per-charge averaging | **Reuse** | `SpectralAveraging/` — same call as 5a, with a single run's scan list |
| NNLS solver | **Reuse (from 5d)** | Same solver as the global fit — the basis is locked, only the `y` vector changes per run |
| Basis reuse across runs | **Build (trivial)** | Pass the `B` matrix from 5d unchanged into each per-run solve |
| `PerRunQuant` record | **Build (small)** | `{ Species: (M, z_set), Dict<Run, (intensity, explained_fraction, residual)> }` |
| Per-run explained fraction as quant confidence | **Build (trivial)** | `||B·x_r|| / ||y_r||` per run |
| Results container | **Adapt** | `FlashLFQ/FlashLFQResults.cs` is peptide/protein oriented; build `TopDownResults` wrapping the `PerRunQuant` table |
| Output writers (.tsv, .parquet) | **Build** | Parquet is not in mzLib; needs `Parquet.Net` or similar in the new engine project (not in mzLib itself) |

**This replaces the old Section 8 (envelope-level quantification).** The naive "sum raw m/z peaks across the envelope" approach is a special case of this NNLS fit with a single-column basis, and gets co-eluting proteoforms wrong. Moving quant into Step 3 means the FDR step (old Section 7) now runs over `PerRunQuant` records directly rather than over a separately-computed intensity matrix.

### 5f. Database / Library Matching (Step 3f)
Same as the previous plan — moved from old Section 6 intact; see Section 6 below.

### 5g. MS2 Aggregation + Synthetic mzML Output (Step 3g)

**What it does.** For each `ChargeStateEnvelope`, pull the MS2 scans whose RT falls in the envelope's elution window *and* whose precursor isolation window intersects at least one of the envelope's charge-state m/z bounds. Average them via one of three grouping strategies, emit a synthetic mzML, and write an envelope sidecar .tsv.

| Need | Status | Where / Notes |
|---|---|---|
| MS2 scan enumeration with RT + MsnOrder filter | **Reuse** | `MsDataFile` scan iteration, `MsDataScan.MsnOrder == 2`, `MsDataScan.RetentionTime`, `ScanInfo` for index lookups |
| Precursor isolation window | **Reuse** | `MsDataScan.IsolationRange` (a `DoubleRange`); intersect with per-charge m/z bounds derived from the envelope's `(M, z)` pairs |
| Per-charge m/z window computation | **Reuse** | `Chemistry/ClassExtensions.ToMz` (derive from `M` and `z`) |
| MS2 averaging primitive | **Reuse** | `SpectralAveraging/` — same library used for Step 3a, but now on MS2 scan subsets |
| Three grouping strategies (per charge / per envelope / per envelope region) | **Build** | Pluggable functions that take the pooled MS2 list and return `List<Grouping>` where each `Grouping` is a set of scans to be averaged together |
| Envelope-region bucketing (strategy 3) | **Build (small)** | Split the envelope's charge-state set into 2–3 buckets ordered by charge (low-charge = higher-m/z, high-charge = lower-m/z) |
| Synthetic `MsDataScan` from averaged peaks | **Build (small)** | Construct an `MsDataScan` with `MsnOrder = 2`, averaged peak arrays, synthetic scan number, precursor m/z/charge metadata derived from the grouping |
| Synthetic `MsDataFile` assembly | **Reuse** | `MassSpectrometry/MsDataFile` — construct with the list of synthetic `MsDataScan`s |
| **mzML writer** | **Reuse** | `Readers/MzML/MzmlMethods.cs` — existing mzML writer. Also reference `SpectralAveraging/AveragedSpectraWriter.cs` for the prior art on writing averaged spectra into mzML (same pattern we need). **No MGF writer needed.** |
| Envelope sidecar .tsv writer | **Build** | CsvHelper (already a mzLib dep) — one row per envelope with `(M, z_set, RT range, PerRunQuant, explained_fraction, synthetic_mzml_scan_numbers_per_strategy)` |
| Synthetic-scan → envelope mapping | **Build (small)** | Dict maintained during the averaging loop; emitted into the sidecar |

Key point: the MGF writer concern evaporates. mzLib already writes mzML (via `MzmlMethods` and `AveragedSpectraWriter`), and mzML is a superset of MGF for downstream consumers — any top-down identification tool that accepts `.mgf` also accepts `.mzML`. Output format choice was a false constraint.

### Why no more `ConsensusDeconvolutionAlgorithm` subclass
The earlier design called for a new `ConsensusDeconvolutionAlgorithm : DeconvolutionAlgorithm` that took an RT cluster as input. The per-charge averaging approach obsoletes it: once you average per charge and get a real `MzSpectrum`, you can feed it to `IsoDecAlgorithm` unchanged. The only genuinely new glue is the per-charge averaging loop (5a), parsimony grouping (5c), and the NNLS forward fit (5d). All three are small.

---

## 6. Step 3 — Proteoform Database and Matching

**What it does.** Load proteins, expand into proteoforms (PTMs + terminal truncations), match deconvolved masses, and optionally match aggregated MS2 fragments.

| Need | Status | Where / Notes |
|---|---|---|
| FASTA / UniProt XML loading | **Reuse** | `UsefulProteomicsDatabases/ProteinDbLoader.cs`, `Loaders.cs`, `PtmListLoader.cs` |
| `Protein`, `Modification`, `ProteolysisProduct` types | **Reuse** | `Proteomics/` + `Omics/` |
| PTM data | **Reuse** | `PtmListLoader`, Unimod via `UnimodLoader.cs` |
| Proteoform enumeration (terminal truncations + combinatorial PTMs) | **Build** | mzLib's digestion logic is peptide-oriented. Proteoform-level expansion (intact protein + N/C truncations + combinatorial modifications) does not exist here — GPTMD-style expansion lives in MetaMorpheus, not mzLib. New. |
| Intact mass matching (search a mass against the proteoform list within tolerance) | **Build** | No direct analog. Conceptually simple — sorted list + binary search — but a dedicated `ProteoformMassIndex` should live in the new project. |
| Predicted fragment spectra for top-down | **Build** | mzLib has fragment-ion generation for peptides (`Proteomics/Fragmentation/` if present, otherwise Omics), but not top-down-style (long c/z series, internal fragments). Will need a top-down fragment predictor. |
| Cross-run MS2 aggregation | **Build** | No analog. New. |
| Tolerance types | **Reuse** | `MzLibUtil/Tolerance.cs` (PPM + absolute) |

---

## 7. Step 4 — Two-Stage FDR (PIP-ECHO-Inspired)

**What it does.** Stage 1 scores per-run membership of each run within an RT cluster; Stage 2 runs target–decoy competition at the cluster-identification level with PTM-aware decoys.

| Need | Status | Where / Notes |
|---|---|---|
| Scorer with multiple feature inputs | **Adapt** | `FlashLFQ/MBR/MbrScorer.cs` + `MbrScorerFactory.cs` — already provides a multi-feature MBR scorer in the Smith Lab's PIP-ECHO lineage. Reusable as the scaffolding; features to score need to be swapped for the top-down set (envelope completeness, intra-bin isotopologue fit, RT deviation, intensity consistency) |
| PEP / probability calibration | **Reuse** | `FlashLFQ/MBR/PEP/` — investigate contents; PIP-ECHO-style PEP machinery lives here |
| Random-RT null model (Stage 1 calibration) | **Build** | The PIP-ECHO random-RT competition has to be implemented against the new cross-run feature group structure. No exact analog. |
| Target–decoy FDR | **Reuse** | FlashLFQ already has target–decoy FDR plumbing in the MBR scorer; reusable for Stage 2 once targets/decoys are proteoforms rather than peptides |
| Sequence-reversed decoy generation | **Reuse** | `UsefulProteomicsDatabases/DecoyGeneration/` — reversed-sequence decoys for proteins exist |
| **PTM-aware decoy generation** | **Build** | Mirror PTMs across the reversed sequence (or shuffle while preserving PTM count) — does not exist in mzLib. Open scientific question as well as an engineering one. |

---

## 8. Step 5 — Envelope-Level Quantification

**What it does.** Sum the raw m/z peak intensities that make up a proteoform's envelope (all charge states, all isotopologues) across its elution profile, weighted by Stage 1 membership scores.

| Need | Status | Where / Notes |
|---|---|---|
| Peak integration across an elution window | **Adapt** | `FlashLFQ/ChromatographicPeak.cs` + the envelope-summation logic in `FlashLfqEngine.cs` does MS1 peak integration for peptides. Adapt to sum across the proteoform envelope's full charge-state set rather than a single peptide's isotopologues. |
| XIC apex and peak-boundary detection | **Reuse** | `MassSpectrometry/PeakIndexing/ExtractedIonChromatogram.cs` — `FindPeakBoundaries`, `RemovePeaks`; operates on `IIndexedPeak` so it plugs directly into the raw m/z index |
| Per-run weighting by membership score | **Build** | New. Trivial arithmetic once Stage 1 scores exist. |
| Results container | **Adapt** | `FlashLFQ/FlashLFQResults.cs` is peptide/protein-oriented. New `TopDownResults` with a proteoform-level intensity matrix. |
| Output writers (.tsv, .parquet) | **Build** | mzLib writes .tsv in various places via `CsvHelper`, but no `.parquet` support. Parquet is a stated requirement of the DEA Tool integration (see `Start-Up/DIA-Engine.md`). |

---

## 9. Supporting Utilities

| Need | Status | Where |
|---|---|---|
| Tolerance (ppm / absolute) | **Reuse** | `MzLibUtil/Tolerance.cs` |
| Double range / interval math | **Reuse** | `MzLibUtil/DoubleRange.cs` |
| m/z ↔ mass conversions | **Reuse** | `Chemistry/ClassExtensions.cs` — `ToMass`, `ToMz` |
| Element / isotope tables | **Reuse** | `Chemistry/PeriodicTable`, loaded via `UsefulProteomicsDatabases/Loaders.cs` |
| Chemical formula arithmetic | **Reuse** | `Chemistry/ChemicalFormula.cs` |
| Exceptions | **Reuse** | `MzLibUtil/MzLibException.cs` |

---

## Summary Matrix

| Algorithm step | Reuse | Adapt | Build |
|---|---|---|---|
| 1a. File I/O & scan access | `Readers/`, `MsDataFile/Scan`, `ScanInfo` | — | — |
| 1b. Raw m/z indexing (fine = default; thick = 10× coarsening) | `PeakIndexingEngine` (default 0.01 m/z), `IndexedMassSpectralPeak`, `BinarySearchForIndexedPeak` | — | `ThickIndexView` (10-bin coarsening wrapper) |
| 1c. Identification-free RT alignment (thick index) | `GetXic`, peak splines | `RetentionTimeCalibDataPoint`, `XICGroups` (drop ID dependency) | Anchor selection + warp fitter |
| 2a. Rough feature detection + cross-run grouping | **`GetAllXics` (1D flood-fill)**, **`FindPeakBoundaries`**, **`CutPeak`**, `matchedPeaks` masking, `PeakSpline/` smoothing, `GetXic` for box intensity reads | `GetAllXics` stop-condition extension | **m/z-axis growth, `FeatureBox`, consensus signal, box matching, RT co-elution grouping** |
| 2b. Cross-run `FeatureBox` matching | `GetXic`, `RtInfo`, `PEP/`, MathNet `Binomial.CDF` | `MbrScorer` scaffolding (swap feature set), `RetentionTimeCalibDataPoint` | Box-level peak-count primitive, shifted null probes, binomial tail score |
| 2c. Charge-state envelope grouping (log-mass template match) | `Chemistry/Constants` (proton mass), `System.Math.Log`, `Tolerance` | — | Log-mz transform, template precomputation, binned sparse sweep, harmonic suppression, `ChargeStateEnvelope` record |
| 3a. Per-charge averaged spectra | `SpectralAveraging/`, `MsDataScan.TotalIonCurrent`, `GetXic`, `ScanInfo` | — | Per-charge averaging loop |
| 3b. Per-charge deconvolution | **`IsoDecAlgorithm` (default)**, `ClassicDeconvolutionAlgorithm` (fallback), `Deconvoluter`, envelope types | — | — |
| 3c. Cross-charge parsimony | `Tolerance` | — | Ppm grouping, `ParsimonyCandidateSet` record |
| 3d. Raw-peak NNLS explained-fraction fit | `Averagine`, `IsotopicDistribution`, `MzSpectrum.YArray` | Investigate `SharpLearning.Optimization` for NNLS | Basis assembly, NNLS solver fallback (Lawson-Hanson), explained-fraction, residual extraction |
| 3e. Proteoform DB + matching | `ProteinDbLoader`, `Loaders`, `PtmListLoader`, `Tolerance` | — | Proteoform enumeration, intact-mass index, top-down fragment predictor, MS2 cross-run aggregation |
| 4. Two-stage FDR | `MbrScorer` scaffolding, `PEP/`, sequence-reversal decoys | Feature set swap, target–decoy to proteoforms | Random-RT null model, **PTM-aware decoy generation** |
| 5. Envelope-level quantification | `ExtractedIonChromatogram` boundaries | `ChromatographicPeak` integration, `FlashLFQResults` | Envelope sum across all z, membership weighting, parquet writer |

---

## Prototyping Order (Suggested)

1. **Raw m/z indexing** — smallest, highest reuse, unblocks everything else. Call `PeakIndexingEngine.InitializeIndexingEngine(scans)` at its default 0.01 m/z resolution, then write the `ThickIndexView` 10-bin coarsening wrapper on top. No modifications to mzLib itself.
2. **Rough feature detection + cross-run matching (Step 2a/2b)** — the scientific core and the biggest "build from scratch" item. Prototype before committing to the rest; if this step doesn't produce clean feature groups on real data, the whole design needs revisiting.
3. **Identification-free RT alignment** — can be bolted on once Step 2 reveals what reproducible features look like.
4. **Charge-state envelope grouping (Step 2c)** — cheap to implement once `FeatureGroup`s exist; validates that the cross-run feature list actually carries enough signal to recover `(M, {z_i})` without deconvolution.
5. **Consensus deconvolution** — only meaningful once Step 2c produces RT clusters with candidate charge assignments.
6. **Proteoform database + intact mass matching** — can proceed in parallel with 2–5 because it only needs a deconvolved mass list as its input contract.
7. **Two-stage FDR** — last, once there are real scores to calibrate against.
8. **Envelope-level quantification + output** — trivial arithmetic on top of (5) and (7); wire up at the end.

---

## Out-of-Scope Dependencies (Not in mzLib)

These live outside mzLib and need their own sourcing:

- **GPTMD-style proteoform expansion** (MetaMorpheus) — port or rewrite
- **Top-down fragment ion prediction** (c/z, internal fragments) — no mzLib module
- **Parquet writer** — add a new NuGet dep (`Parquet.Net` or `ParquetSharp`) to the new engine project; does not belong in mzLib itself
- **PTM-aware decoy strategy** — open research question, no implementation anywhere
