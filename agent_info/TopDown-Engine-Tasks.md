# Top-Down Engine: Task Decomposition & Governing Checklist

Work-plan companion to `TopDown-Engine-Algorithm.md` and `TopDown-Engine-Components.md`. Decomposes the engine into discrete tasks sized for independent (or subagent-dispatched) execution.

**Task format.** Each task is written so a cold-start contributor (human or subagent) can pick it up without re-deriving context. Every task specifies:
- **ID** (stable, used in dependency refs)
- **Inputs** (what must already exist)
- **Outputs** (what the task produces)
- **Files to touch** (pointers, not hard constraints)
- **Acceptance criteria** (how to tell it's done)
- **Deps** (other task IDs that must complete first)
- **Doc refs** (algorithm + components doc sections)

**Dispatchability.** Tasks marked `[agent-OK]` are self-contained and suitable for subagent dispatch — clear inputs, outputs, and acceptance criteria. Tasks marked `[needs-review]` require user judgment or scientific calls that a subagent shouldn't make alone.

**Project location.** The top-down engine lives in a new project outside mzLib; mzLib changes should be minimized. When a task requires touching mzLib (e.g., new public API), flag it explicitly — these trigger the MetaMorpheus integration CI and need extra care.

---

## Governing Checklist

Master milestones. Each corresponds to a group of tasks below. Tick off as they complete.

### Infrastructure / Preamble
- [ ] **M0.** New top-down engine project scaffolded and building against mzLib
- [ ] **M1.** Test harness project scaffolded with simulator placeholder

### Step 1 — Dual-Resolution Indexing
- [ ] **M2.** `ThickIndexView` wrapper over `PeakIndexingEngine` working and unit-tested

### Step 2 — RT Alignment
- [ ] **M3.** Identification-free RT aligner working on synthetic data
- [ ] **M4.** RT aligner validated on real data vs. known landmarks

### Step 3 — Feature Detection & Grouping
- [ ] **M5.** 2D flood-fill feature detection (RT + m/z growth) producing `FeatureBox`es
- [ ] **M6.** Binomial cross-run box matcher producing `FeatureGroup`s
- [ ] **M7.** FLASHDeconv log-mass envelope grouping producing `ChargeStateEnvelope`s

### Step 4 — Deconvolution + Quant
- [ ] **M8.** Per-charge averaging + IsoDec producing candidate mass lists
- [ ] **M9.** Parsimony reconciliation producing `ParsimonyCandidateSet`s
- [ ] **M10.** NNLS explained-fraction fit working on global averaged spectra
- [ ] **M11.** Per-run NNLS quantification producing `PerRunQuant` tables
- [ ] **M12.** Proteoform DB matching wired up

### Step 5 — MS2 Aggregation
- [ ] **M13.** MS2 pool + averaging strategies + synthetic mzML writer + sidecar

### Step 6 — FDR
- [ ] **M14.** Stage 1 per-run membership scorer
- [ ] **M15.** Stage 2 target-decoy + PTM-aware decoy generation

### Testing Suite (parallel stream)
- [ ] **M16.** Real-data ingestion (raw + multi-tool search results)
- [ ] **M17.** Consensus feature extraction from multi-tool intersection
- [ ] **M18.** Empirical characterization module
- [ ] **M19.** MVP simulator (parametric)
- [ ] **M20.** Engine test harness with recovery metrics

---

## Tasks — Infrastructure

### T0.1 — Scaffold the top-down engine project [agent-OK]
**Inputs:** mzLib repo at current HEAD.
**Outputs:** New folder under `mzLib/mzLib/` (e.g., `TopDownEngine/`) with a `.csproj` targeting `net8.0`, x64, referencing `MassSpectrometry`, `FlashLFQ`, `Chemistry`, `MzLibUtil`, `UsefulProteomicsDatabases`, `SpectralAveraging`. Empty `Engine.cs` entry point.
**Files to touch:** `mzLib/TopDownEngine/TopDownEngine.csproj`, `mzLib/mzLib.sln` (add project).
**Acceptance:** `cd mzLib && dotnet build --no-restore ./TopDownEngine/TopDownEngine.csproj` succeeds on Windows CI. Do NOT add to `mzLib.nuspec` yet — this is a private-to-the-engine project until the API stabilizes.
**Deps:** none.
**Doc refs:** Components §1–2.

### T0.2 — Scaffold the test harness project [agent-OK]
**Inputs:** T0.1 complete.
**Outputs:** `mzLib/TopDownEngine.Test/` targeting `net8.0-windows`, NUnit 4, referencing `TopDownEngine` + test-data copies from existing `Test/` conventions.
**Acceptance:** `dotnet test ./TopDownEngine.Test/` runs (zero tests pass) on Windows CI.
**Deps:** T0.1.
**Doc refs:** Testing §A.

---

## Tasks — Step 1: Dual-Resolution Indexing (→ M2)

### T1.1 — `ThickIndexView` wrapper [agent-OK]
**Inputs:** A built `PeakIndexingEngine` (fine index, 0.01 m/z bins).
**Outputs:** A new class in `TopDownEngine/Indexing/ThickIndexView.cs` that wraps a `PeakIndexingEngine` and exposes `GetIndexedPeak`, `GetXic`, `GetXicByScanIndex`, `GetBinsInRange` at an effective 0.1 m/z bin width by coarsening every 10 fine bins. Implementation: stride over the fine `IndexedPeaks[]` array; binary searches on `zeroBasedScanIndex` within each coarsened bin reuse the existing `BinarySearchForIndexedPeak` primitive.
**Acceptance:**
- Returns the same peak set as querying the fine index at the same `(m/z ± 0.05, RT window)` coordinates (modulo expected bin-boundary effects).
- Unit test with a hand-constructed peak list that crosses a 0.1 Th boundary.
- No modifications to `IndexingEngine<T>` or `PeakIndexingEngine`.
**Deps:** T0.1.
**Doc refs:** Algorithm Step 1; Components §2.

### T1.2 — Dual-index builder entry point [agent-OK]
**Inputs:** T1.1 complete.
**Outputs:** `TopDownEngine/Indexing/DualIndexBuilder.cs` with one method: `(PeakIndexingEngine fine, ThickIndexView thick) Build(MsDataFile file)`. Calls `PeakIndexingEngine.InitializeIndexingEngine(scans)` with `BinsPerDalton = 100` and wraps the result.
**Acceptance:** Unit test on a small `.mzML` fixture confirms both indices return expected peaks for known m/z/RT queries.
**Deps:** T1.1.
**Doc refs:** Algorithm Step 1; Components §2.

---

## Tasks — Step 2: RT Alignment (→ M3, M4)

### T2.1 — Anchor selection from thick index [agent-OK]
**Inputs:** N `ThickIndexView`s (one per file).
**Outputs:** `TopDownEngine/Alignment/AnchorSelector.cs` — picks the top-K m/z bins by cross-run reproducibility (bin present in ≥ `min_files_for_anchor` with intensity above threshold). Returns `AnchorBin[]`.
**Acceptance:** Unit test with 3 synthetic files where two share a set of anchor bins + noise; selector picks the shared set.
**Deps:** T1.1, T1.2.
**Doc refs:** Algorithm Step 2.

### T2.2 — XIC extraction for anchors [agent-OK]
**Inputs:** T2.1.
**Outputs:** Helper that calls `IndexingEngine.GetXic(...)` for each anchor in each file, smooths via an existing `PeakSpline` variant, and returns extrema lists.
**Acceptance:** Unit test on a synthetic file where a known Gaussian peak at a known RT produces one extremum at the expected RT.
**Deps:** T2.1.
**Doc refs:** Components §3.

### T2.3 — `IdentificationFreeRtAligner` [agent-OK]
**Inputs:** T2.1, T2.2.
**Outputs:** `TopDownEngine/Alignment/IdentificationFreeRtAligner.cs` — consumes anchor extrema across files, matches them, fits a per-file RT warp (piecewise linear initially; monotone spline as a flag). Adapted from `FlashLFQ/IsoTracker/XICGroups` with the ID dependency stripped. Returns `Dictionary<SpectraFileInfo, Func<double, double>>` (per-file RT warp).
**Acceptance:**
- Unit test with synthetic file pair offset by a known linear RT drift recovers the drift within 1 scan.
- Unit test with a non-linear drift produces a non-linear warp that reduces RMS RT error below a threshold.
**Deps:** T2.2.
**Doc refs:** Algorithm Step 2; Components §3.

### T2.4 — Validate on real data [needs-review]
**Inputs:** T2.3 complete, real raw files from the testing suite available.
**Outputs:** Notebook or test fixture that runs the aligner on the two real raw files and reports RMS RT error against manually-picked landmark peaks.
**Acceptance:** User review of the result plot.
**Deps:** T2.3, T16.1 (raw ingestion).
**Doc refs:** Testing §B.

---

## Tasks — Step 3: Feature Detection & Grouping (→ M5, M6, M7)

### T3.1 — `FeatureBox` record [agent-OK]
**Inputs:** none.
**Outputs:** `TopDownEngine/Features/FeatureBox.cs` — `{ MzRange, RtRange, SeedIntensity, TotalIntensity, SourceFile, PeakCount }`. Immutable.
**Acceptance:** Compiles; simple round-trip test.
**Deps:** T0.1.
**Doc refs:** Algorithm 3a; Components §4.

### T3.2 — 2D flood-fill via `GetAllXics` + m/z growth [needs-review]
**Inputs:** T3.1, a built fine `PeakIndexingEngine`.
**Outputs:** `TopDownEngine/Features/FloodFillDetector.cs`. Calls `fineIndex.GetAllXics(...)` to get 1D XICs, then for each XIC grows in the m/z direction by walking ±1 fine bin and gluing adjacent-bin XICs whose RT range overlaps monotonically-decreasingly. Emits `FeatureBox[]`.
**Notes for implementer:**
- The 1D flood-fill is already `IndexingEngine.GetAllXics` at `IndexingEngine.cs:198` — do NOT reimplement it. Only add the m/z growth wrapper.
- Two growth strategies to implement behind a flag: concentric ring (re-query `GetBinsInRange` with growing tolerance) and seeded bin walk (glue adjacent single-bin XICs). Default TBD after benchmarking.
- Noise-floor stop: pass through as a parameter, default = percentile of fine-index intensity distribution.
**Acceptance:**
- Unit test on synthetic data with three known features at known `(m/z, RT)` recovers all three boxes with bounds within a tolerance.
- Unit test with noise-only input produces zero boxes at reasonable parameters.
- Integration test flagged for review before merge — this is where scientific judgment matters most.
**Deps:** T3.1.
**Doc refs:** Algorithm 3a; Components §4.

### T3.3 — Binomial peak-count primitive [agent-OK]
**Inputs:** A `PeakIndexingEngine` and a `(RT window, m/z window)` query.
**Outputs:** Helper in `TopDownEngine/Features/BinomialScorer.cs` that returns `int CountPeaksInBox(PeakIndexingEngine index, DoubleRange rtRange, DoubleRange mzRange)`. Thin wrapper over `GetXic` / `GetBinsInRange`.
**Acceptance:** Unit test with hand-built peak arrays verifies inclusive/exclusive boundary handling and scan-index filtering.
**Deps:** T1.2.
**Doc refs:** Algorithm 3b; Components §4b.

### T3.4 — Shifted null probes + binomial tail scoring [agent-OK]
**Inputs:** T3.3 + `MathNet.Numerics.Distributions.Binomial`.
**Outputs:** `BinomialScorer.ScoreMatch(donorBox, acceptorIndex, rtWarp, nullShifts)` returning `{ k_true, p_null, pValue }`.
**Notes:**
- Implement RT-only shifts, m/z-only shifts, and combined as separate methods — see open question in algorithm doc.
- Use `MathNet.Numerics.Distributions.Binomial.CDF`; the p-value is `1 − CDF(k_true − 1; n, p_null)`.
**Acceptance:**
- Unit test: a clearly-matched box (high k_true, low k_null) yields a low p-value.
- Unit test: a random-noise box yields a p-value > 0.1.
**Deps:** T3.3.
**Doc refs:** Algorithm 3b; Components §4b.

### T3.5 — `FeatureGroup` record + cross-run matcher [agent-OK]
**Inputs:** T3.4, T2.3 (RT warp).
**Outputs:** `FeatureGroup { DonorBox, Matches: Dict<SpectraFileInfo, (k_true, pValue, intensity)>, TotalSupport }`. A `CrossRunFeatureMatcher` that iterates donor boxes and populates `FeatureGroup`s.
**Acceptance:**
- Integration test on synthetic 3-run data with a known matched feature: all three runs appear in the `FeatureGroup`.
- Integration test with a donor-only feature (noise in other runs): other runs not matched.
**Deps:** T3.4, T2.3.
**Doc refs:** Algorithm 3b; Components §4b.

### T3.6 — RT-cluster formation [agent-OK]
**Inputs:** T3.5.
**Outputs:** `TopDownEngine/Features/RtClusterer.cs` — groups `FeatureGroup`s by RT proximity (window configurable). Returns `List<RtCluster>`.
**Acceptance:** Unit test with three synthetic features at RTs `[10.0, 10.2, 15.0]` produces two clusters.
**Deps:** T3.5.
**Doc refs:** Algorithm 3c preamble.

### T3.7 — Log-mz transform + template precomputation [agent-OK]
**Inputs:** none.
**Outputs:** `TopDownEngine/Envelopes/LogMassTransform.cs` with:
- `static double LogMz(double mz)` → `Math.Log(mz - Constants.ProtonMass)` (use `Chemistry.Constants` for the proton mass).
- `static double[] BuildTemplate(int zMax)` → `{ -Math.Log(1), -Math.Log(2), …, -Math.Log(zMax) }`.
**Acceptance:** Unit test: known `(M, z)` round-trips through the log transform to a value within `1e-12` of `log(M) - log(z)`.
**Deps:** T0.1.
**Doc refs:** Algorithm 3c; Components §4c; `FLASHDeconv-LogMass-Summary.md`.

### T3.8 — Binned sparse template sweep [needs-review]
**Inputs:** T3.7.
**Outputs:** `TopDownEngine/Envelopes/EnvelopeDetector.cs` — takes an `RtCluster` and a ppm tolerance, performs the binned sparse sweep: bin log-mz values at width `2.5e-6 × tol_ppm`, enumerate candidate anchors, compute `k_hits`, return a ranked list of candidate envelopes. Excludes harmonic suppression (next task).
**Acceptance:**
- Unit test with a synthetic cluster of `(M, z ∈ {5,6,7,8,9})` centroids recovers the correct `M` as the top candidate.
- Unit test with a cluster of unrelated centroids produces no candidate above `min_charges` threshold.
**Deps:** T3.7.
**Doc refs:** Algorithm 3c step 4; Components §4c.

### T3.9 — Harmonic suppression [agent-OK]
**Inputs:** T3.8.
**Outputs:** Extension of `EnvelopeDetector` — for each candidate, check colliding candidates at `M × q` and `M / q` for `q ∈ {2, 3, 5, 7, 11}`; suppress the weaker.
**Acceptance:** Unit test: a synthetic cluster containing both a z=1 ghost and a true `2M` envelope retains only the `2M` candidate.
**Deps:** T3.8.
**Doc refs:** Algorithm 3c step 5; `FLASHDeconv-LogMass-Summary.md`.

### T3.10 — Greedy envelope assignment → `ChargeStateEnvelope` [agent-OK]
**Inputs:** T3.9.
**Outputs:** `ChargeStateEnvelope { MonoisotopicMass, Charges[], MatchedFeatureGroups[], TemplateScore }`. Greedy loop that consumes the highest-scoring candidate's feature groups then repeats until exhausted.
**Acceptance:** Unit test: a cluster with two co-eluting proteoforms recovers both envelopes and partitions the feature groups correctly.
**Deps:** T3.9.
**Doc refs:** Algorithm 3c steps 6–7.

---

## Tasks — Step 4: Deconvolution + Quant (→ M8–M12)

### T4.1 — Per-charge averaged spectrum builder [agent-OK]
**Inputs:** A `ChargeStateEnvelope` + its member `FeatureBox`es + the raw files.
**Outputs:** `TopDownEngine/Deconvolution/PerChargeAverager.cs` — for each charge's `FeatureBox`, select the top-N highest-TIC MS1 scans inside its RT range, average via `SpectralAveraging`, return one `MzSpectrum` per charge.
**Acceptance:** Unit test on synthetic data with known peaks produces an averaged spectrum whose peak intensities are ≈ N × single-scan intensity (within noise).
**Deps:** T3.10.
**Doc refs:** Algorithm 4a; Components §5a.

### T4.2 — IsoDec per-charge deconvolution wrapper [agent-OK]
**Inputs:** T4.1 + existing `IsoDecAlgorithm`.
**Outputs:** `TopDownEngine/Deconvolution/PerChargeIsoDec.cs` — wraps `Deconvoluter.Deconvolute(spectrum, IsoDecDeconvolutionParameters)` and returns `List<(double mass, double intensity, double score)>` per charge.
**Acceptance:** Unit test: averaged spectrum of a known proteoform at a known charge yields the known monoisotopic mass within ppm tolerance.
**Deps:** T4.1.
**Doc refs:** Algorithm 4b; Components §5b.

### T4.3 — Parsimony reconciler [agent-OK]
**Inputs:** T4.2.
**Outputs:** `TopDownEngine/Deconvolution/ParsimonyReconciler.cs` — groups per-charge masses within ppm tolerance, counts charge support, emits `ParsimonyCandidateSet[]` sorted by support.
**Acceptance:** Unit test: input with mass `M` recovered at z=5,6,7 produces a candidate with support count = 3.
**Deps:** T4.2.
**Doc refs:** Algorithm 4c; Components §5c.

### T4.4 — Isotopologue basis generator [agent-OK]
**Inputs:** `Averagine` + `IsotopicDistribution` from `Chemistry`.
**Outputs:** `TopDownEngine/Deconvolution/BasisGenerator.cs` — given `(M, z)`, returns an array of `(mz, relative_intensity)` pairs sampled onto a target m/z grid. Also a method to assemble a matrix `B` from many `(M, z)` pairs.
**Acceptance:** Unit test: `(M = 10000, z = 10)` produces the expected averagine envelope peaks at `(M + 10·m_proton) / 10 + n/10` for n = 0…k.
**Deps:** T0.1.
**Doc refs:** Algorithm 4d; Components §5d.

### T4.5 — NNLS solver [needs-review]
**Inputs:** `SharpLearning.Optimization` (already a mzLib dep) — **investigate first** whether it exposes NNLS. If not, implement Lawson-Hanson active-set NNLS (~80 lines).
**Outputs:** `TopDownEngine/Math/Nnls.cs` with `static double[] Solve(double[,] B, double[] y)`.
**Acceptance:**
- Unit tests against known NNLS solutions (textbook examples).
- Benchmark against a reference implementation (e.g., SciPy fixtures checked in as golden data).
- **Review gate:** confirm with user which source (SharpLearning vs. custom) before committing.
**Deps:** T0.1.
**Doc refs:** Algorithm 4d; Components §5d.

### T4.6 — Global NNLS fit + explained fraction [agent-OK]
**Inputs:** T4.3, T4.4, T4.5.
**Outputs:** `TopDownEngine/Deconvolution/EnvelopeFitter.cs` — takes per-charge averaged spectra + parsimony candidates, assembles `B` via T4.4, solves NNLS via T4.5, computes `explained_fraction = ||B·x|| / ||y||`, returns `EnvelopeFitResult { Coefficients, ExplainedFraction, ResidualPeaks }`.
**Acceptance:** Unit test on a synthetic envelope with known species produces `explained_fraction > 0.9` and recovers intensities within 5%.
**Deps:** T4.3, T4.4, T4.5.
**Doc refs:** Algorithm 4d; Components §5d.

### T4.7 — Residual re-parsimony loop [agent-OK]
**Inputs:** T4.6.
**Outputs:** Extension of `EnvelopeFitter`: if `explained_fraction < threshold`, pull residual peaks above noise, re-run parsimony against the residuals, refit. Configurable max iterations (default 2).
**Acceptance:** Unit test: synthetic envelope with two species where parsimony initially finds only one — second pass picks up the second.
**Deps:** T4.6.
**Doc refs:** Algorithm 4d step 5.

### T4.8 — Per-run NNLS quantification [agent-OK]
**Inputs:** T4.6 (locked basis `B`).
**Outputs:** `TopDownEngine/Quantification/PerRunQuantifier.cs` — for each run, compute per-run `y_r` (average top-N scans from that run only, per charge), solve NNLS against the **same `B`**, emit `PerRunQuant { Species → Dict<Run, (intensity, explained_fraction, residual)> }`.
**Acceptance:**
- Unit test: two-run synthetic data with 2× intensity ratio recovers the ratio within 5%.
- Unit test: a run missing a species reports ≈ 0 intensity for it.
**Deps:** T4.6.
**Doc refs:** Algorithm 4e; Components §5e.

### T4.9 — Proteoform DB match [needs-review]
**Inputs:** T4.8.
**Outputs:** `TopDownEngine/Identification/ProteoformMatcher.cs` — loads a protein DB via `ProteinDbLoader`, expands proteoforms (terminal truncations + combinatorial PTMs), sorts by monoisotopic mass, binary-searches each parsimony candidate within ppm tolerance.
**Notes:** Proteoform enumeration is non-trivial — see components doc §6. Defer the full combinatorial PTM expansion to a later iteration; first version handles intact-mass matching only.
**Acceptance:**
- Unit test: known protein + known mass recovers the protein.
- **Review gate:** discuss PTM expansion strategy with user before broad rollout.
**Deps:** T4.8.
**Doc refs:** Algorithm 4f; Components §6.

---

## Tasks — Step 5: MS2 Aggregation (→ M13)

### T5.1 — MS2 scan selector [agent-OK]
**Inputs:** `ChargeStateEnvelope` + raw files.
**Outputs:** `TopDownEngine/Ms2/Ms2ScanSelector.cs` — filters MS2 scans by `MsnOrder == 2`, RT ∈ envelope union range, and `IsolationRange` intersects any charge's m/z window.
**Acceptance:** Unit test on synthetic data with MS2 scans at known m/z recovers only the expected set.
**Deps:** T3.10.
**Doc refs:** Algorithm 5a; Components §5g.

### T5.2 — Three MS2 averaging strategies [agent-OK]
**Inputs:** T5.1 + `SpectralAveraging`.
**Outputs:** `TopDownEngine/Ms2/Ms2Averager.cs` with `AveragePerCharge`, `AveragePerEnvelope`, `AveragePerRegion` methods. All three return `List<(GroupingKey, MzSpectrum averaged)>`.
**Acceptance:** Unit tests for each strategy; verify group counts match expectations.
**Deps:** T5.1.
**Doc refs:** Algorithm 5b; Components §5g.

### T5.3 — Synthetic MS2 `MsDataScan` constructor [agent-OK]
**Inputs:** T5.2.
**Outputs:** Helper that builds `MsDataScan` objects from averaged peak arrays with correct precursor m/z, charge (or charge set), RT, and synthetic scan number.
**Acceptance:** Resulting `MsDataFile` round-trips through `MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra` and back via `MzML` reader with matching peak arrays.
**Deps:** T5.2.
**Doc refs:** Algorithm 5c; Components §5g.

### T5.4 — Envelope sidecar .tsv writer [agent-OK]
**Inputs:** T4.8 (PerRunQuant) + T5.3 (synthetic scan numbers).
**Outputs:** `TopDownEngine/Output/EnvelopeSidecarWriter.cs` using `CsvHelper`. One row per envelope with `M`, charges, RT range, per-run quant, explained fraction, synthetic scan numbers per strategy.
**Acceptance:** Round-trip test: write + read produces the same table.
**Deps:** T4.8, T5.3.
**Doc refs:** Algorithm 5c; Components §5g.

---

## Tasks — Step 6: FDR (→ M14, M15)

### T6.1 — Stage 1 membership scorer [needs-review]
**Inputs:** T4.8 + T3.10.
**Outputs:** `TopDownEngine/Fdr/MembershipScorer.cs` — features: RT deviation, envelope completeness, intra-envelope isotopologue consistency, intensity consistency, per-run explained fraction. Scoring adapted from `FlashLFQ/MBR/MbrScorer` scaffolding.
**Notes:** Feature weights are a research question — ship with equal weights initially, document as TBD.
**Acceptance:** Unit test: synthetic high-confidence envelope scores higher than synthetic low-confidence envelope. **Review gate** on feature set and weight calibration approach.
**Deps:** T4.8.
**Doc refs:** Algorithm 6 Stage 1; Components §7.

### T6.2 — Random-RT null model [agent-OK]
**Inputs:** T6.1.
**Outputs:** Null-model generator that permutes RT labels and re-scores, producing a null distribution for calibration.
**Acceptance:** Unit test verifies the null distribution is flat in p-value under a no-signal input.
**Deps:** T6.1.
**Doc refs:** Algorithm 6 Stage 1.

### T6.3 — Sequence-reversal decoy generation [agent-OK]
**Inputs:** `UsefulProteomicsDatabases/DecoyGeneration/` (already in mzLib).
**Outputs:** Wrapper in `TopDownEngine/Fdr/DecoyGenerator.cs` that produces reversed-sequence decoys at the proteoform level.
**Acceptance:** Unit test: reversed decoy has the same length, composition, and monoisotopic mass as the target.
**Deps:** T0.1.
**Doc refs:** Algorithm 6 Stage 2; Components §7.

### T6.4 — PTM-aware decoy generation [needs-review]
**Inputs:** T6.3.
**Outputs:** Extended decoy generator that preserves PTM count and mirrors PTM positions across the reversed sequence.
**Notes:** **This is an open research question.** Multiple viable strategies exist; pick one with user input, document the choice, and leave the interface pluggable.
**Acceptance:** **Review gate** on strategy choice before implementation.
**Deps:** T6.3.
**Doc refs:** Algorithm 6 Stage 2 + Open Questions.

### T6.5 — Target-decoy FDR calculation [agent-OK]
**Inputs:** T6.1, T6.3 (or T6.4).
**Outputs:** `TopDownEngine/Fdr/TwoStageFdr.cs` — combines Stage 1 scores with Stage 2 target-decoy competition, emits `FdrResult { QValue, PEP }` per envelope.
**Acceptance:** Unit test on synthetic target+decoy score distributions reproduces a known q-value curve.
**Deps:** T6.1, T6.3.
**Doc refs:** Algorithm 6.

---

## Tasks — Testing Suite (→ M16–M20)

See `TopDown-Engine-Testing.md` for full context. Tasks here are summary pointers.

### T16.1 — Real-data ingestion [agent-OK]
**Outputs:** Parser pipeline for the two raw files plus the supplied search-tool results. Uses existing `Readers/ExternalResults/ResultFiles/TopPICSearchResultFile.cs` and `ToppicPrsm.cs`. MetaMorpheus parser verified; ProSight and pTop parsers deferred per testing doc.
**Deps:** T0.2.
**Doc refs:** Testing §A.

### T17.1 — Consensus feature extractor [agent-OK]
**Outputs:** `TopDownEngine.Test/ConsensusFeatureExtractor.cs`. k-of-n intersection of search results within ppm + RT tolerance.
**Deps:** T16.1.
**Doc refs:** Testing §B.

### T18.1 — Empirical characterization module [agent-OK]
**Outputs:** Measures distributions per Testing §C (XIC widths, isotopologue depth, charge breadth, noise floor, neighborhood density, ppm drift, scan correlation). Emits a parameter .json.
**Deps:** T17.1, T1.2.
**Doc refs:** Testing §C.

### T19.1 — MVP simulator [needs-review]
**Outputs:** Parametric simulator per Testing §D. Gaussian RT curves, averagine envelopes, Gaussian noise, fixed ppm jitter. Emits synthetic `.mzML` + ground-truth manifest.
**Notes:** Realism target is a user decision — parametric vs. resampling. Start parametric. **Review gate** on realism criteria before the engine is tuned against it.
**Deps:** T18.1.
**Doc refs:** Testing §D.

### T20.1 — Engine test harness [agent-OK]
**Outputs:** Runner that invokes the engine on simulated + real data, reports recovery metrics per Testing §E.
**Deps:** T19.1, and enough engine tasks to run the pipeline (at minimum T4.8).
**Doc refs:** Testing §E.

---

## Dispatch Notes for Subagents

When dispatching a subagent to a task:

1. **Brief it with the task block verbatim plus the relevant doc section.** The task blocks are written to be self-contained, but the doc sections carry the "why" — hand both over.
2. **Mention the project layout conventions.** The engine is a new project at `mzLib/TopDownEngine/`, tests at `mzLib/TopDownEngine.Test/`, all builds run from `mzLib/`. Target framework `net8.0` for the engine, `net8.0-windows` for tests.
3. **Flag the MetaMorpheus integration CI.** Any change that touches existing mzLib public APIs triggers MetaMorpheus rebuild in CI. Subagents should default to adding new types in the new engine project, not modifying mzLib proper, unless the task explicitly requires it.
4. **Tell it which doc section to consult.** The `Doc refs` line on each task points at the relevant algorithm + components sections. Handing those over lets the subagent resolve design questions without guessing.
5. **Require acceptance-criteria tests.** Every task block's `Acceptance` line is the definition of done. A subagent that returns code without the acceptance tests is not done.
6. **`[needs-review]` tasks stop at draft.** A subagent assigned a `needs-review` task should produce the implementation and tests, then hand it back for user review before landing. Do not let a subagent make the scientific call itself.

---

## Task Dependency Graph (Top-Level)

```
T0.1 ──► T0.2 ──► T16.1 ──► T17.1 ──► T18.1 ──► T19.1 ──► T20.1
  │
  └──► T1.1 ──► T1.2 ──► T2.1 ──► T2.2 ──► T2.3 ──► T2.4
                   │
                   ├──► T3.1 ──► T3.2 ──┐
                   │                    ├──► T3.5 ──► T3.6 ──► T3.8 ──► T3.9 ──► T3.10
                   │       T3.3 ──► T3.4┘          T3.7 ┘
                   │                                        │
                   │                                        ├──► T4.1 ──► T4.2 ──► T4.3
                   │                                        │                       │
                   │                       T4.4 ──┐         │                       │
                   │                       T4.5 ──┼──► T4.6 ──► T4.7                │
                   │                              ┘     │                           │
                   │                                    ├──► T4.8 ──► T4.9          │
                   │                                    │     │                     │
                   │                                    │     └──► T5.1 ──► T5.2 ──► T5.3 ──► T5.4
                   │                                    │
                   │                                    └──► T6.1 ──► T6.2
                   │                                         T6.3 ──► T6.4 ──► T6.5
```

The critical path through the engine is: **T0.1 → T1.1 → T1.2 → T2.1 → T2.2 → T2.3 → T3.1 → T3.2 → T3.5 → T3.6 → T3.8 → T3.10 → T4.1 → T4.2 → T4.3 → T4.6 → T4.8**, which gets you end-to-end from raw data to a quant table. Everything else (T5.x, T6.x, T4.7, T4.9, and the testing suite) can be parallelized or deferred once the critical path is proven on synthetic data.
