# Top-Down Engine: Testing and Simulation Suite

Preliminary/parallel work stream for the top-down search engine described in `TopDown-Engine-Algorithm.md` and `TopDown-Engine-Components.md`. This document is **not** an algorithm step — it is the testing and ground-truth infrastructure we need in place before (and alongside) the engine itself is built.

---

## Context

The top-down engine defers deconvolution, groups cross-run features in raw m/z space, then does per-charge averaging + IsoDec + parsimony + NNLS to recover proteoforms. Every step of that pipeline — feature grouping, RT alignment, charge-state envelope assembly, per-charge deconvolution, parsimony, NNLS quant — has its own failure modes. Without ground truth, we can't tell whether a bad result means the algorithm is wrong, the parameters are wrong, or the data is noisier than expected.

**The problem with "just test on real data":** Real top-down MS1 spectra contain proteoforms we don't know about. Any search tool's ID list is incomplete and partially wrong, so "real data + a search tool's result" is a noisy partial ground truth. It's good for *validation* but bad for unit-level *development* feedback.

**The problem with "just test on fake data":** Fake data that doesn't look like real data tests only the toy cases. We need simulated data whose noise, peak density, envelope shape, isotopologue depth, charge-state breadth, and RT profiles match the distributions in actual top-down runs, or the engine will be tuned for a world it won't encounter.

**The plan:**
1. Load Alex's two real raw files plus the search results from several different top-down tools (TopPIC, ProSight, MetaMorpheus, pTop, …).
2. Take the intersection — proteoforms that **all** (or almost all) tools identified — as a high-confidence feature set.
3. Characterize that set empirically: RT widths, envelope completeness, charge-state breadth, isotopologue depth, intensity distributions, noise floor around each feature, neighborhood peak density.
4. Build a parametric simulator that, given a target proteoform list and empirical noise parameters, emits synthetic mzML files that look indistinguishable from the real ones at the distribution level.
5. Use the simulator as the ground-truth dev dataset throughout engine development. Use the real files + search-tool consensus as the final validation check.

---

## Resources We Already Have

- **Two real .raw files** (Alex's local data). Thermo format — readable via `Readers/ThermoRawFileReader`.
- **Search results from multiple top-down tools.** Formats vary — TopPIC emits `.csv`/`.txt`, ProSight has its own, MetaMorpheus writes TSV, pTop has its own. Each reports proteoform IDs with monoisotopic mass, charge states observed, RT, and (usually) some fragment match info.
- **mzLib reader infrastructure** — `Readers/`, `MsDataFile`, `MsDataScan`, scan iteration, peak arrays, isolation windows, RT lookup.
- **Feature extraction primitives** — `PeakIndexingEngine` + `ExtractedIonChromatogram` for building XICs at known m/z.
- **SpectralAveraging** for averaging scans.
- **mzML writing** via `Readers/MzML/MzmlMethods.cs` and `SpectralAveraging/AveragedSpectraWriter.cs` — needed for the simulator's output.

## Resources We Need to Build

- Start with MetaMorpheus results only.
- Add ProSightPD once its SQL-backed parser exists.
- Parsers for each additional top-down search tool's result format (if we don't already have one that covers the tools on hand).
- A consensus-feature extractor that intersects IDs across tools at the proteoform level, tolerating mass/RT tolerances.
- An empirical characterization module (histograms + fit distributions).
- The simulator itself.
- A test harness that runs the engine against simulated data and reports recovery rate, precision, and quantitative accuracy.

---

## Testing Suite: Components

### A. Real-Data Ingestion

**Goal.** Load the two raw files and all search-tool results into a common in-memory representation.

| Need | Status | Notes |
|---|---|---|
| Raw file reader | **Reuse** | `Readers/ThermoRawFileReader` |
| TopPIC result parser | **Reuse** | `Readers/ExternalResults/ResultFiles/TopPICSearchResultFile.cs` and `Readers/ExternalResults/IndividualResultRecords/ToppicPrsm.cs` already exist. Keep tests focused on contract coverage and edge cases. |
| ProSight result parser | **Build** | ProSight has its own format; no existing mzLib parser. |
| MetaMorpheus result parser | **Reuse or adapt** | mzLib has `Readers/ExternalResults/` coverage for MetaMorpheus-style outputs; verify the exact fields needed for top-down consensus. |
| pTop result parser | **Build** | Add if pTop is one of the tools on hand. |
| Unified `TopDownSearchResult` record | **Build (small)** | Normalize across tools: `{ MonoisotopicMass, Charges[], RtApex, RtRange, Sequence, Ptms[], Score, SourceTool }` |

**Test expectations for this layer.**
- Parser tests must prove we can round-trip or normalize each supported result format into `TopDownSearchResult` without silently dropping fields.
- Fixture tests should include one small gold file per parser and one malformed file per parser.
- Comparison tests should assert tolerance-based equivalence, not string equality, for mass and RT fields.

### B. Consensus Feature Extraction

**Goal.** Find the set of proteoforms that all (or ≥ k out of n) search tools identified on the same raw file, within a ppm mass tolerance and an RT tolerance. These are the high-confidence features we'll use to characterize the real data.

| Need | Status | Notes |
|---|---|---|
| Cross-tool ID intersection | **Build** | Group `TopDownSearchResult` records by (mass within ppm, RT within seconds); keep groups with support from ≥ k tools. |
| Per-group RT/mass/charge reconciliation | **Build** | Take median mass, union of reported charges, intersection or union of RT ranges (configurable) |
| Emit as a `ConsensusFeature` list | **Build (small)** | Container with the ground-truth proteoform info for one raw file |

Use this list as both (a) the empirical characterization input and (b) the "known answer key" for validating engine output against real data.

**Test expectations for this layer.**
- Exact-group tests on tiny synthetic inputs: same mass/RT within tolerance must cluster; borderline cases must not.
- Stability tests: reordering inputs must not change consensus membership.
- Support-threshold tests: `k-of-n` behavior must be explicit and parameterized.
- Conflict tests: discordant charge lists, PTM annotations, or RT ranges must reconcile deterministically.

### C. Empirical Characterization

**Goal.** For each consensus feature, measure the MS1 signal in the raw file around it and record the distributions we need to reproduce in the simulator.

| Measurement | Why it matters | How |
|---|---|---|
| **Per-charge XIC width (FWHM, asymmetry)** | Drives the RT axis of the flood-fill | Extract XIC at each known charge state's m/z, fit peak boundaries with `ExtractedIonChromatogram.FindPeakBoundaries` |
| **Isotopologue depth per charge** | Drives whether we see full envelopes or just the brightest few peaks | Count how many isotopologues are above noise at each charge |
| **Charge-state breadth per proteoform** | Drives how wide the log-mass template match has to sweep | Count distinct charges observed per consensus feature |
| **Per-charge envelope intensity ratio** | Drives the NNLS basis realism | Normalize intensities across charge states for each feature; build an empirical distribution |
| **Baseline noise floor around features** | Drives the `numPeakThreshold` / intensity stop in `GetAllXics` | Measure median peak intensity in an m/z neighborhood around each feature, excluding the feature itself |
| **Neighborhood peak density** | Drives how many false flood-fill seeds we have to discard | Count peaks per m/z unit in the neighborhood |
| **Mass accuracy / ppm drift per scan** | Drives the binomial null model for cross-run matching | Compare observed centroid m/z vs. theoretical across scans; fit ppm error distribution |
| **Scan-to-scan correlation** (signal stability within one envelope) | Drives the "real feature vs. noise" test | Measure XIC smoothness for known features vs. neighborhood |

Output: a JSON or .tsv of distributions + summary statistics that the simulator reads as its parameter set.

**Test expectations for this layer.**
- Each measured statistic should have a unit test with a handcrafted feature and a known expected value.
- Distribution outputs should be reproducible under fixed seed and insensitive to input ordering.
- Summary files should validate schema, units, and non-negativity constraints.

### D. Simulator

**Goal.** Given a target proteoform list and the empirical parameter set from C, emit a synthetic `.mzML` that looks distributionally indistinguishable from the real files.

Inputs:
- Target proteoform list: `[(sequence or formula, abundance, RT apex, RT width, charge distribution)]`
- Empirical parameter set from C
- Optional: random seed, instrument resolution, desired number of scans

Outputs:
- Synthetic `.mzML` written via `MzmlMethods` / `AveragedSpectraWriter`
- Ground-truth manifest (.tsv / .json): per-proteoform, per-scan, which peaks came from which proteoform at what intensity

Building blocks:

| Need | Status | Notes |
|---|---|---|
| Predicted isotopologue pattern per proteoform | **Reuse** | `Averagine` / `IsotopicDistribution` |
| RT elution curve generator (Gaussian / EMG) | **Build (small)** | Parametric curve; parameters drawn from distribution C |
| Per-scan peak-list assembly | **Build** | For each simulated scan time, sum contributions from all active proteoforms at their current elution intensity |
| Noise injection (baseline + per-peak) | **Build** | Poisson or Gaussian noise with parameters from C; also inject random noise peaks at the observed density |
| Mass accuracy perturbation | **Build (small)** | Per-peak ppm jitter with distribution from C |
| Synthetic `MsDataScan` construction | **Reuse pattern** | Same pattern as Step 3g's synthetic MS2 output — build `MsDataScan` from peak arrays with correct metadata |
| Synthetic `MsDataFile` + mzML write | **Reuse** | `MzmlMethods` writer |
| Ground-truth manifest writer | **Build (small)** | CsvHelper-based .tsv with peak provenance |

**Test expectations for this layer.**
- Seeded determinism: same seed, same mzML and manifest.
- Conservation: simulated peak intensity should equal the sum of proteoform contributions plus injected noise within a bounded tolerance.
- Envelope realism: isotope spacing, charge breadth, and RT width distributions should match the empirical parameter file within tolerance bands.
- Manifest integrity: every simulated peak should be traceable to a proteoform or marked as noise.
- Failure-mode fixtures: absent RT curves, zero-abundance proteoforms, and extreme charge states should fail cleanly.

### E. Engine Test Harness

**Goal.** Run the engine against synthetic (and real) data and report recovery metrics.

| Metric | What it measures |
|---|---|
| **Feature recovery rate** | Fraction of ground-truth proteoforms that the engine finds a matching `ChargeStateEnvelope` for |
| **Mass accuracy** | ppm error between engine-reported `M` and ground-truth `M` |
| **Charge-state recovery** | Fraction of ground-truth charges recovered per proteoform |
| **Quant accuracy** | Regression between engine `PerRunQuant` and ground-truth abundance |
| **False discovery rate** | Fraction of engine-reported envelopes that don't match any ground-truth proteoform |
| **Envelope explained fraction** | Engine-reported NNLS explained fraction — validates whether the internal quality metric correlates with actual correctness |

The harness should run on **both** the simulated data (where ground truth is exact) and the real data filtered to the consensus features (where ground truth is partial but independent from the engine). Disagreement between "works on sim, fails on real" points to simulator-realism gaps; disagreement between "works on real, fails on sim" points to simulator drift from the data's actual distribution.

**Test expectations for this layer.**
- Stage-gated tests: feature recovery, then charge recovery, then quantification, then real-data sanity.
- Metric assertions should use thresholds, not exact counts, except for tiny fixtures.
- Each run should emit a machine-readable metrics file plus a human-readable summary.
- Regression tests should pin a small set of synthetic fixtures so algorithm changes cannot silently degrade recovery.

---

## Proposed Test Taxonomy

1. **Parser contract tests**
- Verify each supported search-tool file loads, normalizes, and round-trips.
- Verify malformed headers, missing fields, and duplicated alternative IDs fail predictably.

2. **Consensus logic tests**
- Use tiny synthetic result sets to prove mass/RT clustering, support thresholds, and deterministic reconciliation.
- Include one real-data smoke test per file that only checks the consensus feature count and field validity.

3. **Characterization tests**
- Validate the per-feature measurement routines against hand-computed examples.
- Validate histogram binning, quantiles, and fitted parameters.

4. **Simulator tests**
- Treat the simulator like a compiler: given a fixed seed and input manifest, the output must be reproducible.
- Add statistical tests for the generated distributions, plus provenance checks for the manifest.

5. **Engine acceptance tests**
- Run the full pipeline on simulated fixtures and assert recovery/precision/quant thresholds.
- Run a narrow real-data validation set and assert the engine stays within consensus-derived bounds.

6. **Performance/regression tests**
- Pin runtime and memory envelopes on one medium synthetic dataset so the test suite catches accidental quadratic behavior.
- Record the number of recovered features, envelopes, and quantifiable species for a known seed.

---

## Why This Is a Separate Document

- **Parallel to algorithm development.** The engine can be prototyped against toy synthetic data while the testing suite catches up, but production-quality development needs the real testing suite in place.
- **Different file layout.** The testing suite's parsers, characterization scripts, simulator, and test harness aren't part of the engine's code path — they're a sibling project (possibly a new `Test/TopDownEngine/` subfolder plus a small simulator library).
- **Different iteration cadence.** The simulator parameters will get retuned as we learn more about the real data; that's independent from algorithm changes.

---

## Prototyping Order

1. **Real-data ingestion + consensus feature extraction** (A + B). Cheapest and unblocks everything else. If we can't get consistent consensus features across the on-hand search tools, the ground-truth premise is shaky and we should re-scope before building the simulator.
2. **Empirical characterization** (C). Builds the distribution parameter set from the consensus features. Output is a static file that the simulator consumes.
3. **Simulator — minimum viable version** (D). First cut: Gaussian RT curves, averagine envelopes, Gaussian noise, fixed ppm jitter. Enough to start running the engine against a known answer key. Refine the realism as the test harness surfaces gaps.
4. **Engine test harness** (E). Runs in parallel with D — first version uses the MVP simulator, later versions validate against real data with the consensus feature set as the partial key.
5. **Realism refinements** in the simulator (return to D) as the harness flags distributional mismatches.

## Acceptance Gates

- Parser layer: all supported formats load and normalize on the first try.
- Consensus layer: membership is stable under reordering and tolerance boundaries are explicit.
- Simulator MVP: seeded output is deterministic and distributionally plausible.
- Harness MVP: synthetic fixtures recover the majority of ground-truth proteoforms with bounded mass and quant error.
- Real-data gate: consensus-filtered real files produce sane counts, no crashes, and metrics that are directionally consistent with the simulator.

---

## Open Questions

- **Which search tools' results are available?** The intersection-quality depends heavily on which tools and how concordant they are. If we only have two tools that disagree a lot, "consensus" may need to mean "either tool with high score + manually curated."
- **Do we trust the search tools' reported charge-state lists, or re-derive them from the raw data around the reported mass?** Different tools report charges differently; re-deriving from XIC evidence is safer.
- **How much RT drift between the two raw files?** Affects whether the consensus extraction has to do its own alignment before intersection.
- **Simulator realism target.** Parametric (fit distributions) or non-parametric (resample peaks from real data)? Parametric is easier to control; resampling is more realistic but harder to debug. Start parametric, add resampling as an option if needed.
- **Where does the simulator live?** Inside mzLib as a new project (`TopDownSimulation`), or in the new top-down engine project? Probably the latter, to avoid bloating mzLib with a test-only dependency.
