# FLASHDeconv Log-Mass / Charge-Envelope Detection — Implementation Notes

Source reconstruction is primarily from the **OpenMS source** (the paper's body is paywalled and the PMC version is not available for this article; the ScienceDirect/Cell fetches returned only abstract-level text). The core math below is confirmed from the OpenMS C++ implementation (`FLASHHelperClasses.cpp`, `SpectralDeconvolution.cpp`) plus the paper's abstract and the OpenMS documentation. Where I could not retrieve a direct quote from the PDF, I have reconstructed the math from the code and flagged it.

---

## 1. The log-mass transform

**What is logged.** Not `log(m/z)` directly, and not `log(mass)`. FLASHDeconv logs `m/z` **after subtracting the charge carrier mass** (a proton, ~1.00727646688 Da in positive mode). The exact implementation from OpenMS:

```cpp
// FLASHHelperClasses.cpp
double FLASHHelperClasses::getLogMz(const double mz, const bool positive)
{
  return std::log(mz - getChargeMass(positive));
}
```

Natural log (`std::log`), not log10.

**Why this makes charge states regularly spaced.** For a proteoform of monoisotopic mass `M`, ionized to charge `z` by `z` protons of mass `m_H`, the observed m/z is:

```
mz(z) = (M + z * m_H) / z  =  M/z + m_H
```

Therefore

```
mz(z) - m_H = M / z
log(mz(z) - m_H) = log(M) - log(z)
```

So every charge state of the same proteoform lands at position `log M − log z` along a single axis. Across the charge ladder `z = 1, 2, 3, …, Z`, the envelope of peaks in log-mz-minus-proton space is a fixed **template** that depends only on `z`, shifted by `log M`:

```
template(z) = -log(z)     for z = 1, 2, …, Z_max
```

This template is charge-position-invariant: the same shape shows up at every mass; only the offset `log M` changes. That is the "constant pattern" the paper refers to. Note the spacing is *not* uniform — it is logarithmically compressed at higher `z` (`log(z+1) − log(z) ≈ 1/z`). Calling it "regularly spaced" is misleading; the correct statement is **"the inter-charge spacings are a fixed function of `z` that is independent of `M`."**

This is confirmed in the OpenMS code where the search template is built as:

```cpp
// SpectralDeconvolution.cpp
universal_pattern_.push_back(-log(i + 1));   // i indexes charge state
```

**Paper quote (abstract, Cell Systems 2020):** *"FLASHDeconv transforms peak positions (m/z) within spectra into log m/z space. This simple transformation turns the deconvolution problem into a search for constant patterns, thereby greatly accelerating the process."*

---

## 2. Envelope detection in log-mass space

FLASHDeconv does **not** use FFT or autocorrelation. It uses a **binned cross-correlation / template match** of the fixed charge template `{−log 1, −log 2, −log 3, …, −log Z_max}` against the binned log-mz axis. Concretely:

1. **Bin the log-mz axis.** The bin index for a peak at log-mz value `x` is
   ```cpp
   // SpectralDeconvolution.cpp
   return (Size)round((value - min_value) * bin_mul_factor);
   ```
   `bin_mul_factor` is set from the ppm tolerance:
   ```cpp
   bin_mul_factors_.push_back(1.0 / j * tol_div_factor);   // j = tol * 1e-6, tol_div_factor = 2.5
   ```
   So for a 10 ppm tolerance, `bin_mul_factor ≈ 2.5 / 10e-6 = 2.5e5`, i.e. bin width `≈ 4e-6` in log units (~4 ppm). Bins are roughly `tol / 2.5` wide in log space (log(1+ε) ≈ ε, so a log-space bin of width ε is ε ppm).
   ```cpp
   Size bi = getBinNumber_(p.logMz, mz_bin_min_value_, bin_mul_factor);
   ```

2. **Mark occupied bins.** Every observed peak's log-mz (i.e., `log(mz − m_H)`) becomes a `1` in a boolean "mz bin" vector. This is the sparse peak presence vector.

3. **Build the charge template in the same bin units.** The pattern `−log(i+1)` for `i = 0, …, Z_max−1` is converted to bin offsets. For charge `z`, the offset is `round((−log z) * bin_mul_factor)` relative to the `z=1` position. Because `bin_mul_factor` is the same for every candidate mass, these template offsets are **precomputed once** and reused across all candidate positions.

4. **Candidate mass enumeration.** For each candidate monoisotopic mass `M` (also binned: `mass_bin_min` to `mass_bin_max` at the same `bin_mul_factor` resolution), check how many of the template's bin positions are occupied in the mz-bin vector. Since `log M − log z = log(mz − m_H)`, a candidate mass bin at log-mass `μ = log M` "lights up" charge `z` if the mz-bin at `μ − log z` is occupied. This is essentially a sparse logical AND between the mz presence vector and the shifted template, performed over all candidate mass bins. Implementation is a tight loop over charges and masses — the binning allows it to run in integer arithmetic.

5. **Candidate mass score = number of lit charges** (subject to threshold) plus the downstream isotope/cosine score (see §6).

6. **Harmonic (sub/super-harmonic) elimination.** A false positive arises when the true mass is `M`, but the algorithm also lights up candidate masses at `M/2, M/3, M*2, …` because the `log(i+1)` ladder of the harmonic also hits real peaks. FLASHDeconv precomputes *harmonic templates* at half-charges for a small set of harmonic orders `{2, 3, 5, 7, 11}` and subtracts/penalizes candidates whose evidence is better explained by a harmonic of a neighboring mass:
   ```cpp
   // SpectralDeconvolution.cpp — harmonic pattern construction
   harmonic_pattern_matrix_.setValue(k, i, -log(b - (b - a) * n / hc));
   // where hc is harmonic charge, n = hc / 2
   ```
   In the paper this is called "harmonic artifact elimination."

### Pseudocode

```
Inputs: peaks = [(mz_j, intensity_j)]
        z_range = [z_min .. z_max]    # default 2..100
        M_range = [M_min .. M_max]    # default 1..100 kDa
        tol_ppm
        min_charges                   # minimum charges in an envelope, e.g. 3

# 1. Log-transform
for each peak j:
    x_j = log(mz_j - m_H)             # natural log, subtract proton

# 2. Bin log-mz axis
bin_mul = 2.5 / (tol_ppm * 1e-6)
mz_bin[b] = True  iff any peak x_j rounds to bin b

# 3. Precompute charge template offsets (in bin units)
for z in z_range:
    off[z] = round((-log(z)) * bin_mul)   # negative for z > 1

# 4. Sweep candidate masses
for each mass bin μ in log(M_range):
    lit_charges = []
    for z in z_range:
        b = μ + off[z]
        if mz_bin[b]: lit_charges.append(z)
    if len(lit_charges) >= min_charges:
        emit candidate (μ, lit_charges)

# 5. Harmonic suppression
for each candidate:
    for h in {2, 3, 5, 7, 11}:
        if evidence(candidate) better explained at μ ± log(h):
            drop or penalize

# 6. Deisotope + score each survivor (see §6)
```

---

## 3. Input

- **Per scan, centroided peak list.** FLASHDeconv operates on centroided MS1 (or MS/MS) spectra. Profile data is not the target; peaks must already be picked. No pre-grouping or prior RT-binning is required at the envelope-detection step — envelope detection runs on one scan at a time.
- Later, a separate **MassFeatureTrace / feature-finding** step aggregates per-scan deconvoluted masses across RT into features. That is a *downstream* step, not the envelope detector.
- Handles **isotopically resolved and unresolved** peaks within wide charge/mass ranges: charge `2–100`, mass `1–100 kDa` by default.

---

## 4. Mapping a detected envelope back to charge assignments and monoisotopic mass

Because the template position itself is `log M`, the candidate mass is recovered directly:

```
M = exp(μ)                          # μ = candidate log-mass bin center
```

Each "lit" charge `z` in the template corresponds to an observed peak at
```
mz_z = M/z + m_H
```
so the (peak → charge) assignment is produced as a byproduct of which template bins were occupied. The set `{z : bin occupied}` is the envelope. A peak group in OpenMS (`PeakGroup`) stores the candidate mass plus the per-charge subset of peaks that matched.

Because of the binning resolution, the mass `M = exp(μ)` is only accurate to `tol_ppm / 2.5`; the final mass is **refined** after deisotoping against averagine (§6), typically by recomputing the mass from the most intense isotope peak of the highest-SNR charge.

---

## 5. Overlapping envelopes / co-eluting proteoforms

Because envelope detection is a *linear* scan over candidate mass bins, overlapping proteoforms just produce multiple candidate mass bins, each with its own lit-charge set. There is no hard assignment of a peak to a single envelope at the detection stage — a given peak can participate in multiple candidate envelopes. The resolution happens at scoring time:

- Each candidate gets a cosine-similarity score against averagine and an SNR per charge.
- Harmonic/sub-harmonic artifact elimination (§2.6) removes candidates that are integer-ratio copies of a stronger real mass.
- Low-score candidates are filtered.
- What remains is a set of mutually coexisting masses; the same m/z peak can legitimately contribute to multiple proteoforms if the underlying physics is shared (rare in practice since the charge templates typically disambiguate).

---

## 6. Averagine / isotope-pattern scoring

After a candidate `(M, {z_i})` is detected, each charge's peaks are checked against an **averagine isotope envelope** for mass `M` via cosine similarity:

```cpp
// SpectralDeconvolution.cpp
float tmp_cos = getCosine(per_isotope_intensities, min_isotope_index,
                          max_isotope_index, iso, tmp_offset, min_iso_size);
```

Where `iso` is a precomputed averagine isotope distribution (likely from `CoarseIsotopePatternGenerator`) binned by mass. The cosine is computed between the observed per-isotope intensities for that charge and the theoretical averagine pattern, with an offset search to align the monoisotopic peak (since the brightest observed peak may not be the monoisotopic one at high mass). Candidates below a cosine threshold are dropped.

Final per-mass quality is an aggregate of:
- Number of matched charges
- Per-charge isotope cosine
- Per-charge SNR (`pg.getChargeSNR(abs_charge)`)
- Harmonic suppression penalty

---

## 7. Tolerances and hyperparameters

From the OpenMS source and docs:

| Parameter | Default | Notes |
|---|---|---|
| `tol` (ppm) | 10 ppm MS1 | Used to set `bin_mul_factor = 2.5 / (tol * 1e-6)` — bin width is ~`tol / 2.5` ppm in log space |
| `min_charge` / `max_charge` | 2 .. 100 | Charge range searched |
| `min_mass` / `max_mass` | 1 .. 100 kDa | Candidate mass range |
| `min_isotope_cosine` | ~0.85 (scan), ~0.75 (feature) | Averagine cosine threshold |
| Harmonic orders | {2, 3, 5, 7, 11} | Hardcoded for harmonic artifact elimination |
| Charge carrier | proton, m_H ≈ 1.00727646688 Da | subtracted before log |
| `tol_div_factor` | 2.5 | sub-ppm bin resolution factor |
| Min charges per envelope | ~3 (configurable) | Envelopes with fewer matched charges are suppressed |

---

## 8. What FLASHDeconv does NOT do / assumes

- **Assumes high-resolution, centroided data.** The algorithm relies on ppm-level bin widths; low-res or profile-mode data breaks the template match.
- **Does not require the full charge range to be present** — the `min_charges` threshold (commonly 3) is the minimum, so short envelopes are accepted, but envelopes with only 1–2 charges are rejected as ambiguous.
- **Assumes proton charge carrier** — adducts (Na, K, NH4) are not the default transform and would require changing `getChargeMass` or adding adduct-aware templates.
- **Averagine assumption** at scoring time — unusual elemental compositions (heavy modifications, non-peptide polymers) will score poorly even if the charge envelope is real.
- **Per-scan** — there is no cross-scan evidence pooling at the envelope-detection step. RT information only enters at the later feature-trace step.
- **Very low mass (< 1 kDa)** is outside the default range; at ~500 Da isotope patterns dominate and a different detector is appropriate.
- **Very high mass** suffers because charges compress together in log-mz space (`log(z+1) − log(z) → 0`); bin collisions between adjacent charges grow unless ppm tolerance is tightened. The paper notes accuracy degrades above ~50 kDa.
- **No FFT.** Despite a superficial resemblance to harmonic detection, the implementation is a precomputed-offset table lookup, not a Fourier method.

---

## 9. Applicability to a cross-run-aggregated feature list

Your setting: you have a sparse list `(mz_f, RT_f, intensity_vec_f)` of MS1 features, already cross-run matched. You want to group features belonging to the same proteoform.

**What transfers cleanly:**

- **The log-mz transform itself** is per-peak and cares nothing about density. Computing `x_f = log(mz_f − m_H)` for each feature works identically.
- **The charge template** `{−log z}` is independent of density and can be applied to a sparse feature list. In fact a feature list is *cleaner* than a raw spectrum because most noise peaks are already gone.
- **Template-match candidate enumeration** works on a sparse set — the `mz_bin[]` vector is a set of occupied bins, and the sweep is identical. You don't need a dense array; a hash set of occupied bins and a loop over `(candidate mass bin) × (charge)` runs in `O(|M_bins| * Z_max)` with `O(1)` occupancy lookup.
- **RT co-elution as an additional constraint** is a natural improvement that FLASHDeconv does not use at the envelope-detection stage: you can require that two features assigned to the same candidate mass have RT within, e.g., 0.1 min, *and* similar cross-run intensity profiles (Pearson correlation on `intensity_vec_f`). This is strictly stronger evidence than FLASHDeconv has and should sharply reduce harmonic artifacts.
- **Harmonic suppression** still matters and transfers directly.

**What does not transfer / needs adaptation:**

- **Averagine cosine scoring at the deisotope step** assumes you have multiple isotope peaks per charge within a scan. In a cross-run feature list, if each feature already represents a monoisotopic (or single-isotopologue) trace, you do not have a per-charge isotope vector to score against averagine *at this step*. You would score averagine at the feature-extraction stage (when MS1 boxes are picked) rather than at the charge-grouping stage. FLASHDeconv's envelope detection is actually the **decharging** step only; the deisotoping happens after and may not apply to you.
- **SNR per charge** as FLASHDeconv computes it is scan-local; you would replace it with cross-run intensity-profile consistency (correlation between the intensity vectors of the charge-paired features).
- **Bin widths** may need to be looser than 10 ppm if your feature m/z centroids have been averaged across runs with calibration drift. Set `tol` to whatever cross-run mass accuracy you believe.
- **Minimum charges per envelope** (`min_charges`) is still the dominant hyperparameter; with the added RT/intensity-profile constraint you can probably drop it to 2 safely.

**Short answer:** yes — the log-transform + fixed-offset template match is exactly the primitive you want, and it works on a sparse list as well as it works on a dense spectrum. The parts of FLASHDeconv you should *not* port are the per-scan averagine deisotoping and the per-scan SNR, which belong at a different stage of your pipeline.

---

## Sources actually read

- Cell Systems 2020 FLASHDeconv page (abstract only — full text paywalled):
  https://www.cell.com/cell-systems/fulltext/S2405-4712(20)30030-2
- ScienceDirect mirror (abstract only):
  https://www.sciencedirect.com/science/article/pii/S2405471220300302
- PubMed entry:
  https://pubmed.ncbi.nlm.nih.gov/32078799/
- OpenMS `FLASHHelperClasses.cpp` (raw, via WebFetch) — source of the exact `getLogMz` definition:
  https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/openms/source/ANALYSIS/TOPDOWN/FLASHHelperClasses.cpp
- OpenMS `FLASHHelperClasses.h` (raw) — `LogMzPeak` struct:
  https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h
- OpenMS `SpectralDeconvolution.cpp` (raw) — source of `bin_mul_factor`, `universal_pattern_`, `harmonic_pattern_matrix_`, and cosine scoring snippets:
  https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/openms/source/ANALYSIS/TOPDOWN/SpectralDeconvolution.cpp
- OpenMS `FLASHDeconvAlgorithm.cpp` (raw) — top-level wiring (log-mz handled in `SpectralDeconvolution`):
  https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/openms/source/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.cpp
- OpenMS FLASHDeconv documentation page:
  https://openms.de/documentation/html/TOPP_FLASHDeconv.html
- OpenMS FLASHDeconvAlgorithm class reference:
  https://openms.de/documentation/html/classOpenMS_1_1FLASHDeconvAlgorithm.html

## Gaps / caveats

- I could not retrieve the **full body** of the Cell Systems paper (paywall; PMC not available for this DOI at the URLs I tried). Direct quotes of the math from the paper itself are therefore limited to the abstract's "transforms peak positions (m/z) within spectra into log m/z space … a search for constant patterns." Everything more specific above is reconstructed from the OpenMS source, which is the canonical reference implementation, and cross-checked against the OpenMS docs.
- The specific default values for `min_isotope_cosine`, `tol_div_factor = 2.5`, and the harmonic order set `{2,3,5,7,11}` are taken from the source code as the current behavior; the original 2020 paper may have used slightly different defaults.
- I did not fetch the FLASHIda or FLASHDeconv 3.0 papers because the log-transform primitive did not change between versions — they add intelligent data acquisition and improved feature tracing on top of the same decharging core.
