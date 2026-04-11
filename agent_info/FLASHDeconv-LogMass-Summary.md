# FLASHDeconv Log-Mass Trick — Condensed Summary

Short, actionable summary of the FLASHDeconv charge-state envelope detection primitive, for use in Step 2c of the top-down engine. For full math, pseudocode, source-code citations, and caveats, see `FLASHDeconv-LogMass-Notes.md`.

## Key findings

- **Exact transform: `log(m/z − m_proton)`, natural log.** Not `log(mz)`, not `log(mass)`. The proton subtraction is what collapses the algebra to `log M − log z`. Confirmed from OpenMS `FLASHHelperClasses.cpp` (`getLogMz`).

- **Charge states are NOT uniformly spaced in log space.** The paper's "constant pattern" language is slightly misleading. The template is `{−log 1, −log 2, −log 3, …}`; inter-charge spacing is `log(z+1) − log(z) ≈ 1/z`, which compresses at high charge. What *is* constant is that the template shape is independent of `M`, so it can be precomputed once and shifted across mass candidates.

- **Envelope detection is a binned sparse template match, not FFT / not autocorrelation.** `bin_mul_factor = 2.5 / (tol_ppm × 1e-6)`, bins are ~tol/2.5 ppm wide, the template is precomputed as integer bin offsets, and the sweep is a sparse lookup over (candidate mass bin × charge). Confirmed from `SpectralDeconvolution.cpp` (`universal_pattern_.push_back(-log(i + 1))`, `getBinNumber_`, `harmonic_pattern_matrix_`).

- **Harmonic suppression** is done against integer multiples/divisors of the candidate charge at orders `{2, 3, 5, 7, 11}`, hardcoded.

- **Averagine scoring** is cosine similarity via `getCosine` against a precomputed averagine isotope pattern, with a search over monoisotopic-offset candidates to tolerate non-monoisotopic brightest peaks.

- **The log-transform + template match transfers cleanly to a sparse cross-run feature list.** It is position-agnostic and works on a list of feature m/z values just as well as on a dense spectrum. The averagine deisotoping step does *not* transfer; that belongs at feature extraction upstream (inside the flood-fill in Step 2), not at charge grouping.

## Important caveat

The Cell Systems paper (Jeong et al. 2020) was paywalled during extraction — abstract-only on the publisher page, 403 on the PDF, no PMC entry on the DOI. The math above is **reconstructed from the OpenMS reference implementation** (`SpectralDeconvolution.cpp`, `FLASHHelperClasses.cpp`) rather than quoted directly from the paper. If a discrepancy turns up later between the reference implementation and the published method, the implementation should be treated as authoritative for our purposes.

## Applicability to the top-down engine

- **Use in Step 2c.** The primitive is the template match over log-transformed feature m/z values to group co-eluting cross-run features by implied monoisotopic mass.
- **Input: the list of `FeatureBox` m/z centers from Step 2** (grouped by RT co-elution cluster) — not a dense spectrum, not a per-scan peak list.
- **Output: candidate `(M, {z_i})` assignments** for each RT cluster, ready for Step 3 consensus deconvolution / DB matching.
- **Do not port** the averagine cosine-similarity scoring into Step 2c. Isotope-envelope fit belongs either inside the Step 2 flood-fill or in Step 3 deconvolution, not at charge grouping.
- **Do port** the harmonic-suppression idea — otherwise z=2 envelopes will masquerade as z=1, etc.
