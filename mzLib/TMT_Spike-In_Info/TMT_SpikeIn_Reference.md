# TMT Spike-In Dataset Reference for Quantification Project

## Source

**Paper:** Huang et al., 2020, *Mol Cell Proteomics* 19(10), 1706-1723
**Title:** MSstatsTMT: Statistical Detection of Differentially Abundant Proteins in Experiments with Isobaric Labeling and Multiple Mixtures
**DOI:** 10.1074/mcp.RA120.002105
**Dataset:** SpikeIn-5mix-MS3

---

## Experimental Design Overview

### Sample Composition
- **Background:** 50 µg SILAC HeLa peptides (constant across all channels)
- **Spike-in:** 48 UPS1 (Universal Proteomics Standard) proteins at varying concentrations
- **UPS1 amounts:** 500, 333, 250, and 62.5 fmol, combined with 50 µg SILAC HeLa in duplicate
- **Dilution series:** Produces relative concentrations of 1, 0.667, 0.5, and 0.125 (relative to 500 fmol)
- **Reference sample:** Pool of all four diluted UPS1 peptide samples (286.5 fmol), combined with 50 µg SILAC HeLa in duplicate

### TMT Labeling
- TMT10-plex reagents (channels: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131)
- 10 samples per mixture, labeled and analyzed by LC-MS/MS
- Data acquired using SPS-MS3 (Synchronous Precursor Selection) on Orbitrap Fusion Lumos Tribrid

### Mixtures and Runs
- **5 controlled mixtures**, each profiled in **3 technical replicate MS runs**
- **15 total MS runs** from 5 TMT mixtures
- **Note:** Mixture 1 and Mixture 5 have the exact same design (same channel-to-condition mapping)

---

## Channel-to-Condition Mapping (Supplementary Figure S1.2)

Each value represents the **relative concentration of UPS1 proteins** in that channel. "Ref" = reference channel (pool of all dilutions, used for normalization).

| TMT10plex reagent | 126 | 127N  | 127C  | 128N  | 128C | 129N  | 129C | 130N | 130C  | 131 |
|-------------------|-----|-------|-------|-------|------|-------|------|------|-------|-----|
| **Mixture 1** Run 1-3 | Ref | 0.667 | 0.125 | 0.5   | 1    | 0.125 | 0.5  | 1    | 0.667 | Ref |
| **Mixture 2** Run 4-6 | Ref | 0.5   | 1     | 0.667 | 0.125| 1     | 0.667| 0.125| 0.5   | Ref |
| **Mixture 3** Run 7-9 | Ref | 0.125 | 0.667 | 1     | 0.5  | 0.5   | 0.125| 0.667| 1     | Ref |
| **Mixture 4** Run 10-12| Ref | 1     | 0.5   | 0.125 | 0.667| 0.667 | 1    | 0.5  | 0.125 | Ref |
| **Mixture 5** Run 13-15| Ref | 0.667 | 0.125 | 0.5   | 1    | 0.125 | 0.5  | 1    | 0.667 | Ref |

**Key observations:**
- Channels **126** and **131** are **always reference channels** across all mixtures
- The 8 endogenous channels (127N-130C) contain different UPS1 concentrations
- Conditions are shuffled across channels between mixtures to avoid systematic bias
- Each condition (0.125, 0.5, 0.667, 1) appears **twice per mixture** (duplicate biological replicates)

---

## Raw File Names and Mapping

### File Naming Convention
`161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture{M}_{RR}`

Where `M` = mixture number (1-5), `RR` = technical replicate (01, 02, 03)

### Complete File List (14 files available; Mixture4_01 is missing)

| File Name | Mixture | Tech Rep | Run # |
|-----------|---------|----------|-------|
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01 | 1 | 01 | 1 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02 | 1 | 02 | 2 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03 | 1 | 03 | 3 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_01 | 2 | 01 | 4 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_02 | 2 | 02 | 5 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_03 | 2 | 03 | 6 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_01 | 3 | 01 | 7 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_02 | 3 | 02 | 8 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_03 | 3 | 03 | 9 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_02 | 4 | 02 | 11 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_03 | 4 | 03 | 12 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_01 | 5 | 01 | 13 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_02 | 5 | 02 | 14 |
| 161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_03 | 5 | 03 | 15 |

**Note:** Mixture4_01 (Run 10) is absent from the search results.

---

## Expected Ground Truth

### UPS1 Proteins
- 48 total UPS1 proteins in the standard
- After removing ambiguous protein groups (shared spiked-in + background): **40 UPS1 proteins remaining**
- After removing endogenous "light" species from SILAC-HeLa: **21 truly differentially abundant UPS1 proteins** for SpikeIn-5mix-MS3
- 20 for SpikeIn-5mix-MS2

### Expected Fold Changes (UPS1 proteins)
Pairwise comparisons use ratios of UPS1 concentrations. The conditions are labeled by their UPS1 concentration relative to the highest (500 fmol = 1.0):

| Comparison | Ratio | True Fold Change |
|------------|-------|------------------|
| 1 vs 0.667 | 1/0.667 | 1.5 |
| 1 vs 0.5 | 1/0.5 | 2.0 |
| 1 vs 0.125 | 1/0.125 | 8.0 |
| 0.667 vs 0.5 | 0.667/0.5 | 1.33 |
| 0.667 vs 0.125 | 0.667/0.125 | 5.328 |
| 0.5 vs 0.125 | 0.5/0.125 | 4.0 |

### HeLa Background Proteins
- Should show **no differential abundance** between conditions (constant across all channels)
- Any detected fold change in HeLa proteins represents **false positives**
- Serve as negative controls for specificity assessment

---

## Data Format (AllPSMs.psmtsv)

Tab-delimited file with 67 columns per row.

### Key Columns
| Column # | Name | Description |
|----------|------|-------------|
| 1 | File Name | Raw file name (maps to mixture/run) |
| 2 | Scan Number | MS scan identifier |
| 14 | Base Sequence | Unmodified peptide sequence |
| 15 | Full Sequence | Modified peptide sequence |
| 27 | Accession | Protein accession (UniProt ID) |
| 28 | Name | Protein name |
| 29 | Gene Name | Gene symbol |
| 34 | Decoy | Whether this is a decoy match |
| 40 | Decoy/Contaminant/Target | Classification |
| 56 | PEP | Posterior Error Probability |
| 57 | PEP_QValue | PEP-based q-value |
| **58** | **126** | TMT channel 126 intensity |
| **59** | **127N** | TMT channel 127N intensity |
| **60** | **127C** | TMT channel 127C intensity |
| **61** | **128N** | TMT channel 128N intensity |
| **62** | **128C** | TMT channel 128C intensity |
| **63** | **129N** | TMT channel 129N intensity |
| **64** | **129C** | TMT channel 129C intensity |
| **65** | **130N** | TMT channel 130N intensity |
| **66** | **130C** | TMT channel 130C intensity |
| **67** | **131** | TMT channel 131 intensity |

**Note:** The current search results have all-zero reporter ion intensities. This will need to be resolved (likely a search parameter issue) before the testing harness can use real data.

### Available Search Results

**HeLa+SILAC Search** (`TMT_HELA_SILAC_Search/`):
- Searched against human database
- 115,390 target PSMs (q <= 0.01)
- 3,472 protein groups (1% FDR)
- ~6,000-8,700 PSMs per file

**UPS Search** (`UPS_Search/`):
- Searched against UPS database only
- 2,141 target PSMs (q <= 0.01)
- 47 protein groups (1% FDR)
- ~100-180 PSMs per file

---

## Existing mzLib Quantification Framework

### Project Structure
```
Quantification/
├── Interfaces/
│   ├── ICollapseStrategy.cs      -- Combine sample columns (e.g., tech replicates)
│   ├── INormalizationStrategy.cs  -- Normalize intensity values
│   └── IRollUpStrategy.cs         -- Aggregate rows (PSM→Peptide→Protein)
├── Strategies/
│   ├── Collapse/
│   │   ├── NoCollapse.cs          -- Identity (no-op)
│   │   └── SumCollapse.cs         -- Sum intensities for grouped columns
│   ├── Normalization/
│   │   └── NoNormalization.cs     -- Identity (no-op)
│   └── RollUp/
│       └── SumRollUp.cs           -- Sum rows for each higher-level entity
├── QuantMatrix.cs                 -- Generic matrix wrapper (MathNet DenseMatrix)
├── QuantificationEngine.cs        -- Main processing pipeline
├── QuantificationParameters.cs    -- Strategy configuration container
├── QuantificationResults.cs       -- Output object (mostly stubbed)
└── QuantificationWriter.cs        -- File output (all methods throw NotImplementedException)
```

### Processing Pipeline (QuantificationEngine.Run())
1. **Pivot** - Convert raw PSM list to QuantMatrix (long → wide format)
2. **Normalize PSM Matrix** - Apply SpectralMatchNormalizationStrategy
3. **Roll Up to Peptides** - PSM-to-Peptide mapping + SpectralMatchToPeptideRollUpStrategy
4. **Normalize Peptide Matrix** - Apply PeptideNormalizationStrategy
5. **Collapse Samples** - Combine fractions/tech replicates via CollapseStrategy
6. **Roll Up to Proteins** - Peptide-to-Protein mapping + PeptideToProteinRollUpStrategy
7. **Normalize Protein Matrix** - Apply ProteinNormalizationStrategy

### Implementation Status

| Component | Status | Notes |
|-----------|--------|-------|
| Core Interfaces | Complete | Strategy pattern for all 3 operations |
| QuantMatrix | Complete | Generic matrix with row/column keys |
| QuantificationEngine | Mostly Complete | `GetAllPeptideToProteinMap()` stubbed |
| SumRollUp | Complete | Sums rows for roll-up |
| SumCollapse | Complete | Groups and sums columns |
| NoNormalization | Complete | Identity/no-op |
| NoCollapse | Complete | Identity/no-op |
| QuantificationWriter | Stubbed | All methods throw NotImplementedException |
| QuantificationResults | Minimal | Only `Summary` string property |

### Critical Redesign: Per-File Matrix Construction for TMT Data

The current `Pivot()` method is inefficient for TMT data. It combines all PSMs from every file into one giant matrix with one column per sample across all files. Since each PSM is only observed in one file (with intensities in that file's ~10 channels), this creates a massively sparse matrix where most values are zero/NA.

**Revised pipeline for TMT data:**
1. **Group PSMs by file** - Partition input PSMs into per-file groups
2. **Build per-file matrices** - Each matrix is (PSMs_in_file x channels_in_file), fully dense
3. **File-level normalization** - Normalize within each file (channels are directly comparable within a file)
4. **Roll up to peptides within each file** - Aggregate PSMs to peptide-level intensities per file
5. **Combine per-file peptide matrices** - Merge into a single peptide matrix (peptides x all_samples)
6. **Continue pipeline** - Collapse, protein roll-up, protein-level normalization

This aligns with the MSstatsTMT approach where spectrum-level normalization is applied within runs before cross-run protein summarization. It also means normalization operates on dense data where all values are meaningful.

### What Needs to Be Built/Improved
1. **Normalization strategies** - Currently only NoNormalization exists. Need:
   - Global median normalization (MSstatsTMT approach)
   - Reference channel normalization (local ratio-based)
   - Potentially: quantile normalization, variance stabilizing normalization
2. **Roll-up strategies** - Currently only SumRollUp. Need:
   - Median roll-up
   - Weighted average / isobar-style roll-up
   - Tukey's median polish (MSstatsTMT approach)
3. **Testing harness** - Use this spike-in dataset to evaluate strategies against known ground truth
4. **QuantificationWriter** - Implement output
5. **QuantificationResults** - Define proper results structure

---

## Quantification Workflow Strategies (from Literature)

### Table I Summary from Paper

| Step | Ratio+Median+Limma | Sum+IRS+edgeR | Proteome Discoverer | MSstatsTMT |
|------|-------------------|---------------|--------------------|-----------|
| **Spectrum-level norm** | Local ratio-based (log2 ratio to ref channel) | None | None | Global median (equalize median log2 across all spectra/channels/runs) |
| **Protein summarization** | Median of log2 ratios | Sum of all intensities (not log) | Sum of all intensities (not log) | Tukey's median polish (handles missing values via AFT model) |
| **Protein-level norm** | Global zero median (subtract median log2 ratio) | Global equal sum + Local IRS with ref channel | Global equal sum + Protein scaling to avg 100 | Local normalization with reference channel |

### MSstatsTMT Approach (Most Relevant)
1. **Global median normalization:** log2-transform intensities, equalize median log2 across all spectra, channels, and runs
2. **Protein summarization via Tukey's median polish:** Fits two-way robust additive model, handles missing values with Accelerated Failure Time model
3. **Local protein-level normalization with reference channel:** For each protein, equalize reference channel summaries across runs to their median, then apply corresponding shifts to all channels

---

## Testing Strategy

### Evaluation Metrics (from Paper)
- **True Positives (TP):** UPS1 proteins correctly identified as differentially abundant
- **False Positives (FP):** HeLa proteins incorrectly identified as differentially abundant
- **Sensitivity:** TP / (TP + FN)
- **Specificity:** TN / (TN + FP)
- **Empirical FDR:** FP / (TP + FP)
- **Fold change accuracy:** Compare estimated fold changes to known true fold changes
- **ROC/AUC:** Sensitivity vs 1-specificity at various FDR cutoffs

### What Our Testing Harness Should Measure
1. **For UPS proteins:** How accurately do we recover the known fold changes?
2. **For HeLa proteins:** How close to 1.0 (no change) are the estimated fold changes?
3. **Normalization quality:** Do reference channels equalize properly across runs?
4. **Strategy comparison:** Compare different normalization + roll-up combinations

### Design Note
This is a **controlled mixture without biological variation** (Section 3.6 of supplemental). The appropriate MSstatsTMT model removes the Subject term. The design has:
- M = 5 mixtures
- T = 3 technical replicates per mixture (except Mixture 4 has only 2)
- C = 4 conditions (0.125, 0.5, 0.667, 1.0) + reference
- B = 2 biological replicates per condition per mixture

---

## Future Work: Port FlashLFQ Algorithms

The FlashLFQ project within mzLib contains sophisticated implementations that should eventually be translated into the Quantification framework. **These are longer-term goals; getting the basics working comes first.**

### FlashLFQ IntensityNormalizationEngine (3-stage pipeline)
1. **NormalizeFractions()** - Bounded Nelder-Mead optimization to find per-fraction linear normalization factors minimizing log-intensity error vs. a reference biorep. Factors bounded to [0.3, 3.0]. Uses coarse grid search to seed optimizer.
2. **NormalizeBioreps()** - Median fold-change normalization across biological replicates/conditions. Reference = condition 1, biorep 1. Assumes most peptides don't change between conditions.
3. **NormalizeTechreps()** - Median fold-change normalization between technical replicates within the same condition/biorep/fraction.

### FlashLFQ CalculateProteinResultsMedianPolish + MedianPolish
- Builds a (peptides+1 x samples+1) matrix with log2-transformed intensities (NaN for missing)
- Prunes peptides with only a single valid measurement when others have multiple
- **MedianPolish** is actually a weighted mean polish (robust to missing values):
  - Weights = inverse-squared deviation from median (down-weights outliers)
  - Iterative row/column sweeps extracting effects until convergence
  - Result: overall effect + row effects (peptide ionization efficiency) + column effects (sample abundance shift)
- Protein intensity per sample = 2^(column_effect) * reference_intensity

### Relevance to TMT Quantification
- The median polish algorithm maps directly to protein summarization (Tukey's median polish from MSstatsTMT)
- The normalization stages could inform TMT normalization strategies
- The weighting scheme for handling missing values is particularly relevant for TMT data where missing intensities are common

---

## Critical Issue: Zero Reporter Ion Intensities

The current search results have all-zero TMT reporter ion intensities. Before the testing harness can function with real data, one of the following must happen:
1. Re-search the data with correct parameters to extract reporter ion intensities
2. Use pre-existing reporter ion data from another source
3. Create synthetic test data based on the known experimental design for initial development

For the Wiggum Loop development, option 3 (synthetic data with known properties) may be the most practical starting point, with real data integration following once the intensity issue is resolved.
