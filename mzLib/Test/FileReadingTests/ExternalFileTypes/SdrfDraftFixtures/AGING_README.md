# Design fixtures for sdrf (SDRF-A12, SDRF-A13) and QuantProject (027 §1)

Written 2026-09-23 by aging. Regenerate with `make_fixtures.py` (same folder).

## What is here

- `raw_file_lists/<PXD>_raw_files.tsv`: for each of the twelve datasets in our catalog, **every raw
  file PRIDE lists** (the REST listing, category RAW, which matched the FTP listing in every case),
  with PRIDE's size, our sha256 and whether we searched it. All twelve are whole deposits: every
  listed file was searched (the subset searches from before D38 were replaced).
- `PXD067622.sdrf.tsv`, `PXD049018.sdrf.tsv`: hand-written SDRFs, our reading of the design, one row
  per searched raw file. The treatment cells carry the D31 override columns
  (`comment[treatment source]`, `... source reference`, `... source curator`).

Nothing here was deposited by the authors. Where the record does not say something, the fixture says
`not available` rather than guessing.

## Our reading of each design

The first three are the ones asked for (SDRF-A12). The rest are SDRF-A13's extra fixtures; our
readings of them are shorter because we have only the PRIDE record and the file names.

| PXD | files | design, as we read it | what the file names give |
|---|---:|---|---|
| **PXD067622** | 24 | HeLa, SPRTN-TurboID proximity labelling with streptavidin pulldown. **Construct (WT or CA) x treatment (DMSO, formaldehyde, nocodazole, thymidine) x 3.** The protocol names FA 1 mM 1 h, nocodazole 100 ng/ml 16 h (M phase) and thymidine 2 mM 16 h (S phase); DMSO is read as the vehicle control from the file names. The record never defines `CA`: it is a second SPRTN-TurboID construct, not a treatment. We read the three as biological replicates; the record does not say. | Everything. `WT`/`CA` and `DMSO`/`FA`/`NDC`/`THY` are tokens, and **the trailing index runs 1..24 across the whole study** (WT_DMSO1-3, CA_DMSO4-6, WT_NDC7-9, CA_NDC10-12, WT_THY13-15, CA_THY16-18, WT_FA19-21, CA_FA22-24), not 1..3 per condition. The SDRF numbers biological replicates 1..3 within each condition. This is QuantProject's §3.2 renumbering case. |
| **PXD049018** | 20 | BJ cells, TurboID-NES labels cytoplasmic proteins, then **with or without cGAMP**; nuclei isolated, biotinylated proteins pulled down with streptavidin, run on a gel. **2 samples x 10 bands, no replication.** | Sample (`MSB67868A`, `MSB67869A`) and fraction (`Band_01..10`). **Which sample got cGAMP is not stated anywhere in the record**, so `factor value[treatment]` is `not available`. |
| **PXD058611** | 36 | Mouse liver, H2S donor study. The protocol describes a **global proteome aliquot and a persulfidome pulldown** (DCP-Bio1, streptavidin) per sample, plus **negative controls without DCP-Bio1 in duplicate per group**. So the deposit mixes two assay types and controls. | **Nothing.** Files are `178.raw` .. `213.raw`. This is the no-structure fixture: the right drafter answer is "no structure found". |
| PXD023381 | 20 | Streptavidin pulldown on a gel (19.2% of intensity is streptavidin): 2 samples x 10 bands. | Sample (`MSB40881`/`MSB40882`) and fraction (`Band_01..10`). Which sample is which condition: not in the names. |
| PXD028852 | 20 | Same shape as PXD023381, with a GFP-tagged bait (GFP 1.3% of intensity). | Sample (`MSB42871`/`MSB42872`) and fraction (`Band_01..10`). |
| PXD024803 | 28 | CHD6 immunoprecipitation. | **Two acquisition batches** by date and instrument prefix: `20180519_Q2_UC_..._1..6` (6 files) and `20180824_Q1_AWH_..._sample_1..22` (22). The design within each batch is not in the names. |
| PXD032040 | 15 | IP with controls. | Two batches by prefix: `20161028_AWH_..._IgG1-3` and `_ctrl1-3` (IgG IPs and controls), and `Q1_ColID_292_2718_10..18` (9 files, no condition tokens). |
| PXD050351 | 20 | Kinobead (chemical probe) enrichment. | A letter A..H and a replicate 1..3 (`GA1`..`GH2`): A-D have 3 replicates, E-H have 2. `R1`/`R2` looks like a re-injection of GC1 and GC3. What the letters mean is not in the names. |
| PXD031782 | 12 | Mouse, whole proteome. | Sample numbers only (`S9`..`S20`). |
| PXD027318 | 18 | Human mammary epithelial cells (HMEC), whole-cell lysate, **three acquisition batches**. | Batch by date prefix (`11-20-17`, `02-01-18`, `20210517`). Within them, condition tokens (`young`, `sen`, `E-Luc`, `E-hTERT`, `P1..3` vs `Q1..3`, one `P2alt`), a letter per donor or line (`A`, `C`, `D`, `H`), `D1`/`D2`, and a trailing run index that is not a replicate. |
| PXD032202 | 21 | Eye tissue: control vs PXG (and one PXF set). | A sample number, then a condition spelled several ways (`CONT`, `CONTROL_`, `PXG`, `PXG_`, `PXF_`) and a replicate 1..3 that restarts per block: seven blocks of three. Whether the spelling change marks a batch is our guess, not the record's. |
| PXD036557 | 18 | HGPS and control fibroblasts. A deposited SDRF exists but is a skeleton. | `GM1..GM8` (six lines) x `_a/_b/_c` (three runs each), after a `QE-` run number. The HGPS/control labels are **not recoverable** (our G18). |
