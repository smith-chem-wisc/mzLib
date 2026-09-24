# SDRF design fixtures

Hand-written SDRFs from the `aging` project (`design/fixtures/2026-09-23_sdrf-design/`, commit
`74cd236`, 2026-09-23). One row per searched raw file. Nothing here was deposited by the authors;
where the PRIDE record does not say something, the fixture says `not available`.

- `PXD067622.sdrf.tsv`: 24 files. Construct (WT/CA) x treatment (DMSO/FA/NDC/THY) x 3, one fraction
  each. The file names number replicates 1..24 across the whole study; the SDRF numbers them 1..3
  within each condition.
- `PXD049018.sdrf.tsv`: 20 files. Two streptavidin pulldowns x 10 gel bands, no replication. The
  record says one is +cGAMP and one is -cGAMP but not which, so `factor value[treatment]` is
  `not available` on both, and both carry biological replicate 1.
