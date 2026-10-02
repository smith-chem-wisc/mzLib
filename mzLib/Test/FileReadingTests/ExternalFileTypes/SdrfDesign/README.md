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

Isobaric subsets, cut from `bigbio/sdrf-annotated-datasets` (Apache-2.0) at `4f823dc` (2026-08-06),
every column and cell as deposited; only rows were removed, whole files at a time.

- `PXD008841.subset.sdrf.tsv`: 40 rows, 4 files, TMT10. Fractions 1 and 2 of plex 1 and plex 5 of a
  9-plex bridge design (360 files in all). No column names the plex; it survives only in a file-name
  token spelled two ways (`TMTpool1`, `TMT_pool5`). Channel 131N is the reference, `source name`
  `pool`, in every plex.
- `PXD061609.subset.sdrf.tsv`: 108 rows, 6 files, TMTpro 18. Fractions 1 and 2 of each of three SAX
  arms (`NoSAX`, `HighSalt`, `LowSalt` in `comment[sample preparation batch]`) of ONE plex: all 18
  samples are in every arm, and each arm numbers its fractions from 1.
