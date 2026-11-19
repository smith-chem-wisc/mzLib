# SDRF drafter fixtures: provenance

These are PRIDE deposits that nobody annotated, the kind of deposit `SdrfDrafter` exists for. The
curated corpus the drafter was benchmarked on is the annotated minority, so it can't cover them.

- **`raw_file_lists/`, `PXD049018.sdrf.tsv`, `PXD067622.sdrf.tsv`, `AGING_README.md`:**
  aging, 2026-09-23, file lists from PRIDE's REST listing, designs read by aging from the PRIDE
  record and file names; no publication consulted.
  They are copied byte for byte from the aging project, `design/fixtures/2026-09-23_sdrf-design/`
  at `74cd236`, where `make_fixtures.py` regenerates them. `AGING_README.md` is that folder's
  README, renamed. It gives aging's reading of all twelve designs, and those readings are what the
  tests pin.
- **`records/<PXD>.project.json`:** the PRIDE project record for each pinned deposit. It was
  fetched 2026-09-25 from PRIDE Archive REST v3 through mzLib's own
  `PrideArchiveClient.TryGetProjectAsync`, and serialized with Newtonsoft.Json, so the test reads
  exactly what the drafter receives in production. PRIDE records change, so a re-fetch can differ.

Pinned in `TestSdrfDraftFixtures` so far: PXD049018, PXD067622, PXD058611 (aging's SDRF-A12). The
other nine (SDRF-A13) are here for the next round, and they have no record yet. The hard one is
PXD032202.
