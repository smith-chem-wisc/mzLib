# MslDbBuilder

Standalone .NET 8 console tool for the MetaMorpheus **ManySearchTask** spectral-library (`.msl`)
workstream: converts protein FASTA databases into mzLib `.msl` spectral libraries (the peptide
source / fragmentation cache the parallel search reads instead of digesting FASTAs at search time).

## Dependency

References the **merged mzLib** (PR #1036, branch `mzlib_speclib`) directly via `ProjectReference`
to a local checkout at `E:\GitClones\mzLib`. Edit `MslDbBuilder.csproj` if your mzLib path differs.
It is built against the local `mzlib_speclib` API (`MslWriter`, `MslLibrary`, `RetentionTimePredictorFactory`).

## Build

```
dotnet build MslDbBuilder.csproj -c Release
```

## Modes

- `MslDbBuilder [--verify] [--lean] <outDir> <fasta...>` — FASTA → one `.msl` per database
  (trypsin, **0 missed cleavages**, min length 7, var init-Met, Carbamidomethyl on C / Oxidation on M,
  max 2 mods; target+decoy Reverse; target-wins dedup by full sequence). `--lean` stores no fragments
  (precursor + iRT + sequence only); `--verify` round-trips and reports fragment-count mismatches.
  Peptides with an undefined mass (ambiguous residue `X` → NaN m/z) are skipped.
- `--read <msl...>` — reconstruct peptides and re-fragment in double precision; report parse/mismatch counts.
- `--merge <outDir/out.msl> <msl...>` — combine many `.msl` into one merged file, stamping each entry's
  accession as `db|accession` (the search re-splits per database).
- `--shardbuild <outDir> <numShards> <fasta...|@listfile>` — build N LEAN merged shards directly from
  FASTAs, streaming (`MslWriter.WriteStreaming`), size-balanced, each entry tagged `db|accession` with
  Chronologer iRT. NOTE: per-shard parallel Chronologer construction currently contends on the TorchSharp
  weights temp file (planned Phase1 shared-predictor / Phase2 parallel-writer redesign).
- `--rtcalib <BasePSMs.psmtsv> [ssrcalc]` — regress observed RT vs predicted iRT; report residual (RT window).
- `--checkrt <msl...>` — confirm stored iRT landed in the index.
- `--streamtest <N> [out]` — write N synthetic lean entries via `WriteStreaming` and round-trip read (writer test).
- `--probe <msl>` — open a file with both `LoadIndexOnly` and `Load`; report entry counts and full exceptions.

Companion handoff docs live alongside the wider workstream (`deliverables/RESUME*.md`) and are not part of this repo.
