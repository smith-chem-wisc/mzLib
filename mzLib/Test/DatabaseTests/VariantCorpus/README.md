# Sequence-Variant Corpus

A data-driven, executable specification for how mzLib should read, expand, and digest proteins that carry
**sequence variants** and **modifications** together. Each row is one test case ("node"); `VariantCorpusTests`
runs every row through the real pipeline and compares the result to a hand-authored expected output.

## Why it exists
Variant handling is a combinatorial minefield — a variant can add/remove a modification's target residue,
shift every downstream coordinate (indels), create or destroy a protease cut site, and interact with
post-translational processing. This corpus pins the intended behavior one small, human-checkable case at a
time, so regressions are caught and the intended contract is documented executably.

## How it works
For each row of `cases.tsv`, the harness:
1. builds a UniProt-style XML entry from the row's `base` sequence, `mods`, and `variants`;
2. calls `ProteinDbLoader.LoadProteinXML` (which **expands** the variants into consensus + variant proteins);
3. digests each expanded protein with the row's `protease` (top-down for proteoforms);
4. collects each product's `FullSequence`;
5. asserts the produced set matches `expected/<id>.txt` (count + membership) and is **identical across two
   runs** (determinism).

## Files
- **`cases.tsv`** — one row per case. Columns: `id, layer, tests, base, mods, variants, protease,
  max_isoforms, max_mods, opts, expected_count, verdict, reason`.
  - `mods`: `Name@Pos` (1-based, consensus frame), `;`-separated, `-` = none.
  - `variants`: `OP=orig VAR=var POS=pos SRC=uniprot|vcf`, `;`-separated, `-` = none.
  - `opts`: extensible `key=value;…` tail for secondary knobs (`-` = defaults); lets the schema grow without
    disturbing existing rows.
- **`expected/<id>.txt`** — the expected forms for a case, **one per line, in canonical order**
  ("least-modified first"), in mzLib's real full-sequence notation, e.g. `PEPT[Biological:Phosphorylation on T]IDE`.

## Discipline
- **Expected outputs are the oracle** — authored from the intended behavior, **never** captured from current
  code output (that would bake in bugs).
- **Foundation first, then layers.** L0 (this set) covers substitutions and PTMs under top-down digestion with
  non-binding caps. Later layers add indels, post-translational processing, bottom-up digestion, multi-variant
  combinatorics, VCF depth/zygosity, and read/write round-trips.
- **Canonical test mods** are pinned so expected strings are exact (not database-dependent):
  | key | ModificationType | OriginalId | mass | renders as |
  |---|---|---|---|---|
  | `Phosphorylation` | `Biological` | `Phosphorylation` | 79.966331 | `[Biological:Phosphorylation on T]` |

## Adding a case
Add a row to `cases.tsv` and an `expected/<id>.txt` with the intended forms (canonical order). New knobs go in
`opts` (with a default) so existing rows are untouched. New mods go in the registry in `VariantCorpusTests`.
