# Sequence-Variant Corpus

A data-driven, executable specification for how mzLib should read, expand, and digest proteins that carry
**sequence variants** and **modifications** together. Each `CorpusCase` record is one test case ("node");
`VariantCorpusTests` runs every case through the real pipeline and compares the result to a hand-authored
expected output that lives on the record itself.

## Why it exists
Variant handling is a combinatorial minefield — a variant can add/remove a modification's target residue,
shift every downstream coordinate (indels), create or destroy a protease cut site, and interact with
post-translational processing. This corpus pins the intended behavior one small, human-checkable case at a
time, so regressions are caught and the intended contract is documented executably.

## How it works
The cases live inline as `CorpusCase` records in `VariantCorpusTests.GetCases()` (a `yield return` list fed to
the tests via `[TestCaseSource]`), matching the house style in `Transcriptomics/TestDigestion.cs`. For each
case the harness:
1. builds a UniProt-style XML entry from the case's `Base` sequence, `Mods`, and `Variants`;
2. calls `ProteinDbLoader.LoadProteinXML` (which **expands** the variants into consensus + variant proteins);
3. digests each expanded protein with the case's `Protease` (top-down for proteoforms);
4. collects each product's `FullSequence`;
5. asserts the produced set matches `ExpectedForms` (count + membership) and is **identical across two
   runs** (determinism).

`Corpus_Node_RoundTrip` runs the same pipeline but first serializes the consensus protein through the real
`ProteinDbWriter` and re-reads it, so encode+decode is exercised too — without the writer ever sitting
upstream of the ground truth (the input XML stays hand-authored).

## The `CorpusCase` record
One node per record; input and expected answer together, compiler-checked, no external data files. Fields:
`Id, Layer, Tests, Base, Mods, Variants, Protease, MaxIsoforms, MaxMods, ExpectedCount, Verdict, Reason,
ExpectedForms, Opts`.
- `Mods`: `Name@Pos` (1-based, consensus frame), `;`-separated, `-` = none.
- `Variants`: `OP=orig VAR=var POS=pos SRC=uniprot|vcf`, `;`-separated, `-` = none.
- `ExpectedForms`: the expected forms, **in canonical order** ("least-modified first"), in mzLib's real
  full-sequence notation, e.g. `PEPT[Biological:Phosphorylation on T]IDE`.
- `ExpectedCount`: kept alongside `ExpectedForms` as an independent over/under-generation guard (asserted to
  equal `ExpectedForms.Length`) — mirrors `RnaDigestionTestCase`, which keeps both `DigestionProductCount` and
  `Sequences`.
- `Opts`: optional extensible `key=value;…` tail for secondary knobs (defaults to `-`); lets the schema grow
  without disturbing existing cases.

Why records, not an external `cases.tsv` + `expected/*.txt`: type safety and IDE navigation, no csproj
copy-to-output plumbing or runtime file-not-found fragility, and input next to answer. A later layer that
genuinely explodes to hundreds of forms (L5 throttle) can assert `ExpectedCount` + a representative subset, or
reintroduce an external payload for that layer only.

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
Add a `yield return new CorpusCase(...)` to `GetCases()` with the intended `ExpectedForms` (canonical order)
and a matching `ExpectedCount`. New knobs go in `Opts` (with a default) so existing cases are untouched. New
mods go in the `ModRegistry` in `VariantCorpusTests`.
