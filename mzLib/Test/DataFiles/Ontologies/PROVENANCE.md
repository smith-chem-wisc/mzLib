# Test ontology fixtures

These three files exist so the **required** CI job can run offline. They are **test fixtures, not
production data** — nothing in mzLib's shipping code reads them. Production still downloads the
real ontologies through `UsefulProteomicsDatabases.Loaders`, and the tests that prove that download
path works are the live canaries described at the bottom of this file.

Before these existed, `Loaders.LoadPsiMod` / `LoadUniprot` / `LoadUnimod` downloaded on first use
(they only fetch when the target file is absent), so a UniProt hiccup reddened the required `build`
job and took ~150 XML-parsing, variant, decoy and digestion tests down with it — tests that have
nothing to do with any external service. Tagging those tests `[Category("ExternalService")]` would
have been the wrong fix: that label marks a canary, and moving them out of the blocking job would
have left a large hole in the only check that gates a merge.

## Why the files are trimmed

Each is cut to what the assertions actually reach, not shipped whole. The full ontologies total
7.0 MB; these total 490 KB.

| File | Upstream size | Here | Kept |
|---|---|---|---|
| `PSI-MOD.obo.xml` | 4.27 MB | 215 KB | The 85 of 2001 terms carrying a `FormalCharge` xref, plus `MOD:00046` |
| `unimod_tables.xml` | 2.51 MB | 37 KB | `Phospho` (record_id 21) and `Oxidation` (record_id 35), verbatim |
| `ptmlist.txt` | 244 KB | 244 KB | Everything — untrimmed |

**`PSI-MOD.obo.xml`** — its only consumer is `Loaders.GetFormalChargesDictionary(obo)`, which selects
terms whose `xref_analog` has `dbname == FormalCharge` and ignores every other term. Dropping the
other 1916 terms is therefore invisible to callers: the trimmed file yields a **byte-identical
85-entry dictionary**, verified against the full file before it was committed. `MOD:00046`
(phosphoserine) is kept although it has no formal charge, because `FilesLoading` asserts that it has
none — that assertion can only mean something if the term is present to be checked.

**`unimod_tables.xml`** — the only offline consumer is `ProFormaConverterRealUnimodTests`, which
resolves exactly two modifications. Both records are kept **verbatim**, so the test still validates
against real Unimod definitions rather than hand-built ones, which is its stated purpose. Phospho
keeps all 9 of its specificities, so the `UNIMOD:21` → S/T/Y motif ambiguity the test exists to
exercise is fully intact.

**`ptmlist.txt`** is deliberately **not** trimmed. It feeds every fixture setup, and the test protein
XML databases reference modifications *by name*. A name dropped from this file would not fail —
`ProteinDbLoader` would quietly parse one modification fewer and some unrelated count assertion would
shift, which is far worse than a loud failure. It is also the small one, so trimming is the
high-risk, low-reward case.

## Filenames

These deliberately do **not** reuse the `*2` names (`PSI-MOD.obo2.xml`, `unimod_tables2.xml`,
`ptmlist2.txt`). Those names belong to the live canaries, which depend on the file being *absent* so
that `Loaders.Load*` downloads it. A trimmed `unimod_tables2.xml` sitting in the output directory
would silently satisfy `FilesLoading`'s `Count > 2700` check against a two-modification file.

## Sources and refresh

| File | Source | Version stamp | Fetched |
|---|---|---|---|
| `ptmlist.txt` | `http://uniprot.org/docs/ptmlist.txt` | Release 2026_02 of 10-Jun-2026 | 2026-08-09 |
| `unimod_tables.xml` | `http://www.unimod.org/xml/unimod.xml` | latest `date_time_modified` 2018-08-13 | 2026-08-09 |
| `PSI-MOD.obo.xml` | `https://github.com/smith-chem-wisc/psi-mod-CV/blob/master/PSI-MOD.obo.xml?raw=true` | header date 30:05:2014 | 2026-08-09 |

To refresh, re-fetch the three URLs and re-apply the filters above:

- `PSI-MOD.obo.xml` — keep the header, then every `<term>` whose `xref_analog` contains
  `<dbname>FormalCharge</dbname>`, plus `MOD:00046`. Confirm the resulting formal-charges dictionary
  matches the one built from the full file before committing.
- `unimod_tables.xml` — keep the `<umod:elements>`, `<umod:amino_acids>` and `<umod:mod_bricks>`
  blocks unchanged (they are needed for mass computation), and inside `<umod:modifications>` keep only
  the `Phospho` and `Oxidation` records.
- `ptmlist.txt` — copy as-is.

A refresh is a reviewable change to this repository rather than something that happens to CI
overnight, which is the same reasoning applied to the pinned vocabularies in `UsefulProteomicsDatabases`.

## The canaries these do not replace

Live download coverage still exists, tagged `[Category("ExternalService")]` so it runs in the
dedicated non-blocking job, and routed through `ExternalServiceTestHelper` so a third-party outage
reports **Skipped** with a clear message while a genuine contract break still **Fails**:

`TestUpdateUnimod`, `TestUpdatePsiMod`, `TestUpdatePsiModObo`, `TestUpdateUniprot`, `FilesEqualHash`,
`FilesLoading`, `TestRetrieveUniProtProteome`, `TestDownloadAvailableUniProtProteomes`.
