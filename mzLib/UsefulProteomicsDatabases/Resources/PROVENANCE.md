# Vendored ontology snapshots

These files are pinned third-party controlled vocabularies, embedded into
`UsefulProteomicsDatabases.dll` and read by `ControlledVocabulary`.

They are deliberately **pinned**, not downloaded: resolving a term against a moving target means
the same experiment annotated a year apart gets different accessions, and a corpus built that way
silently stops joining to itself. Refreshing one is a reviewable change to this repository.

| File | Source | Version | Fetched |
|---|---|---|---|
| `psi-ms.obo` | https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo | `data-version: 4.1.258` | 2026-08-07 |
| `pride_cv.obo` | https://raw.githubusercontent.com/PRIDE-Utilities/pride-ontology/master/pride_cv.obo | `data-version: releases/2026-06-19` | 2026-08-07 |

## Refreshing

```bash
curl -sL -o psi-ms.obo   https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
curl -sL -o pride_cv.obo https://raw.githubusercontent.com/PRIDE-Utilities/pride-ontology/master/pride_cv.obo
```

Then update the versions above and **run the tests**. `TestControlledVocabulary` pins specific
accession-to-name mappings on purpose: if a refresh changes what `Q Exactive` or `HCD` resolves to,
that must surface as a failing test and a deliberate decision, not as a silent shift in every file
written afterwards. `AgreesWithTheMzmlReadersOwnAccessionTable` additionally checks the snapshot
still agrees with `Readers.Mzml.DissociationDictionary`.

## Licence

Both are licensed **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/), which requires
attribution to be conveyed to downstream recipients. Because mzLib ships them inside a NuGet
package, that attribution lives in `THIRD-PARTY-NOTICES.md` at the repository root, which is
included in the package.
