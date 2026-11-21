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

# Vendored taxonomy snapshot

`taxonomy_nodes.tsv.gz` and `proteome_taxa.tsv.gz` are a pinned snapshot of the taxonomy above UniProt's
reference proteomes, embedded into `UsefulProteomicsDatabases.dll` and read by `TaxonomyTree.LoadPinned()` and
`ProteomeTaxonMap.LoadPinned()`.

It is **not** the NCBI taxonomy ontology, which `ControlledVocabulary` deliberately leaves out at roughly 600 MB.
It holds only the organisms of the reference proteomes and the ancestors of each: enough to roll any of those
organisms up to its family, order or phylum, in under 2 MB. An organism it lacks is added with
`TaxonomyTree.WithOverlay`; a fuller tree is built with `UniProtTaxonomyClient`.

It is pinned for the same reason as the vocabularies above: so that the same search annotated a year apart places
an organism in the same family. NCBI and UniProt do move and rename taxa between releases.

| File | Source | Version | Fetched | Size | sha256 |
|---|---|---|---|---:|---|
| `taxonomy_nodes.tsv.gz` | https://rest.uniprot.org/proteomes/search?query=reference:true (each proteome's `taxonLineage`), organisms' own ranks from https://rest.uniprot.org/taxonomy/taxonIds | UniProt release `2026_04` | 2026-11-12 | 1,077,630 B | `1333560396a8855a0d340bf172c6096d2b5b393ac159b049bcd966cf9c91dc88` |
| `proteome_taxa.tsv.gz` | the same request | UniProt release `2026_04` | 2026-11-12 | 768,186 B | `9261aebfc8e5d674f91690e0e011198d529f009cda63be563e77e7f8feec17a6` |

Contents: 37,623 reference proteomes and 68,941 taxa under two roots, `cellular organisms` and `Viruses`. The
source, release and query are also in each file's own `#!` header lines.

22 of the 37,623 proteomes are listed by UniProt under a taxon it names by number only, with no name and no
lineage (taxa NCBI has since merged or withdrawn; UP000031433 is one). Each is in `proteome_taxa.tsv.gz` with
that number, and that taxon is deliberately absent from `taxonomy_nodes.tsv.gz`: asking the tree about it answers
"not in the tree" rather than a made-up place.

## Refreshing

Both files are written by mzLib's own client, so a refresh is a few lines against the build being released:

```csharp
using var client = new UniProtTaxonomyClient();
UniProtTaxonomySnapshot snapshot = await client.FetchAsync(UniProtTaxonomyClient.ReferenceProteomesQuery);
TaxonomySnapshotTsv.WriteNodes("taxonomy_nodes.tsv.gz", snapshot.Tree);
TaxonomySnapshotTsv.WriteProteomes("proteome_taxa.tsv.gz", snapshot.Proteomes);
```

It makes about 270 requests and takes a little over two minutes. Refresh both files together, from one run: the
proteome table names taxa that must be in the node table of the same release.

Then update the table above and **run the tests**. `TestPinnedTaxonomy` pins the lineages of human, mouse, yeast,
E. coli K12 and SARS-CoV-2 on purpose: if a refresh moves one of them to a different family, that must surface as a
failing test and a deliberate decision, not as a silent shift in every result grouped by taxon afterwards.

## Licence

UniProt data is licensed **CC BY 4.0** (https://www.uniprot.org/help/license). As for the vocabularies, the
attribution lives in `THIRD-PARTY-NOTICES.md` at the repository root.
