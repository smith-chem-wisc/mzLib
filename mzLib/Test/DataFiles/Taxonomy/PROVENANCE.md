# Test taxonomy fixtures

Three recorded UniProt responses, so the `UniProtTaxonomyClient` tests run offline. They are **test fixtures, not
production data**: nothing in mzLib's shipping code reads them. The taxonomy mzLib ships is
`UsefulProteomicsDatabases/Resources/taxonomy_nodes.tsv.gz`; see the `PROVENANCE.md` beside it.

## Sources

UniProt data is CC BY 4.0 (https://www.uniprot.org/help/license). All three were fetched on 2026-11-12 from UniProt
release 2026_04 (`x-uniprot-release` response header).

| File | Request | Holds |
|---|---|---|
| `uniprot-proteomes-page1.json` | `https://rest.uniprot.org/proteomes/search?query=upid:UP000000625 OR upid:UP000005640 OR upid:UP000464024&size=5` | UP000000625 (E. coli K12) and UP000005640 (human) |
| `uniprot-proteomes-page2.json` | the same request | UP000464024 (SARS-CoV-2) |
| `uniprot-taxa.json` | `https://rest.uniprot.org/taxonomy/taxonIds/9606,83333,2697049?fields=id,rank,parent,scientific_name` | the three organisms' own ranks and parents |

## How they were trimmed

- The one proteome response was split into two files so the tests can serve them as consecutive pages joined by a
  `Link: <...>; rel="next"` header, which is how UniProt pages a longer listing.
- Each proteome keeps only the five properties the client reads: `id`, `taxonomy`, `proteomeType`, `proteinCount` and
  `taxonLineage`. Every value kept is verbatim, including each lineage entry's `rank` and `hidden` flag.
- `uniprot-taxa.json` is the response body unchanged apart from indentation.

## Why these three

| Proteome | What it exercises |
|---|---|
| UP000000625, E. coli K12 | an organism below species (a strain): its species is an ancestor, and UniProt marks that ancestor hidden |
| UP000005640, human | an organism that is itself the species: the proteome listing gives it no rank, so the species lookup only works once the taxonomy endpoint has supplied one |
| UP000464024, SARS-CoV-2 | a virus, under a different root from the other two; its own rank is "no rank" and its species is an ancestor |

## Refreshing

Re-run the two requests, trim as above, and **run the tests**. `UniProtTaxonomyClientTests` pins names and ranks from
these files on purpose: if UniProt renames a taxon or moves one, that should surface as a failing test and a
deliberate update, not pass unnoticed.
