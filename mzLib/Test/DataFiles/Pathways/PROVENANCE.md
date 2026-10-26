# Test pathway fixtures

Two files cut from **Reactome release 96** so the Reactome release tests run offline. They are **test fixtures, not
production data**: nothing in mzLib's shipping code reads them. Callers pass the paths of the real release files to
`ReactomeRelease.Load`.

## Sources

Reactome data is CC0 (https://reactome.org/license, section 1c: "All data in the Reactome database and files derived
from that data are licensed under the Creative Commons Public Domain Dedication (CC0)").

| File | Source | Upstream size | Upstream sha256 | Fetched |
|---|---|---:|---|---|
| `ReactomePathways-96-trimmed.txt` | `https://download.reactome.org/96/ReactomePathways.txt` (23,498 pathways, all species) | 1,583,655 B | `eac3ac074b3849a6d928260f4f71920f346584b0114afa7b60fcc8e40eefd9ec` | 2026-10-25 |
| `ReactomePathwaysRelation-96-trimmed.txt` | `https://download.reactome.org/96/ReactomePathwaysRelation.txt` (23,612 parent-child pairs) | 631,240 B | `72e2da17ef5806d287a0d932cd6d161fea62ef7a567b6e95601b67463c07170e` | 2026-10-25 |

Release 96 is the release UniProt 2026_03's Reactome cross-references carry.

## How they were trimmed

**Kept:** every pathway in the ancestor closure of six seeds. Lines are copied verbatim, in the upstream order, with LF
line endings.

| Seed | Why |
|---|---|
| `R-HSA-77289`, `R-MMU-77289`, `R-RNO-77289` | Mitochondrial Fatty Acid Beta-Oxidation, a 4-level chain under Metabolism in three species. Mouse and rat share the human number and name, so it is the human-equivalent case |
| `R-HSA-5602358` | a chain under Disease, which Reactome does not infer to other species |
| `R-HSA-354192` | Integrin signaling, which has two parents under two roots (Signal Transduction; Hemostasis via Platelet Aggregation) |
| `R-GGA-351470` | a curated chicken pathway with no `R-HSA` twin of the same number, so it has no human equivalent |

**Also kept:** every relation line whose parent and child are both kept.

Because the slice is closed under ancestors, every kept pathway's parents are present. So `ReactomeRelease.Ancestors`
returns the same set here as against the full release.

**Result:** 25 pathways (1,365 B, sha256 `7273f1ca3324ece67f7c6eb97742a2fa6a2ef17af484ed7876dd2d3e9e8756df`), 19 pairs
(501 B, sha256 `bb938866a859265afcd5fa7ae5e80c29d0d2e5fff9bc325ffa0c5331e0d2489c`), and 7 roots.

**Refresh** by re-fetching both URLs for the pinned release and re-applying the rule above. A refresh is a reviewable
change to this repository.
