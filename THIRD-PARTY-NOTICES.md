# Third-party notices

mzLib is distributed under LGPL-3.0-or-later. It additionally **embeds and redistributes** the
third-party data files below. Their licences are separate from mzLib's and are reproduced here
because NuGet consumers receive the compiled assemblies rather than this repository.

---

## PSI-MS controlled vocabulary — `psi-ms.obo`

- **Publisher:** HUPO Proteomics Standards Initiative (PSI), Mass Spectrometry Working Group
- **Source:** https://github.com/HUPO-PSI/psi-ms-CV — https://purl.obolibrary.org/obo/ms.obo
- **Pinned version:** `data-version: 4.1.258`
- **Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0),
  https://creativecommons.org/licenses/by/4.0/
- **Embedded in:** `UsefulProteomicsDatabases.dll`
  (`UsefulProteomicsDatabases/Resources/psi-ms.obo`)

## PRIDE controlled vocabulary — `pride_cv.obo`

- **Publisher:** EMBL-EBI PRIDE Team / PRIDE-Utilities
- **Source:** https://github.com/PRIDE-Utilities/pride-ontology
- **Pinned version:** `data-version: releases/2026-06-19`
- **Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0),
  https://creativecommons.org/licenses/by/4.0/
- **Embedded in:** `UsefulProteomicsDatabases.dll`
  (`UsefulProteomicsDatabases/Resources/pride_cv.obo`)

---

## Pre-existing embedded data

Recorded here for completeness; these predate the notices file.

| File | Source | Embedded in |
|---|---|---|
| `unimod.xml` | Unimod, http://www.unimod.org/ | `Omics.dll` (`Omics/Resources/`) |
| `PSI-MOD.obo.xml` | HUPO-PSI PSI-MOD, https://github.com/HUPO-PSI/psi-mod-CV | `Omics.dll` |
| `ptmlist.txt` | UniProt, https://www.uniprot.org/ | `Omics.dll` |

Fetch procedures for the two vocabularies added alongside this file are in
`mzLib/UsefulProteomicsDatabases/Resources/PROVENANCE.md`.
