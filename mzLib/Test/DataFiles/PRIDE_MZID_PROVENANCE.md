# PRIDE mzIdentML fixtures

The `PXD*.mzid` files next to this note are real mzIdentML output from five different writers,
downloaded from PRIDE Archive and trimmed. They exist because the published HUPO-PSI examples do not
reach the CV spellings and document shapes real writers produce. See `MzidIdentificationsTests` and
`TestMzIdentMLResultFile` for what each one pins. The first five are listed here; the two modification
fixtures and one derived `.gz` are described further down.

| Fixture | PRIDE project | Source file | Writer | Source line endings | SHA-256 of the file as downloaded |
|---|---|---|---|---|---|
| `PXD078927_msgf_1_1_0.mzid` | PXD078927 | `2026/05/PXD078927/MHCI_pep.mzid` | MS-GF+ | LF | `7fcf43b7456ffaa5e105b30332e86ee06f32ba94505666dd55e554682d676477` |
| `PXD000783_scaffold_1_1_0.mzid` | PXD000783 | `2014/06/PXD000783/mascot_daemon_merge_F008897.mzid.gz` | Scaffold (Mascot) | CRLF | `0d96698fefa5288744afccf8778e0e0373830ce2f619a4eb199c140fa985bdd9` |
| `PXD019591_mascotparser_1_1_0.mzid` | PXD019591 | `2021/02/PXD019591/F002080.mzid.gz` | Mascot Parser | CRLF | `17fdb77b853d0a6b1232ffbad5d14acdc0ff3a07ec541c373caf884b8068a8dc` |
| `PXD019733_proteomediscoverer_1_1_0.mzid` | PXD019733 | `2021/06/PXD019733/500to10K_Plate1col3to8.mzid.gz` | Proteome Discoverer | CRLF | `9a5589c6a23ad7fec32d41a21bb834165a633c2e02e69486f617299319b7be53` |
| `PXD070193_xifdr_1_2_0.mzid` | PXD070193 | `2026/06/PXD070193/17804_MCM2-7_sld3-7_Tirs_filter_1linkFDRBB.mzid` | xiFDR | LF | `fe3b65fc290ad53efb8805c34e174db9fb0b222c30d53027c751a3fe5053174a` |

Source paths are relative to `https://ftp.pride.ebi.ac.uk/pride/data/archive/`. The checksums are of
the bytes as served (the `.gz` itself for the three compressed files), taken 2026-09-17.

## The trim rule

- Keep everything up to and including the first `SpectrumIdentificationList`'s opening tag, then that
  list's first two or three `SpectrumIdentificationResult`s.
- From `SequenceCollection`, keep only the `DBSequence`, `Peptide` and `PeptideEvidence` entries those
  results reference.
- Drop the rest of `AnalysisData`, including any `ProteinDetectionList`. Drop the `AnalysisCollection`
  elements that pointed into what was dropped: `ProteinDetection`, and in xiFDR the eight
  `SpectrumIdentification`s whose lists were not kept. Every `*_ref` left in a fixture resolves.
- Close the open elements. Those closing tags are the only text not in the source file.

Every other line is the source line **byte for byte, line ending included**. `.gitattributes` marks
these files `-text` so that git stores and checks out those bytes unchanged. A check against the
downloaded file matches every non-blank fixture line except those closing tags.

The unreferenced `SpectraData` entries are kept on purpose. xiFDR's results reference `SpectraData`
index 5 of 9, which is the shape `Ms2SpectrumID`'s `SpectraData[0]` shortcut gets wrong.

## Modification fixtures

Two more captures come from the same downloads as the Scaffold and Mascot Parser fixtures above. They keep
different results, chosen for the modifications their peptides carry. `TestMzIdentMLResultFile` reads them.

| Fixture | Source file (as above) | Results kept | Modifications reached |
|---|---|---|---|
| `PXD000783_scaffold_mods_1_1_0.mzid` | `mascot_daemon_merge_F008897.mzid.gz` | `Spec_524963`, `Spec_525064`, `Spec_525101`, `Spec_525280`, `Spec_525765` | UNIMOD:21 Phospho, UNIMOD:4 Carbamidomethyl, UNIMOD:35 Oxidation, UNIMOD:7 Deamidated, UNIMOD:1 Acetyl at the N-terminus (location 0, no `residues`) |
| `PXD019591_mascotparser_mods_1_1_0.mzid` | `F002080.mzid.gz` | `SIR_47`, `SIR_61`, `SIR_157`, `SIR_200` | UNIMOD:21 Phospho, UNIMOD:35 Oxidation, UNIMOD:4 Carbamidomethyl, UNIMOD:28 Gln->pyro-Glu at the N-terminus |

Their trim rule differs from the one above in one respect. The kept results are not the first ones:
walking the first `SpectrumIdentificationList` in order, a result is kept when one of its items references
a peptide carrying a (UNIMOD accession, N-terminal or not) pair that no earlier kept result carried. Every
other rule is the same. The referenced sequence entries are kept, `ProteinDetection` is dropped, the closing
tags are the only added text, and the source line endings are kept. The same line-by-line check against the
download passes.

## Derived fixture

`PXD078927_msgf_1_1_0.mzid.gz` is **not** a PRIDE download. It is `PXD078927_msgf_1_1_0.mzid` compressed
with gzip (no file name, modification time 0), so the compressed read path has a fixture whose
uncompressed bytes are pinned. Reading it must give exactly what reading the `.mzid` gives.
`.gitattributes` marks it `binary`.

The live canary, `PrideMzIdentMLLiveTests`, downloads a real compressed file instead:
`2014/07/PXD000710/ma190_19_tandem_pproph.pep.mzid.gz` (X!Tandem + PeptideProphet, 8,563 bytes, SHA-256
`80bee1aa7c6637a8e780d6f7528cae0ac82e031b06f8791879d936a905cae99f` on 2026-09-17).
