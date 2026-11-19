# MassIVE fixtures

## PROXI records (M1a)

Each `<MSV>.proxi.json` is MassIVE's PROXI answer to
`GET https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/<MSV>`, fetched 2026-10-18. Every byte is as served,
with one exception: the `contacts` array is emptied (`[]`), because it carries names and email addresses. Nothing
else is changed. The non-ASCII bytes the charset tests need are in the summaries and references, outside `contacts`.

| file | why it is here |
|---|---|
| `MSV000096226.proxi.json` | an ordinary record: a PXD, keywords, one species group, "no associated manuscript" |
| `MSV000078494.proxi.json` | the string `"null"` for species and instrument; no PXD, no keywords, no modifications |
| `MSV000088323.proxi.json` | Latin-1 bytes (`Rapp\xE9`) under the server's `charset=ISO-8859-1` label |
| `MSV000079927.proxi.json` | UTF-8 bytes (`CAUTI\xC2\xB4s`) under the same label; four species groups |
| `MSV000097834.proxi.json` | a title that begins `"Title: "` |

## Listings (M1b)

`listings/<MSV>.listing.tsv` is the listing `MassiveClient.GetDatasetFilesAsync` (M1a, `88bc5ca9f`) returned from
`massive-ftp.ucsd.edu` on 2026-11-10: one row per entry, `path`, `size` (empty for a link), `is_link` and `modified`
(as LIST shows it). Paths only, nothing personal. They are the run rule's, the SDRF signal's and the design-table
rule's cases:

| deposit | FTP root | case |
|---|---|---|
| MSV000093798 | `/v06` | 192 `.raw` runs; the `ccms_peak` link |
| MSV000080763 | `/x01` | no runs: its spectra are only in `ccms_peak` (unlistable), its originals in a `.rar` |
| MSV000098930 | `/v10` | 508 `raw/RAW/*.zip` runs, and a `search/` zip that is not one; an annotation `.csv` |
| MSV000090777 | `/v05` | SCIEX: 50 runs of `.wiff` + `.wiff.scan` + `.wiff2` + `.timeseries.data`; 4 `.DS_Store`; 27 `ConditionSetup.tsv`; an update |
| MSV000092529 | `/v06` | Bruker: 6 `.d` folders listed file by file |
| MSV000096226 | `/v07` | two `ConditionSetup.tsv` |
| MSV000096205, MSV000100606, MSV000102340, MSV000103504, MSV000102494, MSV000094031, MSV000097296 | `/v07`, `/v12`, `/v13`, `/v15`, `/v13`, `/v07`, `/v09` | the seven deposited SDRFs sdrf knows (sdrf 026); MSV000103504 also has `experiment_annotation.tsv` |

## params.xml (M1b)

`<MSV>.params.xml` is the dataset's `ccms_parameters/params.xml`, fetched over FTPS on 2026-11-10;
`MSV000090777_update.params.xml` is its update's (`updates/2023-01-30_jburton1_42cc6275/ccms_parameters/`). Every
byte is as served except the personal values, which are replaced with fakes that tests check never leave the parser:
every email address becomes `submitter@example.org`, the `dataset.pi` value becomes
`Redacted Name|pi@example.org|Redacted Institution|Redacted Country`, and the `dataset.password` value becomes
`redacted-password`.

## Search answers (M1b)

- `massive_search_zebrafish_p<N>.<status>.json`: MassIVE PROXI
  `datasets?resultType=compact&pageSize=30&pageNumber=<N>&search=zebrafish`, 2026-11-10. Pages 1-3 hold 30, 30 and
  7 rows; page 4 is MassIVE's 404 "No data found". Each row's `contacts` is emptied (emails); the JSON is otherwise
  re-written compactly in the encoding it came in. Page 1 arrived with one Latin-1 byte, but it was inside a contact,
  so emptying `contacts` removed it: all three pages are now UTF-8. (The charset rule is pinned by the PROXI record
  fixtures above, which go through the same decoder.) 67 datasets, 19 with a PXD, first MSV000078304, last MSV000101937.
- `proteomecentral_search_zebrafish_massive_p<N>.200.json`: ProteomeCentral PROXI
  `datasets?resultType=compact&pageSize=5&pageNumber=<N>&search=zebrafish&repository=MassIVE`, 2026-11-10. Pages 1-4
  hold 5, 5, 1 and 0 rows of 11. The "lab head" column (a person's name) is replaced with `Redacted Lab Head`; the
  JSON is otherwise re-written compactly.
