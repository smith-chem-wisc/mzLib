# MassIVE PROXI fixtures

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
