# PRIDE zip tail fixtures

Each `<zip>.tail` is the last 65,577 bytes of a zip in PRIDE project PXD009815, exactly as PRIDE's FTP host served them
over HTTPS on 2026-11-18 to `GET https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD009815/<zip>` with
`Range: bytes=-65577`. Nothing is changed. They hold only the end of the last entry's compressed data, the central
directory (entry names, sizes, CRC-32s, offsets) and the end records.

| file | the zip | answer |
|---|---|---|
| `mzdb_files.zip.tail` | 22,821,363,261 bytes; ZIP64; 41 entries, 32 of them with the ZIP64 extra field | `206`, `Content-Range: bytes 22821297684-22821363260/22821363261`, `ETag: "550425e3d-5d2ba86de0300"` |
| `MaxQuant_Results_50f-25f.zip.tail` | 25,399,698 bytes; 6 entries, all deflated | `206`, `Content-Range: bytes 25334121-25399697/25399698`, `ETag: "1839192-5d2ba967bd080"` |

The tests serve each tail at its real offset in a file of its real length, so the client reads them as it reads PRIDE.
