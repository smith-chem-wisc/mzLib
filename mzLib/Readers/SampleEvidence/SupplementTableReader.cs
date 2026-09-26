using System.Globalization;
using System.IO.Compression;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using CsvHelper;
using CsvHelper.Configuration;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// One table from a paper's supplementary file: its header (found under any title lines), its data rows padded to
    /// one width, and where each row sits in the source, so a value can be cited back to its cell.
    /// </summary>
    /// <param name="File">The file name, without its folder.</param>
    /// <param name="Sheet">The sheet name for a workbook; <c>table N</c> for a document; empty for delimited text.</param>
    /// <param name="Title">The text of any rows above the header (a caption such as "Table S1. Clinical data"), else empty.</param>
    /// <param name="Header">The header cells.</param>
    /// <param name="Rows">The data rows, none empty, each as wide as the widest row.</param>
    /// <param name="RowNumbers">For each data row, its 1-based row number in the source.</param>
    /// <param name="Truncated">True when the source had more rows than the reader's cap.</param>
    internal sealed record SupplementTable(
        string File,
        string Sheet,
        string Title,
        IReadOnlyList<string> Header,
        IReadOnlyList<IReadOnlyList<string>> Rows,
        IReadOnlyList<int> RowNumbers,
        bool Truncated)
    {
        /// <summary>Where a cell sits: <c>file!sheet!R&lt;row&gt;C&lt;column&gt;</c>, 1-based, for a source reference.</summary>
        public string Locator(int row, int column) =>
            (Sheet.Length > 0 ? $"{File}!{Sheet}" : File) + $"!R{RowNumbers[row]}C{column + 1}";
    }

    /// <summary>
    /// Reads the tables in a supplementary file as rows of strings (sdrf design SAMPLE-EVIDENCE.md, E2; D41: in mzLib).
    ///
    /// <para><c>.xlsx</c>/<c>.xlsm</c> and <c>.docx</c> are zipped XML and are read as such, with no spreadsheet
    /// package; <c>.csv</c>, <c>.tsv</c> and <c>.txt</c> (ISA-Tab included) are read quote-aware with the separator
    /// sniffed. Cells are text as written: a whole number stored as <c>3.0</c> becomes <c>3</c>, and nothing else is
    /// reformatted (a date stays the serial number the file stores). A PDF's tables are found from word positions
    /// (<c>SupplementTableReader.Pdf.cs</c>). A format with no reader here (<c>.xls</c>)
    /// gives no tables; a file that claims a format and is not in it throws <see cref="MzLibException"/>.</para>
    /// </summary>
    internal static partial class SupplementTableReader
    {
        /// <summary>Rows read per table before stopping. Sample tables are short; result tables can run to 100,000.</summary>
        internal const int DefaultMaxRows = 5000;

        private const string SheetNs = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
        private const string RelNs = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";
        private const string PackageRelNs = "http://schemas.openxmlformats.org/package/2006/relationships";
        private const string WordNs = "http://schemas.openxmlformats.org/wordprocessingml/2006/main";

        /// <summary>Every table in the file, in document order. Never null.</summary>
        internal static IReadOnlyList<SupplementTable> Read(string path, int maxRows = DefaultMaxRows)
        {
            if (path == null) throw new ArgumentNullException(nameof(path));
            if (maxRows < 1) throw new ArgumentOutOfRangeException(nameof(maxRows), maxRows, "At least one row.");
            string name = Path.GetFileName(path);
            try
            {
                return Path.GetExtension(path).ToLowerInvariant() switch
                {
                    ".xlsx" or ".xlsm" => ReadWorkbook(path, name, maxRows),
                    ".docx" => ReadDocument(path, name, maxRows),
                    ".csv" or ".tsv" or ".txt" => ReadDelimited(path, name, maxRows),
                    ".pdf" => ReadPdf(path, name, maxRows),
                    _ => Array.Empty<SupplementTable>()
                };
            }
            catch (Exception e) when (e is InvalidDataException or XmlException or CsvHelperException)
            {
                throw new MzLibException($"Supplementary file '{name}' could not be read as its extension says: {e.Message}", e);
            }
        }

        // ---------------- tables from raw rows ----------------

        /// <summary>Drops empty rows, pads to one width, finds the header under any title lines.</summary>
        private static SupplementTable? Table(string file, string sheet, List<(int Number, List<string> Cells)> raw, bool truncated)
        {
            var rows = raw.Where(r => r.Cells.Any(c => c.Length > 0)).ToList();
            if (rows.Count < 2) return null;
            int width = rows.Max(r => r.Cells.Count);
            foreach (var r in rows) while (r.Cells.Count < width) r.Cells.Add("");
            int h = HeaderRow(rows.Select(r => r.Cells).ToList());
            var body = rows.Skip(h + 1).ToList();
            if (body.Count == 0) return null;
            string title = string.Join(" ", rows.Take(h).SelectMany(r => r.Cells).Where(c => c.Length > 0));
            return new SupplementTable(file, sheet, title, rows[h].Cells, body.Select(r => (IReadOnlyList<string>)r.Cells).ToList(),
                body.Select(r => r.Number).ToList(), truncated);
        }

        private static readonly Regex Numeric = new(@"^[-+.\deE%]+$", RegexOptions.Compiled);

        /// <summary>
        /// The first of the top 15 rows with at least two filled cells, mostly text, and at least half as full as the
        /// rows under it -- so a one-cell caption above the header is a title, not the header.
        /// </summary>
        private static int HeaderRow(IReadOnlyList<List<string>> rows)
        {
            for (int i = 0; i < Math.Min(15, rows.Count); i++)
            {
                var filled = rows[i].Where(c => c.Length > 0).ToList();
                if (filled.Count < 2 || filled.Count(c => !Numeric.IsMatch(c)) < 0.6 * filled.Count) continue;
                int below = rows.Skip(i + 1).Take(3).Select(r => r.Count(c => c.Length > 0)).DefaultIfEmpty(0).Max();
                if (filled.Count >= 0.5 * below) return i;
            }
            return 0;
        }

        private static string Clean(string s) => Regex.Replace(s, @"\s+", " ").Trim();

        // ---------------- .xlsx ----------------

        private static IReadOnlyList<SupplementTable> ReadWorkbook(string path, string file, int maxRows)
        {
            using var zip = ZipFile.OpenRead(path);
            var shared = SharedStrings(zip);
            var targets = new Dictionary<string, string>(StringComparer.Ordinal);
            var rels = zip.GetEntry("xl/_rels/workbook.xml.rels")
                ?? throw new InvalidDataException("no xl/_rels/workbook.xml.rels");
            using (var r = XmlReader.Create(rels.Open()))
                while (r.Read())
                    if (r.NodeType == XmlNodeType.Element && r.LocalName == "Relationship" && r.GetAttribute("Id") is { } id && r.GetAttribute("Target") is { } t)
                        targets[id] = t.StartsWith('/') ? t.TrimStart('/') : "xl/" + t;

            var workbook = zip.GetEntry("xl/workbook.xml") ?? throw new InvalidDataException("no xl/workbook.xml");
            var sheets = new List<(string Name, string Entry)>();
            using (var r = XmlReader.Create(workbook.Open()))
                while (r.Read())
                    if (r.NodeType == XmlNodeType.Element && r.LocalName == "sheet" && r.NamespaceURI == SheetNs
                        && r.GetAttribute("id", RelNs) is { } rid && targets.TryGetValue(rid, out var entry))
                        sheets.Add((r.GetAttribute("name") ?? "", entry));

            var tables = new List<SupplementTable>();
            foreach (var (name, entryName) in sheets)
            {
                var entry = zip.GetEntry(entryName);
                if (entry == null) continue;
                var (raw, truncated) = SheetRows(entry, shared, maxRows);
                if (Table(file, name, raw, truncated) is { } t) tables.Add(t);
            }
            return tables;
        }

        private static List<string> SharedStrings(ZipArchive zip)
        {
            var strings = new List<string>();
            var entry = zip.GetEntry("xl/sharedStrings.xml");
            if (entry == null) return strings;
            using var r = XmlReader.Create(entry.Open());
            StringBuilder? current = null;
            r.Read();
            while (!r.EOF)
            {
                // ReadElementContentAsString leaves the reader on the NEXT node, so only Read() when nothing was consumed.
                if (r.NodeType == XmlNodeType.Element && r.LocalName == "rPh") { r.Skip(); continue; }   // a reading aid, not the text
                if (r.NodeType == XmlNodeType.Element && r.LocalName == "t" && current != null) { current.Append(r.ReadElementContentAsString()); continue; }
                if (r.NodeType == XmlNodeType.Element && r.LocalName == "si") current = new StringBuilder();
                else if (r.NodeType == XmlNodeType.EndElement && r.LocalName == "si" && current != null) { strings.Add(Clean(current.ToString())); current = null; }
                r.Read();
            }
            return strings;
        }

        private static (List<(int, List<string>)> Rows, bool Truncated) SheetRows(ZipArchiveEntry entry, List<string> shared, int maxRows)
        {
            var rows = new List<(int, List<string>)>();
            using var r = XmlReader.Create(entry.Open());
            List<string>? cells = null;
            int number = 0;
            while (r.Read())
            {
                if (r.NodeType == XmlNodeType.Element && r.LocalName == "row")
                {
                    if (rows.Count >= maxRows) return (rows, true);
                    number = int.TryParse(r.GetAttribute("r"), NumberStyles.None, CultureInfo.InvariantCulture, out int n) ? n : number + 1;
                    cells = new List<string>();
                    rows.Add((number, cells));
                    if (r.IsEmptyElement) cells = null;
                }
                else if (r.NodeType == XmlNodeType.Element && r.LocalName == "c" && cells != null)
                {
                    string reference = r.GetAttribute("r") ?? "";
                    string type = r.GetAttribute("t") ?? "n";
                    int column = ColumnIndex(reference, cells.Count);
                    string value = r.IsEmptyElement ? "" : CellValue(r, type, shared);
                    while (cells.Count < column) cells.Add("");
                    if (cells.Count == column) cells.Add(value); else cells[column] = value;
                }
            }
            return (rows, false);
        }

        /// <summary>The cell's text: shared string, inline string, boolean, or the stored number/text.</summary>
        private static string CellValue(XmlReader r, string type, List<string> shared)
        {
            string v = "", inline = "";
            using (var sub = r.ReadSubtree())
            {
                sub.Read();
                sub.Read();
                while (!sub.EOF)
                {
                    if (sub.NodeType == XmlNodeType.Element && sub.LocalName == "v") { v = sub.ReadElementContentAsString(); continue; }
                    if (sub.NodeType == XmlNodeType.Element && sub.LocalName == "t") { inline += sub.ReadElementContentAsString(); continue; }
                    if (sub.NodeType == XmlNodeType.Element && sub.LocalName == "rPh") { sub.Skip(); continue; }
                    sub.Read();
                }
            }
            return type switch
            {
                "s" => int.TryParse(v, NumberStyles.None, CultureInfo.InvariantCulture, out int i) && i < shared.Count ? shared[i] : "",
                "inlineStr" => Clean(inline),
                "b" => v == "1" ? "TRUE" : "FALSE",
                "n" => Number(v),
                _ => Clean(v)
            };
        }

        /// <summary>A whole number stored as <c>3.0</c> is <c>3</c>; anything else is kept as stored.</summary>
        private static string Number(string v)
        {
            if (double.TryParse(v, NumberStyles.Float, CultureInfo.InvariantCulture, out double d) && d == Math.Floor(d) && Math.Abs(d) < 1e15)
                return ((long)d).ToString(CultureInfo.InvariantCulture);
            return v.Trim();
        }

        /// <summary>Column letters of a reference (<c>C5</c> -> 2); the next column when the cell has no reference.</summary>
        private static int ColumnIndex(string reference, int next)
        {
            int col = 0, k = 0;
            while (k < reference.Length && char.IsAsciiLetter(reference[k])) col = col * 26 + (char.ToUpperInvariant(reference[k++]) - 'A' + 1);
            return k == 0 ? next : col - 1;
        }

        // ---------------- .docx ----------------

        private static IReadOnlyList<SupplementTable> ReadDocument(string path, string file, int maxRows)
        {
            using var zip = ZipFile.OpenRead(path);
            var entry = zip.GetEntry("word/document.xml") ?? throw new InvalidDataException("no word/document.xml");
            var tables = new List<SupplementTable>();
            using var r = XmlReader.Create(entry.Open());
            int k = 0;
            while (r.Read())
            {
                if (r.NodeType != XmlNodeType.Element || r.LocalName != "tbl" || r.NamespaceURI != WordNs) continue;
                k++;
                var raw = new List<(int, List<string>)>();
                bool truncated = false;
                using (var tbl = r.ReadSubtree())
                {
                    tbl.Read();
                    int depth = tbl.Depth;
                    while (tbl.Read())
                    {
                        // Only the table's own rows: a nested table's rows belong to its cell's text.
                        if (tbl.NodeType == XmlNodeType.Element && tbl.LocalName == "tr" && tbl.Depth == depth + 1)
                        {
                            if (raw.Count >= maxRows) { truncated = true; break; }
                            raw.Add((raw.Count + 1, RowCells(tbl)));
                        }
                    }
                }
                if (Table(file, $"table {k}", raw, truncated) is { } t) tables.Add(t);
            }
            return tables;
        }

        private static List<string> RowCells(XmlReader tr)
        {
            var cells = new List<string>();
            using var row = tr.ReadSubtree();
            row.Read();
            int depth = row.Depth;
            while (row.Read())
            {
                if (row.NodeType != XmlNodeType.Element || row.LocalName != "tc" || row.Depth != depth + 1) continue;
                var text = new StringBuilder();
                using var tc = row.ReadSubtree();
                tc.Read();
                while (!tc.EOF)
                {
                    if (tc.NodeType == XmlNodeType.Element && tc.LocalName == "t") { text.Append(tc.ReadElementContentAsString()); continue; }
                    if (tc.NodeType == XmlNodeType.EndElement && tc.LocalName == "p") text.Append(' ');
                    tc.Read();
                }
                cells.Add(Clean(text.ToString()));
            }
            return cells;
        }

        // ---------------- delimited text ----------------

        private static IReadOnlyList<SupplementTable> ReadDelimited(string path, string file, int maxRows)
        {
            var head = File.ReadLines(path).Take(20).ToList();
            string delimiter = new[] { "\t", ",", ";" }.OrderByDescending(d => head.Sum(l => l.Split(d).Length - 1)).First();
            var config = new CsvConfiguration(CultureInfo.InvariantCulture)
            {
                Delimiter = delimiter,
                HasHeaderRecord = false,
                BadDataFound = null,
                MissingFieldFound = null,
                DetectColumnCountChanges = false,
                // A tab-separated file is not quoted CSV: a stray quote there is text, not a field boundary.
                Mode = delimiter == "\t" && !head.Any(l => l.StartsWith('"')) ? CsvMode.NoEscape : CsvMode.RFC4180
            };
            var raw = new List<(int, List<string>)>();
            bool truncated = false;
            using var reader = new StreamReader(path, Encoding.UTF8, detectEncodingFromByteOrderMarks: true);
            using var csv = new CsvParser(reader, config);
            while (csv.Read())
            {
                if (raw.Count >= maxRows) { truncated = true; break; }
                raw.Add((csv.Row, csv.Record?.Select(Clean).ToList() ?? new List<string>()));
            }
            var table = Table(file, "", raw, truncated);
            return table == null ? Array.Empty<SupplementTable>() : new[] { table };
        }
    }
}
