using UglyToad.PdfPig;
using UglyToad.PdfPig.Content;
using UglyToad.PdfPig.Core;
using UglyToad.PdfPig.Exceptions;

namespace Readers
{
    internal static partial class SupplementTableReader
    {
        /// <summary>Pages read per PDF. Supplementary PDFs with sample tables are short; a 300-page atlas is not one.</summary>
        internal const int MaxPdfPages = 200;

        /// <summary>Lines further apart than this many line heights end a table (padded rows run to about 3).</summary>
        private const double RowGap = 4.0;

        /// <summary>
        /// A gap wider than this many line heights starts a new cell. A word space is about 0.3; PXD007160's channel map
        /// puts its first sample ID 0.97 of a line height after the machine number, so 1.0 fused two columns.
        /// </summary>
        private const double CellGap = 0.8;

        private sealed record PdfWord(double Left, double Right, string Text);

        private sealed record PdfCell(double Left, double Right, List<PdfWord> Words)
        {
            public string Text => string.Join(" ", Words.Select(w => w.Text));
        }

        private sealed record PdfLine(double Y, double Height, List<PdfCell> Cells);

        /// <summary>
        /// Tables laid out in a PDF's text, found from word positions (PdfPig): words become lines by their height,
        /// a line splits into cells at a gap wider than a character's height, and a table is a run of at least three
        /// lines of two or more cells. Its columns are where those cells overlap across the run. A caption line
        /// just above ("Table S1 ...") is the table's title. Scanned pages (no text) give nothing.
        /// Measured need (sdrf G36): 28 of 670 deposits with supplement tables have a sample table in a PDF, 21 only
        /// there -- PXD007160's TMT channel map among them.
        /// </summary>
        private static IReadOnlyList<SupplementTable> ReadPdf(string path, string file, int maxRows)
        {
            var tables = new List<SupplementTable>();
            try
            {
                using var document = PdfDocument.Open(path);
                int pageNumber = 0;
                foreach (var page in document.GetPages())
                {
                    if (++pageNumber > MaxPdfPages) break;
                    var lines = Lines(page);
                    int k = 0;
                    foreach (var (start, found, caption) in Blocks(lines))
                    {
                        var columns = Columns(found);
                        if (columns.Count < 2) continue;
                        var block = WithHeaderAbove(lines, start, found, columns, caption);
                        var raw = new List<(int Number, List<string> Cells)>();
                        if (caption != null) raw.Add((0, new List<string> { caption }));
                        bool truncated = block.Count > maxRows;
                        var body = new List<(int Number, List<string> Cells)>();
                        for (int i = 0; i < Math.Min(block.Count, maxRows); i++)
                            body.Add((i + 1, Assign(block[i], columns)));
                        raw.AddRange(MergeWrappedHeader(body, columns.Count));
                        var table = Table(file, $"page {pageNumber} table {k + 1}", raw, truncated);
                        if (table == null) continue;
                        tables.Add(table);
                        k++;
                    }
                }
            }
            catch (Exception e) when (e is PdfDocumentFormatException or PdfDocumentEncryptedException)
            {
                throw new InvalidDataException(e.Message, e);
            }
            return tables;
        }

        /// <summary>Words grouped into lines, top of the page first, each line split into cells at wide gaps.</summary>
        private static List<PdfLine> Lines(Page page)
        {
            var words = page.GetWords().Where(w => !string.IsNullOrWhiteSpace(w.Text) && w.BoundingBox.Height > 0)
                .OrderByDescending(w => w.BoundingBox.Bottom).ToList();
            var groups = new List<List<Word>>();
            double bottom = 0, top = 0;
            foreach (var w in words)
            {
                var box = w.BoundingBox;
                var last = groups.LastOrDefault();
                // Same line when the word overlaps the line's vertical span by half its height or more, so a superscript
                // (the 2 of m2) joins its line instead of making one of its own.
                if (last != null && Math.Min(top, box.Top) - Math.Max(bottom, box.Bottom) >= 0.5 * Math.Min(box.Height, top - bottom))
                {
                    last.Add(w);
                    bottom = Math.Min(bottom, box.Bottom);
                    top = Math.Max(top, box.Top);
                }
                else
                {
                    groups.Add(new List<Word> { w });
                    (bottom, top) = (box.Bottom, box.Top);
                }
            }
            var lines = new List<PdfLine>();
            foreach (var g in groups)
            {
                var sorted = g.OrderBy(w => w.BoundingBox.Left).ToList();
                double height = sorted.Select(w => w.BoundingBox.Height).OrderBy(h => h).ElementAt(sorted.Count / 2);
                var cells = new List<PdfCell>();
                var cellWords = new List<PdfWord> { new(sorted[0].BoundingBox.Left, sorted[0].BoundingBox.Right, sorted[0].Text) };
                double left = sorted[0].BoundingBox.Left, right = sorted[0].BoundingBox.Right;
                for (int i = 1; i < sorted.Count; i++)
                {
                    var b = sorted[i].BoundingBox;
                    if (b.Left - right > CellGap * height)
                    {
                        cells.Add(new PdfCell(left, right, cellWords));
                        cellWords = new List<PdfWord>();
                        left = b.Left;
                    }
                    cellWords.Add(new PdfWord(b.Left, b.Right, sorted[i].Text));
                    right = Math.Max(right, b.Right);
                }
                cells.Add(new PdfCell(left, right, cellWords));
                lines.Add(new PdfLine(g.OrderByDescending(w => w.BoundingBox.Height).First().BoundingBox.Bottom, height, cells));
            }
            return lines;
        }

        private static readonly System.Text.RegularExpressions.Regex Caption =
            new(@"^(supplementa(l|ry) )?table\b", System.Text.RegularExpressions.RegexOptions.IgnoreCase | System.Text.RegularExpressions.RegexOptions.Compiled);

        /// <summary>Runs of three or more close lines with two or more cells, each with the caption line above it, if any.</summary>
        private static IEnumerable<(int Start, List<PdfLine> Block, string? Caption)> Blocks(List<PdfLine> lines)
        {
            int i = 0;
            while (i < lines.Count)
            {
                if (lines[i].Cells.Count < 2) { i++; continue; }
                int j = i + 1;
                // One-cell lines inside a run are a wrapped cell (up to three), kept while the run continues after them.
                while (j < lines.Count && Close(lines[j - 1], lines[j]) && (lines[j].Cells.Count >= 2 || ResumesAfter(lines, j)))
                    j++;
                var block = lines.GetRange(i, j - i);
                if (block.Count(l => l.Cells.Count >= 2) >= 3)
                {
                    string? caption = null;
                    // The caption is usually above; some journals set it below (PXD007160's "Supplementary table 1").
                    foreach (int c in Enumerable.Range(1, 3).Select(d => i - d).Concat(Enumerable.Range(0, 3).Select(d => j + d)))
                    {
                        if (c < 0 || c >= lines.Count) continue;
                        string text = string.Join(" ", lines[c].Cells.Select(x => x.Text));
                        if (Caption.IsMatch(text)) { caption = text; break; }
                    }
                    yield return (i, block, caption);
                }
                i = j;
            }
        }

        /// <summary>
        /// Up to two one-cell lines just above a table, each inside a single column, are the upper half of wrapped
        /// header cells ("Pediatric" over "Age Group") and join the table. Prose and captions run across columns.
        /// </summary>
        private static List<PdfLine> WithHeaderAbove(List<PdfLine> lines, int start, List<PdfLine> block, List<(double Left, double Right)> columns, string? caption)
        {
            var result = new List<PdfLine>(block);
            for (int i = start - 1; i >= Math.Max(0, start - 2); i--)
            {
                var line = lines[i];
                if (line.Cells.Count != 1 || !Close(line, result[0])) break;
                var c = line.Cells[0];
                if (Caption.IsMatch(c.Text)) break;
                // The second line of a caption ("... used for" over "analysis", PXD037923) is not a header.
                if (i > 0 && Close(lines[i - 1], line) && Caption.IsMatch(string.Join(" ", lines[i - 1].Cells.Select(x => x.Text)))) break;
                if (columns.Count(k => Math.Min(c.Right, k.Right) - Math.Max(c.Left, k.Left) > 0) != 1) break;
                result.Insert(0, line);
            }
            return result;
        }

        private static bool Close(PdfLine above, PdfLine below) => above.Y - below.Y < RowGap * Math.Max(above.Height, below.Height);

        /// <summary>Whether a run of one-cell lines starting at <paramref name="j"/> ends, within three lines, in a table line.</summary>
        private static bool ResumesAfter(List<PdfLine> lines, int j)
        {
            for (int k = j + 1; k < Math.Min(lines.Count, j + 4); k++)
            {
                if (!Close(lines[k - 1], lines[k])) return false;
                if (lines[k].Cells.Count >= 2) return true;
            }
            return false;
        }

        /// <summary>
        /// Column spans: the union of overlapping cell spans, left to right, taken from the table's fullest lines (at
        /// the most cells on a line), so a header cell spanning a group of columns, or one fused cell,
        /// cannot join columns. Every other line's cells are then placed by overlap.
        /// </summary>
        private static List<(double Left, double Right)> Columns(List<PdfLine> block)
        {
            int most = block.Max(l => l.Cells.Count);
            var base_ = Merge(block.Where(l => l.Cells.Count == most).SelectMany(l => l.Cells).Select(c => (c.Left, c.Right)));
            // A column only the header or wrapped lines fill (PXD037923's "Colonoscopy Location") is added where it
            // overlaps none of the fullest lines' columns.
            var extra = block.Where(l => l.Cells.Count >= 2 && l.Cells.Count < most).SelectMany(l => l.Cells)
                .Select(c => (c.Left, c.Right))
                .Where(c => !base_.Any(k => Math.Min(c.Right, k.Right) - Math.Max(c.Left, k.Left) > 0));
            return Merge(base_.Concat(extra));
        }

        private static List<(double Left, double Right)> Merge(IEnumerable<(double Left, double Right)> spans)
        {
            var merged = new List<(double Left, double Right)>();
            foreach (var s in spans.OrderBy(s => s.Left))
            {
                if (merged.Count > 0 && s.Left <= merged[^1].Right)
                    merged[^1] = (merged[^1].Left, Math.Max(merged[^1].Right, s.Right));
                else
                    merged.Add(s);
            }
            return merged;
        }

        /// <summary>
        /// A header whose cells wrap ("Pediatric" over "Age Group") arrives as a sparse line above a fuller one; each
        /// such line is joined, cell by cell, into the line below it. A fuller line joins too when it fills columns the line below leaves empty (a staggered header). Stops (after
        /// at most five lines) at the first line that is at least three
        /// quarters full, or that is mostly numbers (data).
        /// </summary>
        private static List<(int Number, List<string> Cells)> MergeWrappedHeader(List<(int Number, List<string> Cells)> rows, int width)
        {
            int top = 0;
            while (top + 1 < rows.Count && top < 5)
            {
                var cells = rows[top].Cells;
                int filled = cells.Count(c => c.Length > 0);
                var next = rows[top + 1].Cells.Where(c => c.Length > 0).ToList();
                if (next.Count == 0 || next.Count(c => Numeric.IsMatch(c)) > next.Count / 2) break;
                // A staggered header (PXD037923: "Age | Height" over "Subject | Gender | weight") fills columns the line
                // above left empty; a data line repeats the header's columns.
                int shared = Enumerable.Range(0, width).Count(k => cells[k].Length > 0 && rows[top + 1].Cells[k].Length > 0);
                bool complements = shared <= 1 && next.Count < width;
                if (filled >= 0.75 * width && !complements) break;
                for (int k = 0; k < width; k++)
                    if (cells[k].Length > 0)
                        rows[top + 1].Cells[k] = rows[top + 1].Cells[k].Length > 0 ? cells[k] + " " + rows[top + 1].Cells[k] : cells[k];
                top++;
            }
            // A header's second half can sit BELOW its main line ("Age (years)" over "Age Group", or units: "Age" over
            // "(years)"): up to two sparse text lines under the header join it when each of their cells sits under a
            // header cell, or the line after them is fuller.
            for (int n = 0; n < 2 && top + 2 < rows.Count; n++)
            {
                var header = rows[top].Cells;
                var below = rows[top + 1].Cells;
                var filled = Enumerable.Range(0, width).Where(k => below[k].Length > 0).ToList();
                if (filled.Count == 0 || filled.Count >= 0.75 * width || filled.Any(k => Numeric.IsMatch(below[k]))) break;
                bool underHeader = filled.All(k => header[k].Length > 0);
                bool fullerAfter = rows[top + 2].Cells.Count(c => c.Length > 0) > filled.Count;
                if (!underHeader && !fullerAfter) break;
                foreach (int k in filled)
                    header[k] = header[k].Length > 0 ? header[k] + " " + below[k] : below[k];
                rows.RemoveAt(top + 1);
            }
            return rows.Skip(top).ToList();
        }

        /// <summary>A line's cells placed in the columns they overlap most; two cells in one column are joined.</summary>
        private static List<string> Assign(PdfLine line, List<(double Left, double Right)> columns)
        {
            var cells = Enumerable.Repeat("", columns.Count).ToList();
            foreach (var piece in line.Cells.SelectMany(c => Pieces(c, columns, line.Height)))
            {
                int best = Best(piece.Left, piece.Right, columns);
                cells[best] = cells[best].Length == 0 ? piece.Text : cells[best] + " " + piece.Text;
            }
            return cells;
        }

        /// <summary>
        /// A cell that crosses two or more columns is split at its wider internal gaps (over a third of a line height,
        /// wider than a word space) and each piece placed on its own: neighbouring sample IDs set 4 pt apart
        /// (PXD007160) read as one cell otherwise. A cell inside one column stays whole.
        /// </summary>
        private static IEnumerable<PdfCell> Pieces(PdfCell cell, List<(double Left, double Right)> columns, double height)
        {
            if (columns.Count(k => Math.Min(cell.Right, k.Right) - Math.Max(cell.Left, k.Left) > 0) < 2 || cell.Words.Count < 2)
            {
                yield return cell;
                yield break;
            }
            var words = new List<PdfWord> { cell.Words[0] };
            for (int i = 1; i < cell.Words.Count; i++)
            {
                var w = cell.Words[i];
                if (w.Left - words[^1].Right > height / 3 && Best(w.Left, w.Right, columns) != Best(words[0].Left, words[^1].Right, columns))
                {
                    yield return new PdfCell(words[0].Left, words.Max(x => x.Right), words);
                    words = new List<PdfWord>();
                }
                words.Add(w);
            }
            yield return new PdfCell(words[0].Left, words.Max(x => x.Right), words);
        }

        /// <summary>The column a span overlaps most (the nearest when it overlaps none).</summary>
        private static int Best(double left, double right, List<(double Left, double Right)> columns)
        {
            int best = 0;
            double bestOverlap = double.MinValue;
            for (int k = 0; k < columns.Count; k++)
            {
                double overlap = Math.Min(right, columns[k].Right) - Math.Max(left, columns[k].Left);
                if (overlap > bestOverlap) { bestOverlap = overlap; best = k; }
            }
            return best;
        }
    }
}
