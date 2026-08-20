using System.Text;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// An SDRF-Proteomics (Sample and Data Relationship Format) document -- the HUPO-PSI standard
    /// for describing the sample-to-data-file relationship of a proteomics experiment. One row per
    /// data file; columns carry sample characteristics, acquisition and search metadata, and the
    /// experimental factor values. Specification v1.1.0 (2026-01).
    ///
    /// This reader is deliberately hand-rolled rather than CsvHelper-driven, which departs from
    /// every other reader in this project. The house pattern maps FIXED columns onto a record type
    /// with [Name(...)] attributes, and SDRF cannot be expressed that way: its column set is
    /// open-ended, its names are data, and names repeat within a single file. The nearest existing
    /// precedent in mzLib for dynamic headers is SpectrumMatchFromTsv's hand-rolled parsedHeader.
    ///
    /// Three further properties of real files were measured across the 1,236-file curated corpus at
    /// bigbio/sdrf-annotated-datasets and each one is load-bearing here:
    ///
    ///   1. THERE IS NO QUOTING OR ESCAPING. SDRF is strictly tab-delimited; a cell cannot contain
    ///      a tab or a newline, so no escape mechanism is needed or defined. Two corpus files
    ///      (PXD026824, PXD028735) contain literal double-quote characters INSIDE cell values, e.g.
    ///      "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y;MT=Variable" with the quotes part of the data. A CSV
    ///      parser with default quoting would strip them and the file would not round-trip. Splitting
    ///      on '\t' with no escape handling is both simpler and more faithful.
    ///
    ///   2. LINE ENDINGS ARE PRESERVED, NOT CHOSEN. The document remembers the ending it was read
    ///      with and reproduces it, defaulting to LF in memory. It does NOT pin one.
    ///      An earlier version pinned CRLF, on a survey reporting "all 1,236 corpus files are
    ///      CRLF" -- which was an artefact of measuring a working tree on a machine with
    ///      core.autocrlf=true. Upstream the corpus is LF (`git ls-files --eol` reports i/lf for
    ///      every file; the repository has no .gitattributes). Pinning either ending makes this
    ///      class normalise the one thing it exists not to normalise, and makes byte-identical
    ///      round-tripping depend on the reader's git configuration instead of on the reader.
    ///      No corpus file carries a byte-order mark, so none is written.
    ///
    ///   3. ROWS MAY BE SHORT. See <see cref="SdrfRow.Cells"/>.
    ///
    /// Together these make a byte-identical round trip achievable for every file in the corpus,
    /// which is the regression test this type exists to pass. Reporting a file's defects is the
    /// separate concern of the validator: a reader that repaired what it read could not round-trip,
    /// and silently rewriting a curated annotation is worse than reporting it.
    /// </summary>
    public class SdrfDocument : ResultFile<SdrfRow>, IResultFile
    {
        /// <summary>
        /// SDRF has no escape mechanism, so a cell may never contain a tab or a newline.
        /// Enforced on write, where a violation would silently corrupt the column structure.
        /// </summary>
        private const string CellSeparator = "\t";

        /// <summary>
        /// The line ending this document was read with, reproduced verbatim on write.
        ///
        /// NOT pinned to a platform or to a survey of the corpus. An earlier version pinned CRLF on
        /// the strength of "all 1,236 curated files are CRLF" -- which was an artefact of measuring
        /// a working tree on a machine with core.autocrlf=true. The upstream bytes are LF
        /// (git ls-files --eol reports i/lf for every file in bigbio/sdrf-annotated-datasets, which
        /// has no .gitattributes). Pinning either ending makes the writer NORMALISE the very thing
        /// this class exists not to normalise, and makes a byte-identical round trip depend on the
        /// reader's git configuration rather than on the reader.
        /// </summary>
        private string _lineEnding = DefaultLineEnding;

        /// <summary>
        /// What a document built IN MEMORY is written with. CRLF, matching lab practice and the
        /// Windows/Excel workflow these files are edited in, so everything we author is uniform.
        ///
        /// This is a choice, not a constraint -- both endings are valid SDRF. It is deliberately
        /// separate from the round-trip behaviour above: a file we READ keeps its own endings
        /// whatever they are, so preserving other people's bytes and keeping our own corpus
        /// consistent are not in tension.
        /// </summary>
        public const string DefaultLineEnding = "\r\n";

        /// <summary>
        /// Whether the source file ended with a line break. Reproduced on write so that a file
        /// lacking one round-trips too. Every corpus file has one; a hand-edited file may not.
        /// </summary>
        private bool _endsWithNewline = true;

        private SdrfHeader? _header;

        public override SupportedFileType FileType => SupportedFileType.Sdrf;

        public override Software Software { get; set; }

        public SdrfDocument(string filePath) : base(filePath, Software.Unspecified) { }

        /// <summary>
        /// Constructor used to initialize from the factory method. Required: FileReader.ReadResultFile
        /// creates result files through Activator.CreateInstance.
        /// </summary>
        public SdrfDocument() : base() { }

        /// <summary>
        /// Builds a document in memory, for writing. FilePath is left empty, which is what stops the
        /// base class lazy-loading over the top of these rows.
        /// </summary>
        public SdrfDocument(SdrfHeader header, IEnumerable<SdrfRow> rows) : base()
        {
            _header = header ?? throw new ArgumentNullException(nameof(header));
            Results = rows?.ToList() ?? throw new ArgumentNullException(nameof(rows));
        }

        /// <summary>
        /// The ordered column names. Triggers a load for the same reason Results does, so callers
        /// can inspect the columns without first touching a row.
        /// </summary>
        public SdrfHeader Header
        {
            get
            {
                if (_header is null && File.Exists(FilePath))
                    LoadResults();
                return _header ?? new SdrfHeader(Array.Empty<string>());
            }
        }

        public override void LoadResults()
        {
            // Read the whole file as text rather than streaming lines: SDRF documents are small
            // (the largest in the corpus is 5,798 rows) and reading the raw text is what lets the
            // line-ending and trailing-newline handling below stay explicit instead of being
            // silently normalized by StreamReader.ReadLine.
            string text = File.ReadAllText(FilePath, Encoding.UTF8);

            _lineEnding = DetectLineEnding(text);
            _endsWithNewline = text.EndsWith("\n", StringComparison.Ordinal)
                               || text.EndsWith("\r", StringComparison.Ordinal);

            var lines = SplitLines(text);
            if (lines.Count == 0)
                throw new MzLibException($"SDRF file is empty: '{FilePath}'");

            _header = new SdrfHeader(lines[0].Split(CellSeparator));

            var rows = new List<SdrfRow>(lines.Count - 1);
            for (int i = 1; i < lines.Count; i++)
                rows.Add(new SdrfRow(_header, lines[i].Split(CellSeparator)));

            Results = rows;
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            var header = Header;

            // Validate EVERYTHING before opening the stream. File.Create truncates, so validating
            // inside the write loop meant a bad cell in row 500 left 499 rows already written over
            // whatever was there -- the guard destroyed the file it was refusing to write, which is
            // exactly what this type's stated ethic forbids. Worse, WriteResults(FilePath) truncated
            // the source and then tripped ResultFile's lazy load, which re-read the now-empty file
            // and threw "SDRF file is empty" with the original already gone.
            // Materialise the rows ONCE, before anything is truncated, and never touch Results
            // again. ResultFile.Results reloads whenever the list is EMPTY -- not merely unloaded --
            // so a header-only document re-reads from disk on every access. The previous version
            // read Results again after File.Create had truncated the file, which on an in-place
            // write re-read the emptied source and threw "SDRF file is empty" with the original
            // already gone. That is the exact failure the comment below claimed to have fixed.
            var rows = Results.ToList();

            RejectUnrepresentable(header, "column name", outputPath);
            foreach (var row in rows)
                RejectUnrepresentable(row.Cells, "value", outputPath);

            // UTF8Encoding(false) -- no byte-order mark. Encoding.UTF8 emits one, and not one of the
            // 1,236 corpus files has it; adding a BOM would break a byte-identical round trip on the
            // very first byte.
            using var writer = new StreamWriter(File.Create(outputPath), new UTF8Encoding(false));

            writer.Write(string.Join(CellSeparator, header));
            foreach (var row in rows)
            {
                writer.Write(_lineEnding);
                writer.Write(string.Join(CellSeparator, row.Cells));
            }
            if (_endsWithNewline)
                writer.Write(_lineEnding);
        }

        /// <summary>
        /// A tab, CR or LF inside a cell or a column name has no representation in this format --
        /// there is no quoting to fall back on -- so it would silently shift every following column.
        /// Fail loudly instead of writing a file that reads back wrong.
        ///
        /// The header is checked too. It was not, originally, and that was the more dangerous of the
        /// two omissions: column names in an authored document are built from run metadata, and a
        /// tab in one emits N+1 columns, leaving every row one narrower than the header and the
        /// whole document misaligned.
        /// </summary>
        private static void RejectUnrepresentable(IEnumerable<string> values, string what, string outputPath)
        {
            foreach (var value in values)
            {
                if (value is null) continue;
                if (value.IndexOf('\t') < 0 && value.IndexOf('\n') < 0 && value.IndexOf('\r') < 0)
                    continue;

                throw new MzLibException(
                    $"An SDRF {what} cannot contain a tab, carriage return or line feed; the format " +
                    $"defines no escape mechanism. Offending {what} for '{outputPath}': '{value}'.");
            }
        }

        /// <summary>
        /// The dominant line ending of the source text. CRLF wins if present at all, since a file
        /// with any CRLF is a CRLF file with, at most, a stray bare LF.
        /// </summary>
        private static string DetectLineEnding(string text)
        {
            int crlf = text.IndexOf("\r\n", StringComparison.Ordinal);
            if (crlf >= 0) return "\r\n";
            if (text.IndexOf('\n') >= 0) return "\n";
            if (text.IndexOf('\r') >= 0) return "\r";
            return DefaultLineEnding;
        }

        /// <summary>
        /// Splits on CRLF or LF, dropping the final empty fragment produced by a trailing newline.
        ///
        /// Blank lines are NOT skipped in the middle of a document -- a blank line is a row of one
        /// empty cell, and dropping it would lose a row on write. Only the single trailing newline
        /// that every well-formed file ends with is absorbed.
        /// </summary>
        private static List<string> SplitLines(string text)
        {
            var lines = new List<string>();
            if (text.Length == 0)
                return lines;

            int start = 0;
            for (int i = 0; i < text.Length; i++)
            {
                char c = text[i];
                if (c != '\n' && c != '\r')
                    continue;

                lines.Add(text.Substring(start, i - start));

                // A bare CR is a line break in its own right. Treating only LF as a break meant an
                // old-Mac or stray-CR file collapsed into ONE line: a two-line file read back as a
                // single column named "a\rb" with zero rows, and the validator then reported
                // "NoRows" -- a confident and completely wrong diagnosis of a file full of data.
                if (c == '\r' && i + 1 < text.Length && text[i + 1] == '\n')
                    i++;

                start = i + 1;
            }

            // A file not ending in a line break still has a final line.
            if (start < text.Length)
                lines.Add(text.Substring(start));

            return lines;
        }
    }
}
