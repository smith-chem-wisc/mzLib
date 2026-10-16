namespace Readers
{
    /// <summary>
    /// The sample columns of one input SDRF row, carried into <see cref="SdrfBuilder"/> output cell for cell.
    ///
    /// <para><b>Why a copy and not the per-field inputs.</b> <see cref="SdrfSample"/>'s fields cannot carry somebody
    /// else's row unchanged: <see cref="SdrfSample.Organism"/> is a <see cref="CvParam"/>, and a round trip through
    /// <see cref="SdrfCell.TryParseTerm"/> drops keys it does not know and reorders the rest;
    /// <see cref="SdrfSample.BiologicalReplicate"/> is a number, so <c>pooled</c> cannot be said; the dictionaries hold
    /// one value per column, so a repeated column loses all but its first; and the builder sorts columns, so the
    /// input's order is lost. A search that re-analyses a curated SDRF must hand its statistics the curator's sample
    /// facts, not a projection of them (QuantProject 036/037, SDRF-D1). So the cells are never parsed at all.</para>
    ///
    /// <para><b>What is a sample column.</b> <c>source name</c>; every <c>characteristics[...]</c> column, organism
    /// and biological replicate included; <c>material type</c>; every <c>factor value[...]</c> column; and every
    /// <c>comment[...]</c> column that is not one the builder writes itself
    /// (<see cref="SdrfBuilder.BuiltInCommentColumns"/>), such as <c>comment[sample preparation batch]</c>,
    /// <c>comment[file uri]</c> or the provenance columns. <c>assay name</c>, <c>technology type</c> and the builder's
    /// own comments are the assay half, which the search writes.</para>
    ///
    /// <para><b>What is not carried.</b> Two kinds of input column reach the output only if the caller sees to it.
    /// <c>comment[proteomexchange accession number]</c> is one of the builder's own columns, written only when
    /// <see cref="SdrfBuilderOptions.ProteomeXchangeAccession"/> is set, so a caller carrying a deposited SDRF passes
    /// its accession there. And a column that is none of the above is not copied: a misspelling
    /// (<c>characteristic[mouse line]</c>), a header with no name, or a case variant (<c>Material Type</c>, 6 of the
    /// 12,298 curated SDRFs), because column names are ordinal throughout (<see cref="SdrfHeader"/>). The input's
    /// header minus the built document's header names every such column, whatever the options, so a caller that
    /// must report them compares the two.</para>
    ///
    /// <para>Public for one named cross-assembly caller (D19): MetaMorpheus's SDRF writer, which copies the input
    /// SDRF's sample columns into the SDRF a search writes.</para>
    /// </summary>
    public sealed class SdrfSampleCopy
    {
        private const string SourceNameColumn = "source name";
        internal const string MaterialTypeColumn = "material type";

        private SdrfSampleCopy(IReadOnlyList<KeyValuePair<string, string>> cells, string sourceName)
        {
            Cells = cells;
            SourceName = sourceName;
        }

        /// <summary>
        /// Every sample column of the row, in the input header's order, repeats kept, values byte for byte. Column
        /// names are as the header spells them.
        /// </summary>
        public IReadOnlyList<KeyValuePair<string, string>> Cells { get; }

        /// <summary>The copied <c>source name</c>, which <see cref="SdrfSample.SourceName"/> must equal.</summary>
        public string SourceName { get; }

        /// <summary>
        /// Takes the sample columns of one row of <paramref name="document"/>.
        ///
        /// A row shorter than the header (PXD059974 has 17 such rows) gives <c>""</c> for the cells it does not
        /// reach. That is not invented here into a reserved word. The builder refuses a blank sample fact
        /// (<c>source name</c>, <c>characteristics[...]</c>, <c>material type</c>, <c>factor value[...]</c>), naming
        /// the column, because a blank there is not a statement anybody made; a blank <c>comment[...]</c> is written
        /// blank, as the input has it.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either argument is null.</exception>
        /// <exception cref="ArgumentException">The row is positioned against a different header, or the header has
        /// no <c>source name</c> column or more than one.</exception>
        public static SdrfSampleCopy FromRow(SdrfDocument document, SdrfRow row)
        {
            if (document is null) throw new ArgumentNullException(nameof(document));
            if (row is null) throw new ArgumentNullException(nameof(row));

            var header = document.Header;
            if (!row.Header.SequenceEqual(header, StringComparer.Ordinal))
                throw new ArgumentException("The row is positioned against a different header from the document's.", nameof(row));

            int sourceNames = header.Count(h => string.Equals(h, SourceNameColumn, StringComparison.Ordinal));
            if (sourceNames != 1)
                throw new ArgumentException(
                    $"An SDRF has exactly one '{SourceNameColumn}' column; this one has {sourceNames}.", nameof(document));

            var cells = new List<KeyValuePair<string, string>>();
            string sourceName = "";
            for (int i = 0; i < header.Count; i++)
            {
                string column = header[i];
                if (!IsSampleColumn(column)) continue;
                string value = i < row.Cells.Count ? row.Cells[i] : "";
                cells.Add(new KeyValuePair<string, string>(column, value));
                if (string.Equals(column, SourceNameColumn, StringComparison.Ordinal)) sourceName = value;
            }
            return new SdrfSampleCopy(cells, sourceName);
        }

        /// <summary>Ordinal, as column names are throughout (<see cref="SdrfHeader"/>).</summary>
        internal static bool IsSampleColumn(string column) =>
            string.Equals(column, SourceNameColumn, StringComparison.Ordinal)
            || column.StartsWith("characteristics[", StringComparison.Ordinal)
            || string.Equals(column, MaterialTypeColumn, StringComparison.Ordinal)
            || column.StartsWith("factor value[", StringComparison.Ordinal)
            || (column.StartsWith("comment[", StringComparison.Ordinal) && !SdrfBuilder.BuiltInCommentColumns.Contains(column));
    }
}
