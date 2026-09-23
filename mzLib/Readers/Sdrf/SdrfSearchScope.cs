using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// An SDRF restricted to one search's files, and what the restriction could not place.
    /// </summary>
    /// <param name="Document">Rows only for acquisitions the search read, each naming the searched file.</param>
    /// <param name="SearchedWithoutRow">Searched files the SDRF has no row for. Reported, never invented.</param>
    /// <param name="DroppedDataFiles">Acquired files the SDRF describes but this search did not read.</param>
    /// <param name="Ambiguous">Acquisitions two or more searched files claimed; the first, in name order, was kept.</param>
    internal sealed record SdrfSearchScoping(
        SdrfDocument Document,
        IReadOnlyList<string> SearchedWithoutRow,
        IReadOnlyList<string> DroppedDataFiles,
        IReadOnlyList<string> Ambiguous);

    /// <summary>
    /// Restricts an SDRF -- deposited, improved or drafted -- to the files ONE SEARCH read (sdrf D36): the
    /// SDRF that governs a quantification, and then accompanies its results, names exactly those files.
    ///
    /// <para><b>The join (MAP-12, note N8).</b> A search reads <c>X-calib.mzML</c> while the SDRF names the
    /// acquired <c>X.raw</c>. A searched file is joined to a row by STEM, case-insensitively, after removing
    /// the derivative suffixes MetaMorpheus writes (<c>-calib</c>, <c>-averaged</c>, in any order and
    /// number), so a converted <c>.mzML</c> finds its <c>.raw</c>. Where the caller KNOWS a searched file's
    /// acquisition -- MetaMorpheus can read it from the mzML's own <c>sourceFile</c> (D26) -- that name wins
    /// over the suffix rule. The searched name is written to <c>comment[searched data file]</c>, directly
    /// after <c>comment[data file]</c> as <see cref="SdrfBuilder"/> places it; a column an earlier search
    /// wrote is replaced, because it is this search's field.</para>
    ///
    /// <para><b>Nothing else changes.</b> Every kept row is carried cell for cell: a search of replicates 1
    /// and 3 says 1 and 3. Ranking replicates for quantification is the design reader's job (MAP-33), and
    /// this type must not do it quietly on the way.</para>
    ///
    /// <para>Pure. Internal (D19).</para>
    /// </summary>
    internal static class SdrfSearchScope
    {
        private const string DataFile = "comment[data file]";
        private const string SearchedDataFile = "comment[searched data file]";

        private static readonly Regex Derivative = new(@"(?:-calib|-averaged)+$", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <summary>
        /// Keeps the rows whose acquisition the search read. Throws on a null argument, on an empty search,
        /// or on an SDRF with no <c>comment[data file]</c> column.
        /// </summary>
        /// <param name="sdrf">The SDRF to restrict.</param>
        /// <param name="searchedFiles">The files the search read, as paths or names.</param>
        /// <param name="acquiredNameOf">
        /// Optional: searched file (as given in <paramref name="searchedFiles"/>) -> the acquired file it came
        /// from, when the caller knows it. Overrides the suffix rule for those files.
        /// </param>
        public static SdrfSearchScoping Restrict(SdrfDocument sdrf, IEnumerable<string> searchedFiles,
            IReadOnlyDictionary<string, string>? acquiredNameOf = null)
        {
            if (sdrf == null) throw new ArgumentNullException(nameof(sdrf));
            if (searchedFiles == null) throw new ArgumentNullException(nameof(searchedFiles));
            var searched = searchedFiles.Where(s => !string.IsNullOrWhiteSpace(s)).Distinct(StringComparer.OrdinalIgnoreCase).ToList();
            if (searched.Count == 0) throw new ArgumentException("A search reads at least one file.", nameof(searchedFiles));
            var header = sdrf.Header.ToList();
            int dataFile = header.IndexOf(DataFile);
            if (dataFile < 0)
                throw new ArgumentException("The SDRF has no comment[data file] column, so no row can be joined to a searched file.", nameof(sdrf));

            // acquired stem -> the searched file that read it
            var bySearched = searched.OrderBy(s => FileName(s), StringComparer.Ordinal)
                .GroupBy(s => AcquiredStem(s, acquiredNameOf), StringComparer.OrdinalIgnoreCase).ToList();
            var ambiguous = bySearched.Where(g => g.Count() > 1)
                .Select(g => $"{g.Key}: {string.Join(", ", g.Select(FileName))}").ToList();
            var readerOf = bySearched.ToDictionary(g => g.Key, g => FileName(g.First()), StringComparer.OrdinalIgnoreCase);

            // The searched-file column goes directly after comment[data file]; an earlier search's is replaced.
            var oldSearched = header.Select((h, i) => (h, i)).Where(x => x.h == SearchedDataFile).Select(x => x.i).ToHashSet();
            var newHeader = new List<string>();
            for (int i = 0; i < header.Count; i++)
            {
                if (oldSearched.Contains(i)) continue;
                newHeader.Add(header[i]);
                if (i == dataFile) newHeader.Add(SearchedDataFile);
            }
            var sdrfHeader = new SdrfHeader(newHeader);

            var kept = new List<SdrfRow>();
            var dropped = new List<string>();
            var joined = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (var row in sdrf.Results)
            {
                var cells = row.Cells.Concat(Enumerable.Repeat("", Math.Max(0, header.Count - row.Cells.Count))).ToList();
                string acquired = cells[dataFile];
                string stem = Derivative.Replace(SdrfFileNamePattern.Stem(acquired), "");
                if (!readerOf.TryGetValue(stem, out var reader))
                {
                    if (!dropped.Contains(acquired, StringComparer.OrdinalIgnoreCase)) dropped.Add(acquired);
                    continue;
                }
                joined.Add(stem);
                var outCells = new List<string>();
                for (int i = 0; i < header.Count; i++)
                {
                    if (oldSearched.Contains(i)) continue;
                    outCells.Add(cells[i]);
                    if (i == dataFile) outCells.Add(reader);
                }
                kept.Add(new SdrfRow(sdrfHeader, outCells));
            }

            var withoutRow = bySearched.Where(g => !joined.Contains(g.Key)).Select(g => FileName(g.First())).ToList();
            return new SdrfSearchScoping(new SdrfDocument(sdrfHeader, kept), withoutRow, dropped, ambiguous);
        }

        private static string AcquiredStem(string searched, IReadOnlyDictionary<string, string>? acquiredNameOf)
        {
            if (acquiredNameOf != null && acquiredNameOf.TryGetValue(searched, out var acquired) && !string.IsNullOrWhiteSpace(acquired))
                return SdrfFileNamePattern.Stem(FileName(acquired));
            return Derivative.Replace(SdrfFileNamePattern.Stem(FileName(searched)), "");
        }

        /// <summary>The file name of a path written with either separator, on any platform.</summary>
        private static string FileName(string path)
        {
            string p = path.Trim().TrimEnd('/', '\\');
            int cut = Math.Max(p.LastIndexOf('/'), p.LastIndexOf('\\'));
            return cut >= 0 ? p[(cut + 1)..] : p;
        }
    }
}
