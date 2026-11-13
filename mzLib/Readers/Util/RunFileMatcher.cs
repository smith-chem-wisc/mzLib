namespace Readers
{
    /// <summary>
    /// Pairs the run names a search tool wrote into its report with the spectra files a caller holds,
    /// for <see cref="IQuantifiableResultFile.FileNameToFilePath"/>. A report names a run by its spectra
    /// file without directory or extension, and often by a converted copy of the file the caller has
    /// ("sample" searched as "sample_centroid"). A wrong pairing makes FlashLFQ quantify one run's
    /// identifications against another run's spectra, so an unclear pairing is always left out rather than guessed.
    /// </summary>
    internal static class RunFileMatcher
    {
        /// <summary>
        /// Maps each run name to the one spectra file it names.
        /// <list type="number">
        /// <item>A file whose name, without directory or extension, equals the run name, ignoring case, is
        /// an exact match. A run with exactly one exact match takes it.</item>
        /// <item>A run with no exact match takes the file whose name differs from it only by a delimited
        /// suffix (see <see cref="IsSuffixExtensionOf"/>), when exactly one such file is the most specific.</item>
        /// <item>A run with two exact matches, or two equally specific suffix matches, is left out.</item>
        /// <item>A file claimed by two runs is taken from both.</item>
        /// </list>
        /// A run that is left out makes <c>MakeIdentifications</c> throw by name, which is easier to
        /// diagnose than quantifying against the wrong file.
        /// </summary>
        /// <param name="runNames">The report's run names, without directory or extension.</param>
        /// <param name="fullFilePaths">The caller's spectra files.</param>
        /// <returns>Run name to full file path, for the runs that pair with exactly one file.</returns>
        internal static Dictionary<string, string> Match(IEnumerable<string> runNames, IEnumerable<string> fullFilePaths)
        {
            // File-system identity is case-insensitive on Windows, so fold case here rather than let
            // "S1.raw" and "s1.raw" survive as two candidates that could match two runs and then slip
            // past the contested-set cleanup below, which groups on the same comparer.
            List<string> files = fullFilePaths.Distinct(StringComparer.OrdinalIgnoreCase).ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var runName in runNames.Distinct())
            {
                string? exactMatch = null;
                int exactMatchCount = 0;

                string? fuzzyMatch = null;
                int fuzzyShared = -1;
                bool fuzzyAmbiguous = false;

                foreach (var file in files)
                {
                    string fileNameWithoutExtension = StripPathAndExtension(file);

                    // An empty name is a prefix of everything and would match every run
                    if (fileNameWithoutExtension.Length == 0)
                        continue;

                    // An exact match is the best possible pairing, but two entries can strip to the
                    // same name (C:\batch1\QC.raw and C:\batch2\QC.mzML both strip to "QC"), so count
                    // them rather than take the first -- the contested-set cleanup below only catches
                    // one file claimed by two runs, never two files claimed by one run.
                    if (fileNameWithoutExtension.Equals(runName, StringComparison.OrdinalIgnoreCase))
                    {
                        exactMatch = file;
                        exactMatchCount++;
                        continue;
                    }

                    if (!IsSuffixExtensionOf(runName, fileNameWithoutExtension))
                        continue;

                    // Rank by how much of the run name the candidate reproduces, not by the
                    // candidate's own length. When the file name is a prefix of the run name
                    // ("sample" -> run "sample_rep1_centroid") a longer file shares more of the run
                    // and is the likelier pairing. But when the run name is a prefix of the file name
                    // ("sample_rep1" and "sample_rep10" both extend run "sample") every candidate
                    // reproduces the whole run name and none is more specific, so those must tie and
                    // be rejected rather than ranked by their trailing digits.
                    int sharedWithRun = Math.Min(fileNameWithoutExtension.Length, runName.Length);
                    if (sharedWithRun > fuzzyShared)
                    {
                        fuzzyMatch = file;
                        fuzzyShared = sharedWithRun;
                        fuzzyAmbiguous = false;
                    }
                    else if (sharedWithRun == fuzzyShared)
                    {
                        fuzzyAmbiguous = true;
                    }
                }

                // A single exact match is unambiguous and outranks any fuzzy match; two of them give
                // no reason to prefer either.
                if (exactMatchCount == 1)
                    allFiles.Add(runName, exactMatch!);
                else if (exactMatchCount == 0 && fuzzyMatch is not null && !fuzzyAmbiguous)
                    allFiles.Add(runName, fuzzyMatch);
            }

            // One spectra file cannot be the source of two runs. If it were claimed twice, FlashLFQ
            // would merge both runs' identifications into a single SpectraFileInfo and quantify them
            // as one file, so drop the whole contested set rather than pick a winner.
            foreach (var contested in allFiles.GroupBy(pair => pair.Value, StringComparer.OrdinalIgnoreCase).Where(group => group.Count() > 1).ToList())
                foreach (var pair in contested)
                    allFiles.Remove(pair.Key);

            return allFiles;
        }

        /// <summary>
        /// Strips the directory and the extension from a path that a search tool wrote.
        /// Path.GetFileNameWithoutExtension is not used, because tools commonly run on Linux and write
        /// '/'-delimited paths that .NET on Linux would not split, and a Windows path written on Windows
        /// must still split when it is read on Linux.
        /// </summary>
        internal static string StripPathAndExtension(string? path)
        {
            if (string.IsNullOrEmpty(path)) return "";

            int lastSeparator = path.LastIndexOfAny(new[] { '/', '\\' });
            string fileName = lastSeparator >= 0 ? path.Substring(lastSeparator + 1) : path;

            int lastDot = fileName.LastIndexOf('.');
            return lastDot > 0 ? fileName.Substring(0, lastDot) : fileName;
        }

        /// <summary>
        /// True when <paramref name="runName"/> is <paramref name="fileName"/> followed by a
        /// delimited suffix, or the other way round. This is the shape a report takes when the search
        /// ran on a converted copy of the caller's file ("sample" -> run "sample_centroid").
        /// <para>
        /// The suffix has to start at a delimiter. Plain containment would make the file "A1" match
        /// the run "A10", and 96-well plate naming produces exactly that collision.
        /// </para>
        /// </summary>
        private static bool IsSuffixExtensionOf(string runName, string fileName)
        {
            // Either side may carry the extra suffix: the caller may hold the converted file instead
            string longer = runName.Length > fileName.Length ? runName : fileName;
            string shorter = runName.Length > fileName.Length ? fileName : runName;

            // Equal-length names cannot differ by a suffix. Guarding here also keeps the boundary
            // index below in range: an ordinal StartsWith on a strictly longer string always leaves
            // at least one character at shorter.Length.
            if (longer.Length == shorter.Length)
                return false;

            // Ordinal, not culture-aware: file-system identity is ordinal, readers pool these strings
            // with StringComparer.Ordinal, and it is faster inside this O(runs x files) loop. A
            // linguistic comparison can also report a prefix match between strings of unequal
            // code-point length (ICU collapses ignorable code points such as the soft hyphen), which
            // would push longer[shorter.Length] past the end.
            if (!longer.StartsWith(shorter, StringComparison.OrdinalIgnoreCase))
                return false;

            char boundary = longer[shorter.Length];
            return boundary is '_' or '-' or '.' or ' ';
        }
    }
}
