using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Reads psmtsv files into <see cref="LightWeightSpectralMatch"/> records.
    /// Supports two kinds of column-based filters:
    /// <list type="bullet">
    ///   <item><b>rowFilters</b> — rows that fail any filter are skipped, but reading continues.</item>
    ///   <item><b>terminatingFilters</b> — reading stops entirely on the first row that fails any filter.
    ///     The caller is responsible for ensuring the file is sorted by the relevant column.</item>
    /// </list>
    /// </summary>
    public static class LightWeightSpectralMatchTsvReader
    {
        private static readonly char[] Split = { '\t' };

        /// <summary>
        /// Reads a psmtsv file into lightweight spectral match records with optional filtering.
        /// </summary>
        /// <param name="filePath">Path to the psmtsv file.</param>
        /// <param name="warnings">Warnings generated during reading.</param>
        /// <param name="rowFilters">Optional filters that skip non-matching rows.
        ///     Key = column name from <see cref="SpectrumMatchFromTsvHeader"/>.
        ///     Value = function taking the raw cell string, returning true to include the row.</param>
        /// <param name="terminatingFilters">Optional filters that stop reading on first failure.
        ///     Same signature as rowFilters. Caller must ensure the file is sorted by the relevant column.</param>
        /// <returns>A list of <see cref="LightWeightSpectralMatch"/> records that passed all filters.</returns>
        public static List<LightWeightSpectralMatch> ReadTsv(
            string filePath,
            out List<string> warnings,
            Dictionary<string, Func<string, bool>>? rowFilters = null,
            Dictionary<string, Func<string, bool>>? terminatingFilters = null)
        {
            warnings = new List<string>();

            StreamReader reader;
            string headerLine;
            try
            {
                reader = new StreamReader(filePath);
                headerLine = reader.ReadLine() ?? throw new MzLibException("File is empty: " + filePath);
            }
            catch (MzLibException) { throw; }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            Dictionary<string, int> parsedHeader = ParseLightweightHeader(headerLine, rowFilters, terminatingFilters);

            // Resolve filter column indices once
            var resolvedRowFilters = ResolveFilters(rowFilters, parsedHeader);
            var resolvedTerminatingFilters = ResolveFilters(terminatingFilters, parsedHeader);

            var results = new List<LightWeightSpectralMatch>();
            int lineNumber = 1; // 1-based, after header
            bool terminated = false;

            try
            {
                string? line;
                while ((line = reader.ReadLine()) != null)
                {
                    lineNumber++;
                    string[] spl = line.Split(Split);
                    for (int j = 0; j < spl.Length; j++)
                        spl[j] = spl[j].Trim('"');

                    // Check terminating filters first
                    if (resolvedTerminatingFilters != null)
                    {
                        foreach (var (colIndex, filterFunc) in resolvedTerminatingFilters)
                        {
                            if (colIndex < 0 || colIndex >= spl.Length)
                                continue;
                            if (!filterFunc(spl[colIndex]))
                            {
                                terminated = true;
                                break;
                            }
                        }
                        if (terminated) break;
                    }

                    // Check row filters
                    if (resolvedRowFilters != null)
                    {
                        bool skipRow = false;
                        foreach (var (colIndex, filterFunc) in resolvedRowFilters)
                        {
                            if (colIndex < 0 || colIndex >= spl.Length)
                                continue;
                            if (!filterFunc(spl[colIndex]))
                            {
                                skipRow = true;
                                break;
                            }
                        }
                        if (skipRow) continue;
                    }

                    // Parse the line into a lightweight record
                    try
                    {
                        results.Add(new LightWeightSpectralMatch(spl, parsedHeader));
                    }
                    catch (Exception)
                    {
                        warnings.Add("Could not read line: " + lineNumber);
                    }
                }
            }
            finally
            {
                reader.Dispose();
            }

            return results;
        }

        /// <summary>
        /// Parses only the columns needed by <see cref="LightWeightSpectralMatch"/> plus any
        /// columns referenced by the provided filters. Much lighter than the full
        /// <see cref="SpectrumMatchTsvReader.ParseHeader"/> which maps 60+ columns.
        /// </summary>
        private static Dictionary<string, int> ParseLightweightHeader(
            string header,
            Dictionary<string, Func<string, bool>>? rowFilters,
            Dictionary<string, Func<string, bool>>? terminatingFilters)
        {
            var spl = header.Split(Split);

            // Build reverse lookup: column name → index
            var columnIndex = new Dictionary<string, int>(spl.Length);
            for (int i = 0; i < spl.Length; i++)
            {
                columnIndex[spl[i]] = i;
            }

            int IndexOf(string name) => columnIndex.TryGetValue(name, out int idx) ? idx : -1;

            var parsedHeader = new Dictionary<string, int>();

            // Required columns
            parsedHeader[SpectrumMatchFromTsvHeader.FileName] = IndexOf(SpectrumMatchFromTsvHeader.FileName);
            parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanNumber] = IndexOf(SpectrumMatchFromTsvHeader.Ms2ScanNumber);
            parsedHeader[SpectrumMatchFromTsvHeader.FullSequence] = IndexOf(SpectrumMatchFromTsvHeader.FullSequence);
            parsedHeader[SpectrumMatchFromTsvHeader.BaseSequence] = IndexOf(SpectrumMatchFromTsvHeader.BaseSequence);
            parsedHeader[SpectrumMatchFromTsvHeader.Score] = IndexOf(SpectrumMatchFromTsvHeader.Score);
            parsedHeader[SpectrumMatchFromTsvHeader.PrecursorCharge] = IndexOf(SpectrumMatchFromTsvHeader.PrecursorCharge);
            parsedHeader[SpectrumMatchFromTsvHeader.DecoyContaminantTarget] = IndexOf(SpectrumMatchFromTsvHeader.DecoyContaminantTarget);
            parsedHeader[SpectrumMatchFromTsvHeader.GeneName] = IndexOf(SpectrumMatchFromTsvHeader.GeneName);
            parsedHeader[SpectrumMatchFromTsvHeader.OrganismName] = IndexOf(SpectrumMatchFromTsvHeader.OrganismName);
            parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime] = IndexOf(SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime);
            parsedHeader[SpectrumMatchFromTsvHeader.QValue] = IndexOf(SpectrumMatchFromTsvHeader.QValue);
            parsedHeader[SpectrumMatchFromTsvHeader.PEP_QValue] = IndexOf(SpectrumMatchFromTsvHeader.PEP_QValue);

            // Accession — handle legacy column names
            if (columnIndex.ContainsKey(SpectrumMatchFromTsvHeader.Accession))
            {
                parsedHeader[SpectrumMatchFromTsvHeader.Accession] = IndexOf(SpectrumMatchFromTsvHeader.Accession);
                parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass] = IndexOf(SpectrumMatchFromTsvHeader.MonoisotopicMass);
            }
            else
            {
                parsedHeader[SpectrumMatchFromTsvHeader.Accession] = IndexOf(SpectrumMatchFromTsvHeader.ProteinAccession);
                parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass] = IndexOf(SpectrumMatchFromTsvHeader.PeptideMonoMass);
            }

            // TMT/Isobaric reporter ion channels
            foreach (var channelName in SpectrumMatchFromTsvHeader.TmtChannelNames)
            {
                parsedHeader[channelName] = IndexOf(channelName);
            }

            // Add any filter columns that aren't already in the header dict
            AddFilterColumns(rowFilters, columnIndex, parsedHeader);
            AddFilterColumns(terminatingFilters, columnIndex, parsedHeader);

            return parsedHeader;
        }

        private static void AddFilterColumns(
            Dictionary<string, Func<string, bool>>? filters,
            Dictionary<string, int> columnIndex,
            Dictionary<string, int> parsedHeader)
        {
            if (filters == null) return;
            foreach (var columnName in filters.Keys)
            {
                if (!parsedHeader.ContainsKey(columnName))
                {
                    parsedHeader[columnName] = columnIndex.TryGetValue(columnName, out int idx) ? idx : -1;
                }
            }
        }

        /// <summary>
        /// Converts a filter dictionary from column-name-keyed to column-index-keyed for fast per-line evaluation.
        /// </summary>
        private static List<(int colIndex, Func<string, bool> filter)>? ResolveFilters(
            Dictionary<string, Func<string, bool>>? filters,
            Dictionary<string, int> parsedHeader)
        {
            if (filters == null || filters.Count == 0) return null;

            var resolved = new List<(int, Func<string, bool>)>(filters.Count);
            foreach (var (columnName, filterFunc) in filters)
            {
                int colIndex = parsedHeader.TryGetValue(columnName, out int idx) ? idx : -1;
                resolved.Add((colIndex, filterFunc));
            }
            return resolved;
        }
    }
}
