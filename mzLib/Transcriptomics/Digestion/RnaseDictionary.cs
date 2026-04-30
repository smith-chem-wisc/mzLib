using MzLibUtil;
using Omics.Digestion;

namespace Transcriptomics.Digestion
{
    /// <summary>
    /// Provides a centralized, protected-baseline dictionary of RNases used for RNA digestion.
    ///
    /// <para><b>Embedded Resource Architecture:</b></para>
    /// <para>
    /// All default RNase definitions are loaded from the embedded resource (rnases.tsv) compiled
    /// directly into the assembly. When the library is updated, the embedded resource is updated
    /// with it — any local file previously named rnases.tsv is superseded automatically because
    /// the embedded resource is never read from disk.
    /// </para>
    ///
    /// <para><b>Custom Digestion Agent Support:</b></para>
    /// <para>
    /// Downstream consumers (e.g. MetaMorpheus, ProteaseGuru) may supplement the default set with
    /// custom RNases by calling <see cref="LoadAndMergeCustomRnases(IEnumerable{string})"/>.
    /// Names that collide with an embedded entry are not overridden by the merge; the embedded
    /// definition wins and the skipped name is reported via the returned
    /// <see cref="CustomDigestionAgentLoadResult"/> rather than throwing, so callers can warn users
    /// without crashing.
    /// </para>
    ///
    /// <para><b>Design notes:</b></para>
    /// <list type="bullet">
    ///   <item><description>No default custom file is provided by mzLib. Downstream projects supply their own.</description></item>
    ///   <item><description>The <see cref="Dictionary"/> property has a private setter; external code cannot
    ///   replace the whole dictionary, though the dictionary's contents can still be modified directly via the indexer.</description></item>
    ///   <item><description>Duplicate names within a single custom file still throw, because that file is itself
    ///   malformed and the caller should fix it.</description></item>
    /// </list>
    /// </summary>
    public static class RnaseDictionary
    {
        private const string EmbeddedResourceName = "Transcriptomics.Digestion.rnases.tsv";

        static RnaseDictionary()
        {
            Dictionary = LoadEmbeddedRnaseDictionary();
        }

        /// <summary>
        /// The active RNase dictionary. Seeded from the embedded resource at startup.
        /// Custom RNases may be added via <see cref="LoadAndMergeCustomRnases(IEnumerable{string})"/>,
        /// but embedded entries are never replaced.
        /// The setter is private to prevent external code from swapping the entire dictionary.
        /// </summary>
        public static Dictionary<string, Rnase> Dictionary { get; private set; }

        #region Embedded Resource Loading

        /// <summary>
        /// Reads and parses the embedded rnases.tsv resource.
        /// Called once by the static constructor; not intended for direct use by consumers.
        /// Throws <see cref="MzLibException"/> if the resource cannot be found.
        /// </summary>
        private static Dictionary<string, Rnase> LoadEmbeddedRnaseDictionary()
        {
            var assembly = typeof(RnaseDictionary).Assembly;

            using (var stream = assembly.GetManifestResourceStream(EmbeddedResourceName))
            {
                if (stream == null)
                {
                    throw new MzLibException(
                        $"Could not find embedded resource '{EmbeddedResourceName}'. " +
                        $"Available resources: {string.Join(", ", assembly.GetManifestResourceNames())}");
                }

                using (var reader = new StreamReader(stream))
                {
                    string[] lines = reader.ReadToEnd()
                        .Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseRnaseLines(lines);
                }
            }
        }

        #endregion

        #region Custom Digestion Agent Loading

        /// <summary>
        /// Convenience overload for loading a single custom RNase file.
        /// See <see cref="LoadAndMergeCustomRnases(IEnumerable{string})"/> for full details.
        /// </summary>
        /// <param name="path">Path to a custom RNases TSV file.</param>
        /// <returns>
        /// A <see cref="CustomDigestionAgentLoadResult"/> describing which RNases were added and which were skipped.
        /// </returns>
        public static CustomDigestionAgentLoadResult LoadAndMergeCustomRnases(string path)
        {
            ArgumentNullException.ThrowIfNull(path);
            return LoadAndMergeCustomRnases(new[] { path });
        }

        /// <summary>
        /// Loads custom RNases from one or more TSV files and merges new entries into <see cref="Dictionary"/>.
        ///
        /// <para><b>Atomicity guarantee:</b></para>
        /// <para>
        /// All files are parsed and validated before any entry is added to <see cref="Dictionary"/>.
        /// If any file cannot be read or contains malformed data, the dictionary remains completely unchanged.
        /// </para>
        ///
        /// <para><b>Merge rules:</b></para>
        /// <list type="bullet">
        ///   <item><description>
        ///     A custom RNase whose name already exists in the embedded resource is <b>skipped</b>.
        ///     The embedded definition always wins. The skipped name is recorded in
        ///     <see cref="CustomDigestionAgentLoadResult.Skipped"/>.
        ///   </description></item>
        ///   <item><description>
        ///     A custom RNase whose name was already added by an earlier file in the same call is
        ///     likewise skipped (first-encountered wins).
        ///   </description></item>
        ///   <item><description>
        ///     A custom RNase with a genuinely new name is added to the dictionary and recorded in
        ///     <see cref="CustomDigestionAgentLoadResult.Added"/>.
        ///   </description></item>
        ///   <item><description>
        ///     Duplicate names <b>within</b> a single file still throw <see cref="MzLibException"/> because
        ///     that file is malformed and the caller should fix it.
        ///   </description></item>
        /// </list>
        ///
        /// <para><b>Note on Skipped semantics:</b> The returned <see cref="CustomDigestionAgentLoadResult.Skipped"/>
        /// list does not distinguish between collisions with embedded entries, collisions with names added by
        /// an earlier file in the same batch, and collisions with names added by previous calls. Callers that
        /// need to differentiate these cases (for user-facing diagnostics) must compare against the embedded
        /// baseline themselves.</para>
        ///
        /// <para>Files are processed in the order supplied. No default custom file is provided by mzLib;
        /// downstream consumers supply their own paths.</para>
        /// </summary>
        /// <param name="paths">One or more paths to custom RNases TSV files.</param>
        /// <returns>
        /// A <see cref="CustomDigestionAgentLoadResult"/> with the names of added and skipped RNases.
        /// </returns>
        /// <exception cref="FileNotFoundException">Thrown if any path does not exist. Dictionary is not modified.</exception>
        /// <exception cref="MzLibException">Thrown if any file contains duplicate RNase names. Dictionary is not modified.</exception>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="paths"/> is <c>null</c> or contains a <c>null</c> element.</exception>
        public static CustomDigestionAgentLoadResult LoadAndMergeCustomRnases(IEnumerable<string> paths)
        {
            ArgumentNullException.ThrowIfNull(paths);

            var added = new List<string>();
            var skipped = new List<string>();

            // Stage 1: Parse ALL files into a temporary collection.
            // If any file is missing or malformed this will throw before Dictionary is touched.
            var parsedFiles = new List<Dictionary<string, Rnase>>();
            foreach (string path in paths)
            {
                ArgumentNullException.ThrowIfNull(path);
                string[] lines = File.ReadAllLines(path);
                parsedFiles.Add(ParseRnaseLines(lines));
            }

            // Stage 2: Merge into Dictionary only after every file parsed successfully.
            // Build a combined view that also deduplicates across files (first-encountered wins).
            var staged = new Dictionary<string, Rnase>();
            foreach (var customRnases in parsedFiles)
            {
                foreach (var kvp in customRnases)
                {
                    if (Dictionary.ContainsKey(kvp.Key) || staged.ContainsKey(kvp.Key))
                    {
                        skipped.Add(kvp.Key);
                    }
                    else
                    {
                        staged.Add(kvp.Key, kvp.Value);
                        added.Add(kvp.Key);
                    }
                }
            }

            // Commit all staged entries atomically.
            foreach (var kvp in staged)
            {
                Dictionary.Add(kvp.Key, kvp.Value);
            }

            return new CustomDigestionAgentLoadResult(added.AsReadOnly(), skipped.AsReadOnly());
        }

        #endregion

        #region Deprecated API (kept for backward compatibility)

        /// <summary>
        /// Returns a fresh copy of the embedded RNase dictionary.
        /// </summary>
        /// <remarks>
        /// Preserved as a thin wrapper for downstream consumers that previously called the
        /// no-args <c>LoadRnaseDictionary()</c>. New code should read <see cref="Dictionary"/>
        /// directly, which is seeded from the embedded resource at startup.
        /// </remarks>
        [Obsolete("Read RnaseDictionary.Dictionary instead. " +
                  "The dictionary is now seeded from the embedded resource at startup; " +
                  "this method returns a fresh parsed copy of the same embedded data. " +
                  "Scheduled for removal in a future major version.", error: false)]
        public static Dictionary<string, Rnase> LoadRnaseDictionary()
        {
            return LoadEmbeddedRnaseDictionary();
        }

        /// <summary>
        /// Parses a custom RNases TSV file and returns a fresh dictionary of the entries it
        /// contains. Does NOT mutate <see cref="Dictionary"/>.
        /// </summary>
        /// <remarks>
        /// Preserved as a thin wrapper for downstream consumers that previously called
        /// <c>LoadRnaseDictionary(path)</c> and assigned the result themselves. New code should
        /// call <see cref="LoadAndMergeCustomRnases(string)"/>, which merges atomically and reports
        /// collisions via <see cref="CustomDigestionAgentLoadResult"/>.
        /// </remarks>
        [Obsolete("Use LoadAndMergeCustomRnases(path) instead. " +
                  "This method returns a parsed-but-unmerged dictionary; the new API merges into " +
                  "the global Dictionary atomically and reports collisions via CustomDigestionAgentLoadResult. " +
                  "Scheduled for removal in a future major version.", error: false)]
        public static Dictionary<string, Rnase> LoadRnaseDictionary(string path)
        {
            string[] lines = File.ReadAllLines(path);
            return ParseRnaseLines(lines);
        }

        /// <summary>
        /// Rebuilds <see cref="Dictionary"/> from the embedded resource, discarding any custom
        /// RNases that were merged in via <see cref="LoadAndMergeCustomRnases(string)"/>.
        /// </summary>
        /// <remarks>
        /// Preserved for downstream consumers that previously called <c>ResetToDefaults()</c> to
        /// roll back custom additions. New code should track the names returned in
        /// <see cref="CustomDigestionAgentLoadResult.Added"/> and remove them individually.
        /// </remarks>
        [Obsolete("The dictionary is now seeded from the embedded resource at startup. " +
                  "To roll back custom additions, remove the names returned in " +
                  "CustomDigestionAgentLoadResult.Added. Scheduled for removal in a future major version.",
                  error: false)]
        public static void ResetToDefaults()
        {
            Dictionary = LoadEmbeddedRnaseDictionary();
        }

        #endregion

        #region TSV Parsing

        /// <summary>
        /// Parses RNase definitions from TSV-formatted lines.
        /// Comment lines (starting with '#') and blank lines are skipped.
        /// The first non-comment line is treated as the header and used to resolve column positions.
        /// Duplicate names within the same set of lines throw <see cref="MzLibException"/>.
        /// </summary>
        private static Dictionary<string, Rnase> ParseRnaseLines(string[] lines)
        {
            var dict = new Dictionary<string, Rnase>();
            Dictionary<string, int> columnIndices = null;

            // Column name aliases for backward compatibility (first alias is canonical)
            string[] nameAliases = { "name" };
            string[] motifAliases = { "motif", "sequences inducing cleavage" };
            string[] specificityAliases = { "specificity", "cleavage specificity" };

            foreach (string line in lines)
            {
                string trimmedLine = line.Trim().TrimStart('\uFEFF');

                if (string.IsNullOrWhiteSpace(trimmedLine) || trimmedLine.StartsWith("#"))
                    continue;

                string[] fields = trimmedLine.Split('\t');

                if (columnIndices == null)
                {
                    columnIndices = ParseHeaderLine(fields);
                    ValidateRequiredColumn(columnIndices, "Name", nameAliases);
                    ValidateRequiredColumn(columnIndices, "Motif", motifAliases);
                    ValidateRequiredColumn(columnIndices, "Specificity", specificityAliases);
                    continue;
                }

                string name = GetFieldValue(fields, columnIndices, nameAliases);
                string motifField = GetFieldValue(fields, columnIndices, motifAliases);
                string specificityField = GetFieldValue(fields, columnIndices, specificityAliases);

                if (string.IsNullOrWhiteSpace(name))
                    continue; // skip lines without a name

                var motifList = DigestionMotif.ParseDigestionMotifsFromString(motifField);

                CleavageSpecificity cleavageSpecificity = CleavageSpecificity.None;
                if (!string.IsNullOrWhiteSpace(specificityField))
                {
                    cleavageSpecificity = (CleavageSpecificity)Enum.Parse(
                        typeof(CleavageSpecificity), specificityField, true);
                }

                var rnase = new Rnase(name, cleavageSpecificity, motifList);

                if (dict.ContainsKey(rnase.Name))
                {
                    throw new MzLibException(
                        $"More than one RNase named '{rnase.Name}' exists within the same file. " +
                        $"Please ensure all RNase names are unique within a single file.");
                }

                dict.Add(rnase.Name, rnase);
            }

            if (columnIndices == null)
            {
                throw new MzLibException(
                    "RNase file contains no header line. Expected columns: Name, Motif, Specificity");
            }

            return dict;
        }

        private static Dictionary<string, int> ParseHeaderLine(string[] headerFields)
        {
            var columnIndices = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            for (int i = 0; i < headerFields.Length; i++)
            {
                string columnName = headerFields[i].Trim().ToLowerInvariant();
                if (!string.IsNullOrEmpty(columnName))
                    columnIndices[columnName] = i;
            }

            return columnIndices;
        }

        private static void ValidateRequiredColumn(
            Dictionary<string, int> columnIndices,
            string displayName,
            string[] aliases)
        {
            if (!aliases.Any(alias => columnIndices.ContainsKey(alias)))
            {
                throw new MzLibException(
                    $"RNase file header is missing required column '{displayName}'. " +
                    $"Expected one of: {string.Join(", ", aliases)}. " +
                    $"Found columns: {string.Join(", ", columnIndices.Keys)}");
            }
        }

        private static string GetFieldValue(
            string[] fields,
            Dictionary<string, int> columnIndices,
            string[] columnAliases)
        {
            foreach (string alias in columnAliases)
            {
                if (columnIndices.TryGetValue(alias, out int index) && index < fields.Length)
                    return fields[index].Trim();
            }
            return string.Empty;
        }

        #endregion
    }
}