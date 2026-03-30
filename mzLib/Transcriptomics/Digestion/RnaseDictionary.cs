using MzLibUtil;
using Omics.Digestion;

namespace Transcriptomics.Digestion
{
    /// <summary>
    /// Provides a centralized, immutable-baseline dictionary of RNases used for RNA digestion.
    ///
    /// <para><b>Embedded Resource Architecture:</b></para>
    /// <para>
    /// All default RNase definitions are loaded from the embedded resource (rnases.tsv) compiled
    /// directly into the assembly. These definitions are the authoritative source of truth and cannot
    /// be overridden by any custom file. When the library is updated, the embedded resource is updated
    /// with it — any local file previously named rnases.tsv is superseded automatically because the
    /// embedded resource is never read from disk.
    /// </para>
    ///
    /// <para><b>Custom Digestion Agent Support:</b></para>
    /// <para>
    /// Downstream consumers (e.g. MetaMorpheus, ProteaseGuru) may supplement the default set with
    /// custom RNases by calling <see cref="LoadAndMergeCustomRnases(IEnumerable{string})"/>.
    /// Custom RNases are additive only — they can never replace an embedded definition.
    /// Name collisions are reported via the returned <see cref="CustomDigestionAgentLoadResult"/>
    /// rather than throwing, so callers can warn users without crashing.
    /// </para>
    ///
    /// <para><b>Design notes:</b></para>
    /// <list type="bullet">
    ///   <item><description>No default custom file is provided by mzLib. Downstream projects supply their own.</description></item>
    ///   <item><description>The <see cref="Dictionary"/> property has a private setter; external code cannot
    ///   replace the whole dictionary.</description></item>
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
            return LoadAndMergeCustomRnases(new[] { path });
        }

        /// <summary>
        /// Loads custom RNases from one or more TSV files and merges new entries into <see cref="Dictionary"/>.
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
        /// <para>Files are processed in the order supplied. No default custom file is provided by mzLib;
        /// downstream consumers supply their own paths.</para>
        /// </summary>
        /// <param name="paths">One or more paths to custom RNases TSV files.</param>
        /// <returns>
        /// A <see cref="CustomDigestionAgentLoadResult"/> with the names of added and skipped RNases.
        /// </returns>
        public static CustomDigestionAgentLoadResult LoadAndMergeCustomRnases(IEnumerable<string> paths)
        {
            var added = new List<string>();
            var skipped = new List<string>();

            foreach (string path in paths)
            {
                string[] lines = File.ReadAllLines(path);
                var customRnases = ParseRnaseLines(lines);

                foreach (var kvp in customRnases)
                {
                    if (Dictionary.ContainsKey(kvp.Key))
                    {
                        // Name collision — embedded or earlier custom entry wins
                        skipped.Add(kvp.Key);
                    }
                    else
                    {
                        Dictionary.Add(kvp.Key, kvp.Value);
                        added.Add(kvp.Key);
                    }
                }
            }

            return new CustomDigestionAgentLoadResult(added, skipped);
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
            string[] nameAliases        = { "name" };
            string[] motifAliases       = { "motif", "sequences inducing cleavage" };
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

                string name             = GetFieldValue(fields, columnIndices, nameAliases);
                string motifField       = GetFieldValue(fields, columnIndices, motifAliases);
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