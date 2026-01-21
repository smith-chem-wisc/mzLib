using Easy.Common.Extensions;
using MzLibUtil;
using Omics.Digestion;

namespace Transcriptomics.Digestion
{
    /// <summary>
    /// Provides a centralized dictionary of RNases used for RNA digestion.
    /// 
    /// <para><b>Embedded Resource Architecture:</b></para>
    /// <para>
    /// This class loads its default RNase definitions from an embedded resource (rnases.tsv) 
    /// compiled directly into the assembly. An embedded resource is a file that becomes part of the 
    /// compiled DLL/assembly at build time, rather than existing as a separate file on disk.
    /// </para>
    /// 
    /// <para><b>Benefits of Embedded Resources:</b></para>
    /// <list type="bullet">
    ///   <item><description><b>Deployment simplicity:</b> No need to distribute or manage separate data files; 
    ///   the RNase definitions travel with the assembly.</description></item>
    ///   <item><description><b>Version consistency:</b> RNase definitions are always matched to the library 
    ///   version, preventing mismatches between code and data.</description></item>
    ///   <item><description><b>Path independence:</b> Works regardless of where the application is installed 
    ///   or the current working directory.</description></item>
    ///   <item><description><b>Tamper resistance:</b> Users cannot accidentally modify or delete the default 
    ///   RNase definitions.</description></item>
    /// </list>
    /// 
    /// <para><b>Potential Limitations:</b></para>
    /// <list type="bullet">
    ///   <item><description><b>Requires rebuild to modify defaults:</b> Changing the default RNases requires 
    ///   recompiling the library.</description></item>
    ///   <item><description><b>Memory usage:</b> The resource is loaded into memory (though rnases.tsv is small, 
    ///   this is negligible).</description></item>
    /// </list>
    /// 
    /// <para><b>Custom RNase Support:</b></para>
    /// <para>
    /// While defaults come from the embedded resource, this class fully supports user-defined custom RNases 
    /// via <see cref="LoadAndMergeCustomRnases"/>. Users can:
    /// </para>
    /// <list type="bullet">
    ///   <item><description>Add new RNases not in the default set</description></item>
    ///   <item><description>Override built-in RNase definitions with custom cleavage rules</description></item>
    ///   <item><description>Reset to defaults at any time via <see cref="ResetToDefaults"/></description></item>
    /// </list>
    /// <para>
    /// This design provides reliable defaults out-of-the-box while maintaining full flexibility for 
    /// specialized research needs.
    /// </para>
    /// </summary>
    public static class RnaseDictionary
    {
        private const string EmbeddedResourceName = "Transcriptomics.Digestion.rnases.tsv";

        static RnaseDictionary()
        {
            // Load from embedded resource
            Dictionary = LoadRnaseDictionary();
        }

        public static Dictionary<string, Rnase> Dictionary { get; set; }

        /// <summary>
        /// Loads the default RNases from the embedded resource.
        /// </summary>
        /// <returns>Dictionary of RNase name to Rnase object.</returns>
        public static Dictionary<string, Rnase> LoadRnaseDictionary()
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
                    string fileContent = reader.ReadToEnd();
                    // RemoveEmptyEntries skips blank lines and lines with only whitespace,
                    // which is the desired behavior for TSV parsing
                    string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseRnaseLines(lines);
                }
            }
        }

        /// <summary>
        /// Loads RNases from an external file path. Useful for loading custom user-defined RNases.
        /// Returns a new dictionary containing only the RNases from the specified file.
        /// To merge with the main dictionary, use <see cref="LoadAndMergeCustomRnases"/>.
        /// </summary>
        /// <param name="path">Path to the rnases.tsv file.</param>
        /// <returns>Dictionary of RNase name to Rnase object.</returns>
        public static Dictionary<string, Rnase> LoadRnaseDictionary(string path)
        {
            string[] myLines = File.ReadAllLines(path);
            return ParseRnaseLines(myLines);
        }

        /// <summary>
        /// Loads custom RNases from a file and merges them into the main <see cref="Dictionary"/>.
        /// 
        /// Merge rules:
        /// - If an RNase name already exists in the main dictionary, it will be overwritten with the custom definition
        /// - If an RNase name is new, it will be added to the main dictionary
        /// 
        /// This allows users to:
        /// 1. Override built-in RNase definitions with custom cleavage rules
        /// 2. Add entirely new RNases not included in the default set
        /// </summary>
        /// <param name="path">Path to the custom rnases.tsv file.</param>
        /// <returns>List of RNase names that were added or updated.</returns>
        public static List<string> LoadAndMergeCustomRnases(string path)
        {
            var customRnases = LoadRnaseDictionary(path);
            var addedOrUpdated = new List<string>();

            foreach (var kvp in customRnases)
            {
                if (Dictionary.ContainsKey(kvp.Key))
                {
                    // Overwrite existing RNase
                    Dictionary[kvp.Key] = kvp.Value;
                }
                else
                {
                    // Add new RNase
                    Dictionary.Add(kvp.Key, kvp.Value);
                }
                addedOrUpdated.Add(kvp.Key);
            }

            return addedOrUpdated;
        }

        /// <summary>
        /// Resets the dictionary to the default embedded RNases, discarding any custom additions.
        /// </summary>
        public static void ResetToDefaults()
        {
            Dictionary = LoadRnaseDictionary();
        }

        /// <summary>
        /// Parses RNase definitions from TSV-formatted lines.
        /// Lines starting with '#' are treated as comments and skipped.
        /// The header line is parsed to determine column positions.
        /// Supports both old and new column name formats for backward compatibility.
        /// </summary>
        /// <param name="lines">Lines from the RNases file.</param>
        /// <returns>Dictionary of RNase name to Rnase object.</returns>
        private static Dictionary<string, Rnase> ParseRnaseLines(string[] lines)
        {
            Dictionary<string, Rnase> dict = new Dictionary<string, Rnase>();
            Dictionary<string, int> columnIndices = null;

            // Column name aliases for backward compatibility (first name is the canonical name)
            // Old format used longer, more descriptive names; new format uses shorter names
            string[] nameAliases = { "name" };
            string[] motifAliases = { "motif", "sequences inducing cleavage" };
            string[] specificityAliases = { "specificity", "cleavage specificity" };

            foreach (string line in lines)
            {
                // Trim to handle potential BOM or leading/trailing whitespace
                string trimmedLine = line.Trim().TrimStart('\uFEFF'); // \uFEFF is the UTF-8 BOM character

                // Skip empty lines and comment lines
                if (string.IsNullOrWhiteSpace(trimmedLine) || trimmedLine.StartsWith("#"))
                {
                    continue;
                }

                string[] fields = trimmedLine.Split('\t');

                // Check if this is the header line (first non-comment, non-empty line should be header)
                if (columnIndices == null)
                {
                    columnIndices = ParseHeaderLine(fields);

                    // Validate required columns exist (using any of their aliases)
                    ValidateRequiredColumn(columnIndices, "Name", nameAliases);
                    ValidateRequiredColumn(columnIndices, "Motif", motifAliases);
                    ValidateRequiredColumn(columnIndices, "Specificity", specificityAliases);
                    continue;
                }

                // Parse data line using column indices from header
                string name = GetFieldValue(fields, columnIndices, nameAliases);
                string motifField = GetFieldValue(fields, columnIndices, motifAliases);
                string specificityField = GetFieldValue(fields, columnIndices, specificityAliases);

                // Validate that required fields are present
                if (string.IsNullOrWhiteSpace(name))
                {
                    continue; // Skip lines without a name
                }

                List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(motifField);
                
                // Default to 'none' specificity if not specified (for top-down)
                CleavageSpecificity cleavageSpecificity = CleavageSpecificity.None;
                if (!string.IsNullOrWhiteSpace(specificityField))
                {
                    cleavageSpecificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), specificityField, true);
                }

                var rnase = new Rnase(name, cleavageSpecificity, motifList);

                if (dict.ContainsKey(rnase.Name))
                {
                    throw new MzLibException($"More than one RNase named {rnase.Name} exists");
                }

                dict.Add(rnase.Name, rnase);
            }

            if (columnIndices == null)
            {
                throw new MzLibException("RNase file contains no header line. Expected columns: Name, Motif, Specificity");
            }

            return dict;
        }

        /// <summary>
        /// Parses the header line and returns a dictionary mapping column names (lowercase) to their indices.
        /// </summary>
        private static Dictionary<string, int> ParseHeaderLine(string[] headerFields)
        {
            var columnIndices = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            for (int i = 0; i < headerFields.Length; i++)
            {
                string columnName = headerFields[i].Trim().ToLowerInvariant();
                if (!string.IsNullOrEmpty(columnName))
                {
                    columnIndices[columnName] = i;
                }
            }

            return columnIndices;
        }

        /// <summary>
        /// Validates that at least one of the column name aliases exists in the header.
        /// </summary>
        private static void ValidateRequiredColumn(Dictionary<string, int> columnIndices, string displayName, string[] aliases)
        {
            if (!aliases.Any(alias => columnIndices.ContainsKey(alias)))
            {
                throw new MzLibException(
                    $"RNase file header is missing required column '{displayName}'. " +
                    $"Expected one of: {string.Join(", ", aliases)}. " +
                    $"Found columns: {string.Join(", ", columnIndices.Keys)}");
            }
        }

        /// <summary>
        /// Gets a field value by trying multiple column name aliases, returning empty string if none exist.
        /// </summary>
        private static string GetFieldValue(string[] fields, Dictionary<string, int> columnIndices, string[] columnAliases)
        {
            foreach (string alias in columnAliases)
            {
                if (columnIndices.TryGetValue(alias, out int index) && index < fields.Length)
                {
                    return fields[index].Trim();
                }
            }
            return string.Empty;
        }
    }
}