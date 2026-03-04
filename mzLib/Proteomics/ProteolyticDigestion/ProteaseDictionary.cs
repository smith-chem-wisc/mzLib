using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;

namespace Proteomics.ProteolyticDigestion
{
    /// <summary>
    /// Provides a centralized dictionary of proteases used for protein digestion.
    /// 
    /// <para><b>Embedded Resource Architecture:</b></para>
    /// <para>
    /// This class loads its default protease definitions from an embedded resource (proteases.tsv) 
    /// compiled directly into the assembly. An embedded resource is a file that becomes part of the 
    /// compiled DLL/assembly at build time, rather than existing as a separate file on disk.
    /// </para>
    /// 
    /// <para><b>Benefits of Embedded Resources:</b></para>
    /// <list type="bullet">
    ///   <item><description><b>Deployment simplicity:</b> No need to distribute or manage separate data files; 
    ///   the protease definitions travel with the assembly.</description></item>
    ///   <item><description><b>Version consistency:</b> Protease definitions are always matched to the library 
    ///   version, preventing mismatches between code and data.</description></item>
    ///   <item><description><b>Path independence:</b> Works regardless of where the application is installed 
    ///   or the current working directory.</description></item>
    ///   <item><description><b>Tamper resistance:</b> Users cannot accidentally modify or delete the default 
    ///   protease definitions.</description></item>
    /// </list>
    /// 
    /// <para><b>Potential Limitations:</b></para>
    /// <list type="bullet">
    ///   <item><description><b>Requires rebuild to modify defaults:</b> Changing the default proteases requires 
    ///   recompiling the library.</description></item>
    ///   <item><description><b>Memory usage:</b> The resource is loaded into memory (though proteases.tsv is small, 
    ///   this is negligible).</description></item>
    /// </list>
    /// 
    /// <para><b>Custom Protease Support:</b></para>
    /// <para>
    /// While defaults come from the embedded resource, this class fully supports user-defined custom proteases 
    /// via <see cref="LoadAndMergeCustomProteases"/>. Users can:
    /// </para>
    /// <list type="bullet">
    ///   <item><description>Add new proteases not in the default set</description></item>
    ///   <item><description>Override built-in protease definitions with custom cleavage rules</description></item>
    ///   <item><description>Reset to defaults at any time via <see cref="ResetToDefaults"/></description></item>
    /// </list>
    /// <para>
    /// This design provides reliable defaults out-of-the-box while maintaining full flexibility for 
    /// specialized research needs.
    /// </para>
    /// </summary>
    public static class ProteaseDictionary
    {
        private const string EmbeddedProteaseResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";
        private const string EmbeddedProteaseModsResourceName = "Proteomics.ProteolyticDigestion.protease_mods.txt";

        static ProteaseDictionary()
        {
            // Load from embedded resources with embedded modifications
            // This ensures proteases with cleavage modifications (like CNBr) work out-of-the-box
            Dictionary = LoadProteaseDictionaryWithEmbeddedMods();
        }

        public static Dictionary<string, Protease> Dictionary { get; set; }

        #region Embedded Resource Loading

        /// <summary>
        /// Loads the protease dictionary using embedded resources for both proteases and their modifications.
        /// This is the recommended method for applications that want to use mzLib's default proteases
        /// without providing their own modification definitions.
        /// </summary>
        /// <returns>Dictionary of protease name to Protease object, with all embedded modifications applied.</returns>
        public static Dictionary<string, Protease> LoadProteaseDictionaryWithEmbeddedMods()
        {
            var embeddedMods = LoadEmbeddedProteaseMods();
            return LoadProteaseDictionaryFromEmbeddedResource(embeddedMods);
        }

        /// <summary>
        /// Loads protease-specific modifications from the embedded resource (protease_mods.txt).
        /// These are modifications that are applied at cleavage sites (e.g., Homoserine lactone for CNBr).
        /// </summary>
        /// <returns>
        /// List of modifications defined in the embedded protease_mods.txt resource.
        /// Returns an empty list if the resource is not found or cannot be parsed.
        /// </returns>
        public static List<Modification> LoadEmbeddedProteaseMods()
        {
            var assembly = typeof(ProteaseDictionary).Assembly;

            using (var stream = assembly.GetManifestResourceStream(EmbeddedProteaseModsResourceName))
            {
                if (stream == null)
                {
                    // Return empty list if the embedded mods file doesn't exist
                    // This allows backward compatibility with assemblies that don't have the file
                    return new List<Modification>();
                }

                using (var reader = new StreamReader(stream))
                {
                    string fileContent = reader.ReadToEnd();
                    return ParseModificationsFromString(fileContent);
                }
            }
        }

        /// <summary>
        /// Parses modifications from a string in the standard mzLib modification format.
        /// This is a simplified parser for the protease modifications embedded resource.
        /// </summary>
        /// <param name="content">The modification file content as a string.</param>
        /// <returns>List of parsed Modification objects.</returns>
        private static List<Modification> ParseModificationsFromString(string content)
        {
            var modifications = new List<Modification>();

            // Split by modification entry delimiter
            var entries = content.Split(new[] { "//" }, StringSplitOptions.RemoveEmptyEntries);

            foreach (var entry in entries)
            {
                var lines = entry.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

                string id = null;
                string target = null;
                string locationRestriction = null;
                string modificationType = null;
                Chemistry.ChemicalFormula chemicalFormula = null;
                double? monoisotopicMass = null;
                Dictionary<string, IList<string>> databaseReference = null;

                foreach (var line in lines)
                {
                    // Skip comments and empty lines
                    if (string.IsNullOrWhiteSpace(line) || line.TrimStart().StartsWith("#"))
                        continue;

                    if (line.Length < 5)
                        continue;

                    string key = line.Substring(0, 2).Trim();
                    string value = line.Length > 5 ? line.Substring(5).Trim() : string.Empty;

                    switch (key)
                    {
                        case "ID":
                            id = value;
                            break;
                        case "TG":
                            target = value.TrimEnd('.');
                            break;
                        case "PP":
                            locationRestriction = value;
                            break;
                        case "MT":
                            modificationType = value;
                            break;
                        case "CF":
                            try
                            {
                                chemicalFormula = Chemistry.ChemicalFormula.ParseFormula(value.Replace(" ", string.Empty));
                            }
                            catch
                            {
                                // Skip invalid formulas
                            }
                            break;
                        case "MM":
                            if (double.TryParse(value, System.Globalization.NumberStyles.Any,
                                System.Globalization.CultureInfo.InvariantCulture, out double mm))
                            {
                                monoisotopicMass = mm;
                            }
                            break;
                        case "DR":
                            var splitDR = value.TrimEnd('.').Split(new[] { "; " }, StringSplitOptions.None);
                            if (splitDR.Length >= 2)
                            {
                                databaseReference ??= new Dictionary<string, IList<string>>();
                                if (databaseReference.TryGetValue(splitDR[0], out var existingList))
                                {
                                    existingList.Add(splitDR[1]);
                                }
                                else
                                {
                                    databaseReference.Add(splitDR[0], new List<string> { splitDR[1] });
                                }
                            }
                            break;
                    }
                }

                // Create the modification if we have at least an ID
                if (!string.IsNullOrEmpty(id))
                {
                    ModificationMotif motif = null;
                    if (!string.IsNullOrEmpty(target))
                    {
                        ModificationMotif.TryGetMotif(target, out motif);
                    }

                    var mod = new Modification(
                        _originalId: id,
                        _accession: null,
                        _modificationType: modificationType ?? "Protease",
                        _featureType: null,
                        _target: motif,
                        _locationRestriction: locationRestriction,
                        _chemicalFormula: chemicalFormula,
                        _monoisotopicMass: monoisotopicMass,
                        _databaseReference: databaseReference,
                        _taxonomicRange: null,
                        _keywords: null,
                        _neutralLosses: null,
                        _diagnosticIons: null,
                        _fileOrigin: "Embedded:protease_mods.txt"
                    );

                    modifications.Add(mod);
                }
            }

            return modifications;
        }

        /// <summary>
        /// Loads the default proteases from the embedded resource.
        /// </summary>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        private static Dictionary<string, Protease> LoadProteaseDictionaryFromEmbeddedResource(List<Modification> proteaseMods)
        {
            var assembly = typeof(ProteaseDictionary).Assembly;

            using (var stream = assembly.GetManifestResourceStream(EmbeddedProteaseResourceName))
            {
                if (stream == null)
                {
                    throw new MzLibException(
                        $"Could not find embedded resource '{EmbeddedProteaseResourceName}'. " +
                        $"Available resources: {string.Join(", ", assembly.GetManifestResourceNames())}");
                }

                using (var reader = new StreamReader(stream))
                {
                    string fileContent = reader.ReadToEnd();
                    // RemoveEmptyEntries skips blank lines and lines with only whitespace,
                    // which is the desired behavior for TSV parsing
                    string[] lines = fileContent.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseProteaseLines(lines, proteaseMods);
                }
            }
        }

        #endregion

        #region External File Loading

        /// <summary>
        /// Loads proteases from an external file path. Useful for loading custom user-defined proteases.
        /// Returns a new dictionary containing only the proteases from the specified file.
        /// To merge with the main dictionary, use <see cref="LoadAndMergeCustomProteases"/>.
        /// </summary>
        /// <param name="path">Path to the proteases.tsv file.</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        public static Dictionary<string, Protease> LoadProteaseDictionary(string path, List<Modification> proteaseMods = null)
        {
            string[] myLines = File.ReadAllLines(path);
            return ParseProteaseLines(myLines, proteaseMods);
        }

        /// <summary>
        /// Loads custom proteases from a file and merges them into the main <see cref="Dictionary"/>.
        /// 
        /// Merge rules:
        /// - If a protease name already exists in the main dictionary, it will be overwritten with the custom definition
        /// - If a protease name is new, it will be added to the main dictionary
        /// 
        /// This allows users to:
        /// 1. Override built-in protease definitions with custom cleavage rules
        /// 2. Add entirely new proteases not included in the default set
        /// </summary>
        /// <param name="path">Path to the custom proteases.tsv file.</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>List of protease names that were added or updated.</returns>
        public static List<string> LoadAndMergeCustomProteases(string path, List<Modification> proteaseMods = null)
        {
            var customProteases = LoadProteaseDictionary(path, proteaseMods);
            var addedOrUpdated = new List<string>();

            foreach (var kvp in customProteases)
            {
                if (Dictionary.ContainsKey(kvp.Key))
                {
                    // Overwrite existing protease
                    Dictionary[kvp.Key] = kvp.Value;
                }
                else
                {
                    // Add new protease
                    Dictionary.Add(kvp.Key, kvp.Value);
                }
                addedOrUpdated.Add(kvp.Key);
            }

            return addedOrUpdated;
        }

        /// <summary>
        /// Resets the dictionary to the default embedded proteases, discarding any custom additions.
        /// Uses embedded modifications by default.
        /// </summary>
        /// <param name="proteaseMods">
        /// Optional list of modifications to apply to proteases that require them.
        /// If null, embedded modifications will be used.
        /// </param>
        public static void ResetToDefaults(List<Modification> proteaseMods = null)
        {
            if (proteaseMods == null)
            {
                // Use embedded modifications
                Dictionary = LoadProteaseDictionaryWithEmbeddedMods();
            }
            else
            {
                // Use provided modifications
                Dictionary = LoadProteaseDictionaryFromEmbeddedResource(proteaseMods);
            }
        }

        #endregion

        /// <summary>
        /// Gets a protease by name, with backward compatibility for old naming conventions.
        /// Old names like "chymotrypsin (don't cleave before proline)" are automatically 
        /// converted to the new format "chymotrypsin|P".
        /// </summary>
        /// <param name="name">The protease name (old or new format).</param>
        /// <returns>The Protease object.</returns>
        /// <exception cref="KeyNotFoundException">Thrown when the protease is not found after normalization.</exception>
        public static Protease GetProtease(string name)
        {
            // Try exact match first
            if (Dictionary.TryGetValue(name, out var protease))
            {
                return protease;
            }

            // Try normalizing old-style name
            string normalizedName = NormalizeProteaseName(name);
            if (normalizedName != name && Dictionary.TryGetValue(normalizedName, out protease))
            {
                return protease;
            }

            throw new KeyNotFoundException($"Protease '{name}' not found in dictionary. " +
                $"If using an old-style name, ensure it follows the pattern 'name (don't cleave before proline)' " +
                $"which maps to 'name|P'.");
        }

        /// <summary>
        /// Tries to get a protease by name, with backward compatibility for old naming conventions.
        /// Old names like "chymotrypsin (don't cleave before proline)" are automatically 
        /// converted to the new format "chymotrypsin|P".
        /// </summary>
        /// <param name="name">The protease name (old or new format).</param>
        /// <param name="protease">When successful, contains the Protease object; otherwise null.</param>
        /// <returns>True if the protease was found; otherwise false.</returns>
        public static bool TryGetProtease(string name, out Protease protease)
        {
            // Try exact match first
            if (Dictionary.TryGetValue(name, out protease))
            {
                return true;
            }

            // Try normalizing old-style name
            string normalizedName = NormalizeProteaseName(name);
            if (normalizedName != name && Dictionary.TryGetValue(normalizedName, out protease))
            {
                return true;
            }

            protease = null;
            return false;
        }

        /// <summary>
        /// Normalizes old-style protease names to the new format.
        /// <para>
        /// Converts patterns like:
        /// <list type="bullet">
        ///   <item><description>"chymotrypsin (don't cleave before proline)" → "chymotrypsin|P"</description></item>
        ///   <item><description>"Lys-C (don't cleave before proline)" → "Lys-C|P"</description></item>
        /// </list>
        /// </para>
        /// The "|P" suffix indicates the protease has a proline restriction (won't cleave when 
        /// the next residue is proline).
        /// </summary>
        /// <param name="name">The protease name to normalize.</param>
        /// <returns>The normalized protease name, or the original name if no pattern matched.</returns>
        public static string NormalizeProteaseName(string name)
        {
            if (string.IsNullOrEmpty(name))
                return name;

            // Common old-style patterns that indicate proline restriction (case-insensitive)
            // These all map to the "|P" suffix in the new naming convention
            string[] prolineRestrictionPatterns =
            {
                " (don't cleave before proline)",
            };

            foreach (var pattern in prolineRestrictionPatterns)
            {
                int index = name.IndexOf(pattern, StringComparison.OrdinalIgnoreCase);
                if (index >= 0)
                {
                    // Remove the pattern and trim, then add |P
                    string baseName = name.Substring(0, index).Trim();
                    return baseName + "|P";
                }
            }

            return name;
        }

        /// <summary>
        /// Parses protease definitions from TSV-formatted lines.
        /// Lines starting with '#' are treated as comments and skipped.
        /// The header line is parsed to determine column positions.
        /// Supports both old and new column name formats for backward compatibility.
        /// </summary>
        /// <param name="lines">Lines from the proteases file.</param>
        /// <param name="proteaseMods">Optional list of modifications to apply to proteases that require them.</param>
        /// <returns>Dictionary of protease name to Protease object.</returns>
        private static Dictionary<string, Protease> ParseProteaseLines(string[] lines, List<Modification> proteaseMods)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            Dictionary<string, int> columnIndices = null;

            // Column name aliases for backward compatibility (first name is the canonical name)
            // Old format used longer, more descriptive names; new format uses shorter names
            string[] nameAliases = { "name" };
            string[] motifAliases = { "motif", "sequences inducing cleavage" };
            string[] specificityAliases = { "specificity", "cleavage specificity" };
            string[] psiMsAccessionAliases = { "psi-ms accession", "psi-ms accession number" };
            string[] psiMsNameAliases = { "psi-ms name" };
            string[] cleavageModAliases = { "cleavage modification", "cleavage mass shifts" };

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

                // Calculate minimum required field count based on the highest column index of required columns
                int minRequiredColumns = nameAliases.Concat(motifAliases).Concat(specificityAliases)
                    .Select(alias => columnIndices.TryGetValue(alias, out int idx) ? idx + 1 : 0)
                    .Max();

                // Check if we have enough fields to read all required columns
                if (fields.Length < minRequiredColumns)
                {
                    throw new MzLibException(
                        $"Line for protease '{name}' has only {fields.Length} field(s), but the required columns " +
                        $"(Name, Motif, Specificity) extend to column {minRequiredColumns}. " +
                        $"Please ensure the line has enough tab-separated values.");
                }

                // Check if required fields have values
                if (string.IsNullOrWhiteSpace(specificityField))
                {
                    throw new MzLibException(
                        $"Line for protease '{name}' is missing a value for the required 'Specificity' column. " +
                        $"Please ensure all required fields (Name, Motif, Specificity) have values.");
                }

                string psiMsAccessionNumber = GetFieldValue(fields, columnIndices, psiMsAccessionAliases);
                string psiMsName = GetFieldValue(fields, columnIndices, psiMsNameAliases);
                string proteaseModDetails = GetFieldValue(fields, columnIndices, cleavageModAliases);

                List<DigestionMotif> motifList = DigestionMotif.ParseDigestionMotifsFromString(motifField);
                var cleavageSpecificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), specificityField, true);

                Protease protease = CreateProtease(
                    name, motifList, cleavageSpecificity, psiMsAccessionNumber, psiMsName,
                    proteaseModDetails, proteaseMods);

                if (dict.ContainsKey(protease.Name))
                {
                    throw new MzLibException($"More than one protease named {protease.Name} exists");
                }

                dict.Add(protease.Name, protease);
            }

            if (columnIndices == null)
            {
                throw new MzLibException("Protease file contains no header line. Expected columns: Name, Motif, Specificity");
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
                    $"Protease file header is missing required column '{displayName}'. " +
                    $"Expected one of: {string.Join(", ", aliases)}. " +
                    $"Found columns: {string.Join(", ", columnIndices.Keys)}");
            }
        }

        /// <summary>
        /// Gets a field value by trying multiple column name aliases, returning empty string if none exist.
        /// </summary>
        private static string GetFieldValue(string[] fields, Dictionary<string, int> columnIndices, string[] columnAliases)
        {
            // Return the value for the first alias that exists in columnIndices with a valid index
            foreach (string alias in columnAliases)
            {
                if (columnIndices.TryGetValue(alias, out int index) && index < fields.Length)
                {
                    return fields[index].Trim();
                }
            }
            return string.Empty;
        }

        /// <summary>
        /// Creates a Protease object, optionally with an associated modification.
        /// </summary>
        private static Protease CreateProtease(
            string name,
            List<DigestionMotif> motifList,
            CleavageSpecificity cleavageSpecificity,
            string psiMsAccessionNumber,
            string psiMsName,
            string proteaseModDetails,
            List<Modification> proteaseMods)
        {
            // If this protease has an associated modification, look it up
            if (!string.IsNullOrEmpty(proteaseModDetails) && proteaseMods != null)
            {
                Modification proteaseModification = proteaseMods
                    .FirstOrDefault(p => p.IdWithMotif == proteaseModDetails);

                if (proteaseModification != null)
                {
                    return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList, proteaseModification);
                }

                // Modification was specified but not found in the provided list
                throw new MzLibException($"{proteaseModDetails} is not a valid modification");
            }

            // No modification required
            return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
        }
    }
}