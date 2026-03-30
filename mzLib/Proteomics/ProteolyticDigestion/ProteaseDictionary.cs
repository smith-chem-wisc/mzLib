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
    /// Provides a centralized, immutable-baseline dictionary of proteases used for protein digestion.
    ///
    /// <para><b>Embedded Resource Architecture:</b></para>
    /// <para>
    /// All default protease definitions are loaded from embedded resources (proteases.tsv and
    /// protease_mods.txt) compiled directly into the assembly. These definitions are the authoritative
    /// source of truth and cannot be overridden by any custom file. When the library is updated, the
    /// embedded resources are updated with it — any local file previously named proteases.tsv is
    /// superseded automatically because the embedded resource is never read from disk.
    /// </para>
    ///
    /// <para><b>Custom Digestion Agent Support:</b></para>
    /// <para>
    /// Downstream consumers (e.g. MetaMorpheus, ProteaseGuru) may supplement the default set with
    /// custom proteases by calling <see cref="LoadAndMergeCustomProteases(IEnumerable{string}, List{Modification})"/>.
    /// Custom proteases are additive only — they can never replace an embedded definition.
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
    public static class ProteaseDictionary
    {
        private const string EmbeddedProteaseResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";
        private const string EmbeddedProteaseModsResourceName = "Proteomics.ProteolyticDigestion.protease_mods.txt";

        static ProteaseDictionary()
        {
            Dictionary = LoadProteaseDictionaryWithEmbeddedMods();
        }

        /// <summary>
        /// The active protease dictionary. Seeded from the embedded resource at startup.
        /// Custom proteases may be added via <see cref="LoadAndMergeCustomProteases(IEnumerable{string}, List{Modification})"/>,
        /// but embedded entries are never replaced.
        /// The setter is private to prevent external code from swapping the entire dictionary.
        /// </summary>
        public static Dictionary<string, Protease> Dictionary { get; private set; }

        #region Embedded Resource Loading

        /// <summary>
        /// Loads the protease dictionary using both embedded resources (proteases.tsv and protease_mods.txt).
        /// Called once by the static constructor; not intended for direct use by consumers.
        /// </summary>
        private static Dictionary<string, Protease> LoadProteaseDictionaryWithEmbeddedMods()
        {
            var embeddedMods = LoadEmbeddedProteaseMods();
            return LoadProteaseDictionaryFromEmbeddedResource(embeddedMods);
        }

        /// <summary>
        /// Loads protease-specific modifications from the embedded protease_mods.txt resource.
        /// These are modifications applied at cleavage sites (e.g., Homoserine lactone for CNBr).
        /// Returns an empty list if the resource is absent, allowing backward compatibility.
        /// </summary>
        public static List<Modification> LoadEmbeddedProteaseMods()
        {
            var assembly = typeof(ProteaseDictionary).Assembly;

            using (var stream = assembly.GetManifestResourceStream(EmbeddedProteaseModsResourceName))
            {
                if (stream == null)
                    return new List<Modification>();

                using (var reader = new StreamReader(stream))
                {
                    return ParseModificationsFromString(reader.ReadToEnd());
                }
            }
        }

        /// <summary>
        /// Parses modifications from a string in the standard mzLib modification format.
        /// Simplified parser for the embedded protease_mods.txt resource.
        /// </summary>
        private static List<Modification> ParseModificationsFromString(string content)
        {
            var modifications = new List<Modification>();

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
                                    existingList.Add(splitDR[1]);
                                else
                                    databaseReference.Add(splitDR[0], new List<string> { splitDR[1] });
                            }
                            break;
                    }
                }

                if (!string.IsNullOrEmpty(id))
                {
                    ModificationMotif motif = null;
                    if (!string.IsNullOrEmpty(target))
                        ModificationMotif.TryGetMotif(target, out motif);

                    modifications.Add(new Modification(
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
                    ));
                }
            }

            return modifications;
        }

        /// <summary>
        /// Reads and parses the embedded proteases.tsv resource.
        /// Throws <see cref="MzLibException"/> if the resource cannot be found.
        /// </summary>
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
                    string[] lines = reader.ReadToEnd()
                        .Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
                    return ParseProteaseLines(lines, proteaseMods);
                }
            }
        }

        #endregion

        #region Custom Digestion Agent Loading

        /// <summary>
        /// Convenience overload for loading a single custom protease file.
        /// See <see cref="LoadAndMergeCustomProteases(IEnumerable{string}, List{Modification})"/> for full details.
        /// </summary>
        /// <param name="path">Path to a custom proteases TSV file.</param>
        /// <param name="proteaseMods">Optional modifications to apply to proteases that require them.</param>
        /// <returns>
        /// A <see cref="CustomDigestionAgentLoadResult"/> describing which proteases were added and which were skipped.
        /// </returns>
        public static CustomDigestionAgentLoadResult LoadAndMergeCustomProteases(
            string path,
            List<Modification> proteaseMods = null)
        {
            return LoadAndMergeCustomProteases(new[] { path }, proteaseMods);
        }

        /// <summary>
        /// Loads custom proteases from one or more TSV files and merges new entries into <see cref="Dictionary"/>.
        ///
        /// <para><b>Merge rules:</b></para>
        /// <list type="bullet">
        ///   <item><description>
        ///     A custom protease whose name already exists in the embedded resource is <b>skipped</b>.
        ///     The embedded definition always wins. The skipped name is recorded in
        ///     <see cref="CustomDigestionAgentLoadResult.Skipped"/>.
        ///   </description></item>
        ///   <item><description>
        ///     A custom protease whose name was already added by an earlier file in the same call is
        ///     likewise skipped (first-encountered wins).
        ///   </description></item>
        ///   <item><description>
        ///     A custom protease with a genuinely new name is added to the dictionary and recorded in
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
        /// <param name="paths">One or more paths to custom proteases TSV files.</param>
        /// <param name="proteaseMods">
        /// Optional modifications to apply to proteases that require cleavage mods.
        /// When null (the default), the embedded protease_mods.txt modifications are used automatically,
        /// so custom files that reference standard modifications (e.g. "Homoserine lactone on M")
        /// work without the caller having to load any mods file.
        /// </param>
        /// <returns>
        /// A <see cref="CustomDigestionAgentLoadResult"/> with the names of added and skipped proteases.
        /// </returns>
        public static CustomDigestionAgentLoadResult LoadAndMergeCustomProteases(
            IEnumerable<string> paths,
            List<Modification> proteaseMods = null)
        {
            // Fall back to the embedded mods when the caller does not supply their own.
            // This allows custom files to reference the same cleavage modifications as the
            // embedded resource (e.g. "Homoserine lactone on M" for CNBr-style proteases)
            // without requiring every caller to load the mods file themselves.
            var resolvedMods = proteaseMods ?? LoadEmbeddedProteaseMods();

            var added = new List<string>();
            var skipped = new List<string>();

            foreach (string path in paths)
            {
                string[] lines = File.ReadAllLines(path);
                var customProteases = ParseProteaseLines(lines, resolvedMods);

                foreach (var kvp in customProteases)
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

        #region Lookup Helpers

        /// <summary>
        /// Gets a protease by name, with backward compatibility for old naming conventions.
        /// Old names like "chymotrypsin (don't cleave before proline)" are automatically
        /// converted to the new format "chymotrypsin|P".
        /// </summary>
        /// <param name="name">The protease name (old or new format).</param>
        /// <returns>The <see cref="Protease"/> object.</returns>
        /// <exception cref="KeyNotFoundException">Thrown when the protease is not found after normalization.</exception>
        public static Protease GetProtease(string name)
        {
            if (Dictionary.TryGetValue(name, out var protease))
                return protease;

            string normalizedName = NormalizeProteaseName(name);
            if (normalizedName != name && Dictionary.TryGetValue(normalizedName, out protease))
                return protease;

            throw new KeyNotFoundException(
                $"Protease '{name}' not found in dictionary. " +
                $"If using an old-style name, ensure it follows the pattern " +
                $"'name (don't cleave before proline)' which maps to 'name|P'.");
        }

        /// <summary>
        /// Tries to get a protease by name, with backward compatibility for old naming conventions.
        /// </summary>
        /// <param name="name">The protease name (old or new format).</param>
        /// <param name="protease">When successful, the <see cref="Protease"/> object; otherwise null.</param>
        /// <returns>True if found; otherwise false.</returns>
        public static bool TryGetProtease(string name, out Protease protease)
        {
            if (Dictionary.TryGetValue(name, out protease))
                return true;

            string normalizedName = NormalizeProteaseName(name);
            if (normalizedName != name && Dictionary.TryGetValue(normalizedName, out protease))
                return true;

            protease = null;
            return false;
        }

        /// <summary>
        /// Normalizes old-style protease names to the current format.
        /// Converts "chymotrypsin (don't cleave before proline)" → "chymotrypsin|P".
        /// </summary>
        public static string NormalizeProteaseName(string name)
        {
            if (string.IsNullOrEmpty(name))
                return name;

            string[] prolineRestrictionPatterns = { " (don't cleave before proline)" };

            foreach (var pattern in prolineRestrictionPatterns)
            {
                int index = name.IndexOf(pattern, StringComparison.OrdinalIgnoreCase);
                if (index >= 0)
                {
                    string baseName = name.Substring(0, index).Trim();
                    return baseName + "|P";
                }
            }

            return name;
        }

        #endregion

        #region TSV Parsing

        /// <summary>
        /// Parses protease definitions from TSV-formatted lines.
        /// Comment lines (starting with '#') and blank lines are skipped.
        /// The first non-comment line is treated as the header and used to resolve column positions.
        /// Duplicate names within the same set of lines throw <see cref="MzLibException"/>.
        /// </summary>
        private static Dictionary<string, Protease> ParseProteaseLines(string[] lines, List<Modification> proteaseMods)
        {
            var dict = new Dictionary<string, Protease>();
            Dictionary<string, int> columnIndices = null;

            // Column name aliases for backward compatibility (first alias is canonical)
            string[] nameAliases = { "name" };
            string[] motifAliases = { "motif", "sequences inducing cleavage" };
            string[] specificityAliases = { "specificity", "cleavage specificity" };
            string[] psiMsAccessionAliases = { "psi-ms accession", "psi-ms accession number" };
            string[] psiMsNameAliases = { "psi-ms name" };
            string[] cleavageModAliases = { "cleavage modification", "cleavage mass shifts" };

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

                int minRequiredColumns = nameAliases.Concat(motifAliases).Concat(specificityAliases)
                    .Select(alias => columnIndices.TryGetValue(alias, out int idx) ? idx + 1 : 0)
                    .Max();

                if (fields.Length < minRequiredColumns)
                {
                    throw new MzLibException(
                        $"Line for protease '{name}' has only {fields.Length} field(s), but the required columns " +
                        $"(Name, Motif, Specificity) extend to column {minRequiredColumns}. " +
                        $"Please ensure the line has enough tab-separated values.");
                }

                if (string.IsNullOrWhiteSpace(specificityField))
                {
                    throw new MzLibException(
                        $"Line for protease '{name}' is missing a value for the required 'Specificity' column.");
                }

                string psiMsAccessionNumber = GetFieldValue(fields, columnIndices, psiMsAccessionAliases);
                string psiMsName = GetFieldValue(fields, columnIndices, psiMsNameAliases);
                string proteaseModDetails = GetFieldValue(fields, columnIndices, cleavageModAliases);

                var motifList = DigestionMotif.ParseDigestionMotifsFromString(motifField);
                var cleavageSpecificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), specificityField, true);

                Protease protease = CreateProtease(
                    name, motifList, cleavageSpecificity,
                    psiMsAccessionNumber, psiMsName,
                    proteaseModDetails, proteaseMods);

                if (dict.ContainsKey(protease.Name))
                {
                    throw new MzLibException(
                        $"More than one protease named '{protease.Name}' exists within the same file. " +
                        $"Please ensure all protease names are unique within a single file.");
                }

                dict.Add(protease.Name, protease);
            }

            if (columnIndices == null)
            {
                throw new MzLibException(
                    "Protease file contains no header line. Expected columns: Name, Motif, Specificity");
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
                    $"Protease file header is missing required column '{displayName}'. " +
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

        private static Protease CreateProtease(
            string name,
            List<DigestionMotif> motifList,
            CleavageSpecificity cleavageSpecificity,
            string psiMsAccessionNumber,
            string psiMsName,
            string proteaseModDetails,
            List<Modification> proteaseMods)
        {
            if (!string.IsNullOrEmpty(proteaseModDetails))
            {
                // A cleavage modification was specified in the file — it must be resolvable.
                // Passing a null mods list while the file names a modification is itself an
                // error: the caller has no way to supply the required modification.
                Modification proteaseModification = proteaseMods?
                    .FirstOrDefault(p => p.IdWithMotif == proteaseModDetails);

                if (proteaseModification != null)
                {
                    return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName,
                        motifList, proteaseModification);
                }

                throw new MzLibException(
                    $"Protease '{name}' specifies cleavage modification '{proteaseModDetails}' " +
                    $"but it was not found in the provided modifications list. " +
                    $"Ensure the modification is defined and the correct mods list is passed.");
            }

            return new Protease(name, cleavageSpecificity, psiMsAccessionNumber, psiMsName, motifList);
        }

        #endregion
    }
}