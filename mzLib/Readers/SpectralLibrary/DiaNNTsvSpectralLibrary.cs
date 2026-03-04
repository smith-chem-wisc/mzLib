using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Reader and writer for DIA-NN's tab-separated (TSV) spectral library format.
    /// 
    /// Format overview:
    /// - Tab-delimited text file with a header row
    /// - Each row represents ONE fragment ion
    /// - A precursor with N fragments occupies N rows
    /// - Precursor-level columns are repeated on every fragment row
    /// - Modified sequences use UniMod notation: _PEPTM[UniMod:35]IDE_
    /// 
    /// This class supports:
    /// - Reading DIA-NN TSV libraries into List&lt;DiaNNLibraryEntry&gt;
    /// - Writing DiaNNLibraryEntry collections to DIA-NN-compatible TSV
    /// - Converting to/from mzLib List&lt;LibrarySpectrum&gt;
    /// - Recognizing common column name aliases (Spectronaut, OpenSWATH)
    /// 
    /// Column schema documented in architecture doc Section 5.
    /// </summary>
    public class DiaNNTsvSpectralLibrary
    {
        #region Column Name Constants and Aliases

        // Primary DIA-NN column names (used for writing)
        private const string Col_ModifiedPeptide = "ModifiedPeptide";
        private const string Col_StrippedPeptide = "StrippedPeptide";
        private const string Col_PrecursorCharge = "PrecursorCharge";
        private const string Col_PrecursorMz = "PrecursorMz";
        private const string Col_Tr_recalibrated = "Tr_recalibrated";
        private const string Col_IonMobility = "IonMobility";
        private const string Col_ProteinId = "ProteinId";
        private const string Col_ProteinName = "ProteinName";
        private const string Col_Genes = "Genes";
        private const string Col_Proteotypic = "Proteotypic";
        private const string Col_Decoy = "Decoy";
        private const string Col_FragmentMz = "FragmentMz";
        private const string Col_RelativeIntensity = "RelativeIntensity";
        private const string Col_FragmentType = "FragmentType";
        private const string Col_FragmentNumber = "FragmentNumber";
        private const string Col_FragmentCharge = "FragmentCharge";
        private const string Col_FragmentLossType = "FragmentLossType";
        private const string Col_ExcludeFromAssay = "ExcludeFromAssay";

        /// <summary>
        /// Maps common column name aliases to canonical DIA-NN names.
        /// DIA-NN accepts Spectronaut-style and OpenSWATH-style column headers.
        /// Case-insensitive matching is applied during reading.
        /// </summary>
        private static readonly Dictionary<string, string> ColumnAliases = new(StringComparer.OrdinalIgnoreCase)
        {
            // ModifiedPeptide aliases
            { "ModifiedPeptideSequence", Col_ModifiedPeptide },
            { "FullUniModPeptideName", Col_ModifiedPeptide },

            // PrecursorCharge aliases
            { "Charge", Col_PrecursorCharge },
            { "PrecursorZ", Col_PrecursorCharge },

            // Retention time aliases
            { "iRT", Col_Tr_recalibrated },
            { "RetentionTime", Col_Tr_recalibrated },
            { "NormalizedRetentionTime", Col_Tr_recalibrated },
            { "RT", Col_Tr_recalibrated },

            // Intensity aliases
            { "LibraryIntensity", Col_RelativeIntensity },
            { "Intensity", Col_RelativeIntensity },

            // Protein aliases
            { "UniprotID", Col_ProteinId },
            { "ProteinAccession", Col_ProteinId },
            { "Protein", Col_ProteinName },
            { "ProteinDescription", Col_ProteinName },

            // Gene aliases
            { "GeneName", Col_Genes },
            { "Gene", Col_Genes },
        };

        #endregion

        #region Read: TSV → DiaNNLibraryEntry

        /// <summary>
        /// Reads a DIA-NN TSV spectral library file into a list of DiaNNLibraryEntry objects.
        /// 
        /// Groups fragment rows by their precursor key (ModifiedSequence + PrecursorCharge)
        /// to reconstruct full spectra.
        /// 
        /// Handles:
        /// - Column name aliases (Spectronaut, OpenSWATH formats)
        /// - Missing optional columns (IonMobility, ProteinId, Genes, etc.)
        /// - Both tab and comma delimiters (auto-detected from header)
        /// </summary>
        /// <param name="filePath">Path to the DIA-NN TSV library file</param>
        /// <returns>List of DiaNNLibraryEntry, one per unique precursor</returns>
        /// <exception cref="FileNotFoundException">If the file does not exist</exception>
        /// <exception cref="FormatException">If required columns are missing from the header</exception>
        public static List<DiaNNLibraryEntry> ReadTsv(string filePath)
        {
            if (!File.Exists(filePath))
                throw new FileNotFoundException($"DIA-NN TSV library file not found: {filePath}", filePath);

            var entries = new Dictionary<string, DiaNNLibraryEntry>(); // keyed by ModifiedSequence/Charge

            using var reader = new StreamReader(filePath, Encoding.UTF8);

            // Read and parse header line
            string headerLine = reader.ReadLine();
            if (string.IsNullOrEmpty(headerLine))
                throw new FormatException("DIA-NN TSV file is empty or has no header.");

            // Auto-detect delimiter
            char delimiter = headerLine.Contains('\t') ? '\t' : ',';
            string[] rawHeaders = headerLine.Split(delimiter);

            // Resolve column indices, applying aliases for non-canonical names
            var columnIndex = ResolveColumnIndices(rawHeaders);

            // Validate required columns
            ValidateRequiredColumns(columnIndex);

            // Read data rows
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                if (string.IsNullOrWhiteSpace(line))
                    continue;

                string[] fields = line.Split(delimiter);

                // Parse precursor-level data
                string modifiedSequence = GetField(fields, columnIndex, Col_ModifiedPeptide);
                int precursorCharge = GetIntField(fields, columnIndex, Col_PrecursorCharge);
                string precursorKey = modifiedSequence + "/" + precursorCharge;

                // Parse fragment-level data
                var fragment = new DiaNNFragmentIon
                {
                    Mz = GetDoubleField(fields, columnIndex, Col_FragmentMz),
                    Intensity = GetDoubleField(fields, columnIndex, Col_RelativeIntensity),
                    IonType = GetField(fields, columnIndex, Col_FragmentType).ToLowerInvariant()[0],
                    SeriesNumber = GetIntField(fields, columnIndex, Col_FragmentNumber),
                    Charge = GetIntField(fields, columnIndex, Col_FragmentCharge),
                    LossType = GetFieldOrDefault(fields, columnIndex, Col_FragmentLossType, "noloss"),
                    ExcludeFromAssay = GetIntFieldOrDefault(fields, columnIndex, Col_ExcludeFromAssay, 0) == 1,
                };

                // If this precursor already exists, just add the fragment
                if (entries.TryGetValue(precursorKey, out var existingEntry))
                {
                    existingEntry.Fragments.Add(fragment);
                    continue;
                }

                // New precursor — parse all precursor-level fields
                var entry = new DiaNNLibraryEntry
                {
                    ModifiedSequence = modifiedSequence,
                    StrippedSequence = GetFieldOrDefault(fields, columnIndex, Col_StrippedPeptide,
                        DiaNNModificationMapping.GetStrippedSequence(modifiedSequence)),
                    PrecursorMz = GetDoubleField(fields, columnIndex, Col_PrecursorMz),
                    PrecursorCharge = precursorCharge,
                    RetentionTime = GetDoubleField(fields, columnIndex, Col_Tr_recalibrated),
                    IonMobility = GetDoubleFieldOrDefault(fields, columnIndex, Col_IonMobility, 0.0),
                    ProteinId = GetFieldOrDefault(fields, columnIndex, Col_ProteinId, string.Empty),
                    ProteinName = GetFieldOrDefault(fields, columnIndex, Col_ProteinName, string.Empty),
                    GeneName = GetFieldOrDefault(fields, columnIndex, Col_Genes, string.Empty),
                    IsProteotypic = GetIntFieldOrDefault(fields, columnIndex, Col_Proteotypic, 1) == 1,
                    IsDecoy = GetIntFieldOrDefault(fields, columnIndex, Col_Decoy, 0) == 1,
                    QValue = null, // Not typically in TSV; populated from search results
                    Fragments = new List<DiaNNFragmentIon> { fragment },
                };

                entries[precursorKey] = entry;
            }

            return entries.Values.ToList();
        }

        /// <summary>
        /// Reads a DIA-NN TSV file and converts directly to mzLib LibrarySpectrum objects.
        /// Convenience method for workflows that don't need DIA-NN-specific metadata.
        /// </summary>
        /// <param name="filePath">Path to the DIA-NN TSV library file</param>
        /// <returns>List of mzLib LibrarySpectrum objects</returns>
        public static List<LibrarySpectrum> ReadTsvAsLibrarySpectra(string filePath)
        {
            return ReadTsv(filePath)
                .Select(entry => entry.ToLibrarySpectrum())
                .ToList();
        }

        #endregion

        #region Write: DiaNNLibraryEntry → TSV

        /// <summary>
        /// Writes a collection of DiaNNLibraryEntry objects to a DIA-NN-compatible TSV file.
        /// 
        /// Output format:
        /// - Tab-delimited with standard DIA-NN column headers
        /// - One row per fragment ion
        /// - Precursor columns repeated on each fragment row
        /// - Numeric values formatted for DIA-NN compatibility (InvariantCulture)
        /// </summary>
        /// <param name="entries">Library entries to write</param>
        /// <param name="filePath">Output file path</param>
        /// <param name="includeIonMobility">Whether to include the IonMobility column (default true)</param>
        public static void WriteTsv(IEnumerable<DiaNNLibraryEntry> entries, string filePath,
            bool includeIonMobility = true)
        {
            using var writer = new StreamWriter(filePath, false, new UTF8Encoding(false));

            // Write header
            var columns = GetWriteColumns(includeIonMobility);
            writer.WriteLine(string.Join('\t', columns));

            // Write data rows
            foreach (var entry in entries)
            {
                foreach (var fragment in entry.Fragments)
                {
                    var row = FormatRow(entry, fragment, includeIonMobility);
                    writer.WriteLine(string.Join('\t', row));
                }
            }
        }

        /// <summary>
        /// Writes mzLib LibrarySpectrum objects to a DIA-NN-compatible TSV file.
        /// Converts each spectrum to DiaNNLibraryEntry first, then writes.
        /// DIA-NN-specific metadata will use default values.
        /// </summary>
        /// <param name="spectra">mzLib library spectra to write</param>
        /// <param name="filePath">Output file path</param>
        public static void WriteTsvFromLibrarySpectra(IEnumerable<LibrarySpectrum> spectra, string filePath)
        {
            var entries = spectra.Select(DiaNNLibraryEntry.FromLibrarySpectrum);
            WriteTsv(entries, filePath, includeIonMobility: false);
        }

        #endregion

        #region Column Resolution

        /// <summary>
        /// Resolves column header names to their indices, applying aliases for non-canonical names.
        /// Returns a dictionary mapping canonical column names to their 0-based index.
        /// </summary>
        private static Dictionary<string, int> ResolveColumnIndices(string[] headers)
        {
            var result = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            for (int i = 0; i < headers.Length; i++)
            {
                string header = headers[i].Trim();

                // Try direct match first (canonical DIA-NN name)
                if (!result.ContainsKey(header))
                {
                    result[header] = i;
                }

                // Try alias resolution
                if (ColumnAliases.TryGetValue(header, out var canonicalName))
                {
                    if (!result.ContainsKey(canonicalName))
                    {
                        result[canonicalName] = i;
                    }
                }
            }

            return result;
        }

        /// <summary>
        /// Validates that all required columns are present in the resolved column index.
        /// Throws FormatException with a clear message listing missing columns.
        /// </summary>
        private static void ValidateRequiredColumns(Dictionary<string, int> columnIndex)
        {
            string[] requiredColumns = new[]
            {
                Col_ModifiedPeptide,
                Col_PrecursorCharge,
                Col_PrecursorMz,
                Col_Tr_recalibrated,
                Col_FragmentMz,
                Col_RelativeIntensity,
                Col_FragmentType,
                Col_FragmentNumber,
                Col_FragmentCharge,
            };

            var missing = requiredColumns.Where(c => !columnIndex.ContainsKey(c)).ToList();
            if (missing.Count > 0)
            {
                throw new FormatException(
                    $"DIA-NN TSV library is missing required columns: {string.Join(", ", missing)}. " +
                    $"Found columns: {string.Join(", ", columnIndex.Keys)}");
            }
        }

        #endregion

        #region Field Accessors

        private static string GetField(string[] fields, Dictionary<string, int> columnIndex, string columnName)
        {
            if (!columnIndex.TryGetValue(columnName, out int idx) || idx >= fields.Length)
                throw new FormatException($"Required column '{columnName}' not found or row too short.");
            return fields[idx].Trim();
        }

        private static string GetFieldOrDefault(string[] fields, Dictionary<string, int> columnIndex,
            string columnName, string defaultValue)
        {
            if (!columnIndex.TryGetValue(columnName, out int idx) || idx >= fields.Length)
                return defaultValue;
            string value = fields[idx].Trim();
            return string.IsNullOrEmpty(value) ? defaultValue : value;
        }

        private static int GetIntField(string[] fields, Dictionary<string, int> columnIndex, string columnName)
        {
            string value = GetField(fields, columnIndex, columnName);
            if (!int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out int result))
                throw new FormatException($"Cannot parse '{value}' as integer in column '{columnName}'.");
            return result;
        }

        private static int GetIntFieldOrDefault(string[] fields, Dictionary<string, int> columnIndex,
            string columnName, int defaultValue)
        {
            if (!columnIndex.TryGetValue(columnName, out int idx) || idx >= fields.Length)
                return defaultValue;
            string value = fields[idx].Trim();
            if (string.IsNullOrEmpty(value))
                return defaultValue;
            return int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out int result)
                ? result
                : defaultValue;
        }

        private static double GetDoubleField(string[] fields, Dictionary<string, int> columnIndex, string columnName)
        {
            string value = GetField(fields, columnIndex, columnName);
            if (!double.TryParse(value, NumberStyles.Float | NumberStyles.AllowExponent,
                CultureInfo.InvariantCulture, out double result))
                throw new FormatException($"Cannot parse '{value}' as double in column '{columnName}'.");
            return result;
        }

        private static double GetDoubleFieldOrDefault(string[] fields, Dictionary<string, int> columnIndex,
            string columnName, double defaultValue)
        {
            if (!columnIndex.TryGetValue(columnName, out int idx) || idx >= fields.Length)
                return defaultValue;
            string value = fields[idx].Trim();
            if (string.IsNullOrEmpty(value))
                return defaultValue;
            return double.TryParse(value, NumberStyles.Float | NumberStyles.AllowExponent,
                CultureInfo.InvariantCulture, out double result)
                ? result
                : defaultValue;
        }

        #endregion

        #region Write Helpers

        /// <summary>
        /// Returns the ordered list of column names for the TSV header.
        /// </summary>
        private static string[] GetWriteColumns(bool includeIonMobility)
        {
            var columns = new List<string>
            {
                Col_ModifiedPeptide,
                Col_StrippedPeptide,
                Col_PrecursorCharge,
                Col_PrecursorMz,
                Col_Tr_recalibrated,
            };

            if (includeIonMobility)
                columns.Add(Col_IonMobility);

            columns.AddRange(new[]
            {
                Col_ProteinId,
                Col_ProteinName,
                Col_Genes,
                Col_Proteotypic,
                Col_Decoy,
                Col_FragmentMz,
                Col_RelativeIntensity,
                Col_FragmentType,
                Col_FragmentNumber,
                Col_FragmentCharge,
                Col_FragmentLossType,
                Col_ExcludeFromAssay,
            });

            return columns.ToArray();
        }

        /// <summary>
        /// Formats a single TSV row for one fragment of one precursor.
        /// All numeric values use InvariantCulture formatting.
        /// </summary>
        private static string[] FormatRow(DiaNNLibraryEntry entry, DiaNNFragmentIon fragment,
            bool includeIonMobility)
        {
            var fields = new List<string>
            {
                entry.ModifiedSequence,
                entry.StrippedSequence ?? DiaNNModificationMapping.GetStrippedSequence(entry.ModifiedSequence),
                entry.PrecursorCharge.ToString(CultureInfo.InvariantCulture),
                entry.PrecursorMz.ToString("F6", CultureInfo.InvariantCulture),
                entry.RetentionTime.ToString("F4", CultureInfo.InvariantCulture),
            };

            if (includeIonMobility)
                fields.Add(entry.IonMobility.ToString("F4", CultureInfo.InvariantCulture));

            fields.AddRange(new[]
            {
                entry.ProteinId ?? string.Empty,
                entry.ProteinName ?? string.Empty,
                entry.GeneName ?? string.Empty,
                entry.IsProteotypic ? "1" : "0",
                entry.IsDecoy ? "1" : "0",
                fragment.Mz.ToString("F6", CultureInfo.InvariantCulture),
                fragment.Intensity.ToString("F6", CultureInfo.InvariantCulture),
                fragment.IonType.ToString(),
                fragment.SeriesNumber.ToString(CultureInfo.InvariantCulture),
                fragment.Charge.ToString(CultureInfo.InvariantCulture),
                fragment.LossType ?? "noloss",
                fragment.ExcludeFromAssay ? "1" : "0",
            });

            return fields.ToArray();
        }

        #endregion
    }
}
