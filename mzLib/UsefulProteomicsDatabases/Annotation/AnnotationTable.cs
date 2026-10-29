using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.Annotation
{
    /// <summary>
    /// A table read back from one of mzLib's annotation files: its "#!key value" header and its rows.
    /// </summary>
    /// <param name="SourceFileName">The file name the table was read from.</param>
    /// <param name="Header">Every header key and its value, the format line included.</param>
    /// <param name="Rows">The rows, in file order.</param>
    public sealed record AnnotationTable<T>(string SourceFileName, IReadOnlyDictionary<string, string> Header, IReadOnlyList<T> Rows);

    /// <summary>
    /// Reads the annotation tables' shared layout: "#!key value" header lines, the format key first, then
    /// one column header, then rows. Strict throughout -- an unknown or repeated key, a missing required key,
    /// another format version, a column header that is not exactly the schema, or a row of the wrong width
    /// is refused with its line number -- because a table read wrongly is joined wrongly, and nothing
    /// downstream could tell.
    /// </summary>
    internal static class AnnotationTableReader
    {
        /// <summary>Header keys every writer may add: the mzLib build that wrote the file and its release.</summary>
        private static readonly string[] ProducerKeys = { "mzlib_version", "mzlib_release" };

        public static AnnotationTable<T> Read<T>(string path, string formatKey, int formatVersion,
            IReadOnlyCollection<string> requiredKeys, IReadOnlyCollection<string> optionalKeys,
            IReadOnlyList<string> columns, Func<string[], string, int, T> parseRow)
        {
            ArgumentNullException.ThrowIfNull(path);
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Annotation table not found.", path);
            }
            string fileName = Path.GetFileName(path);
            var allowed = new HashSet<string>(requiredKeys.Concat(optionalKeys).Concat(ProducerKeys).Append(formatKey), StringComparer.Ordinal);
            var header = new Dictionary<string, string>(StringComparer.Ordinal);
            var rows = new List<T>();
            bool columnsRead = false;

            using var reader = new StreamReader(path);
            int lineNumber = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                lineNumber++;
                if (!columnsRead)
                {
                    if (line.StartsWith("#!", StringComparison.Ordinal))
                    {
                        string[] parts = line.Substring(2).Split(' ', 2);
                        string key = parts[0];
                        if (header.Count == 0 && key != formatKey)
                        {
                            throw new InvalidDataException($"{fileName} line {lineNumber}: the first header line must be '#!{formatKey}', found '#!{key}'.");
                        }
                        if (!allowed.Contains(key))
                        {
                            throw new InvalidDataException($"{fileName} line {lineNumber}: unknown header key '{key}'.");
                        }
                        if (!header.TryAdd(key, parts.Length == 2 ? parts[1] : ""))
                        {
                            throw new InvalidDataException($"{fileName} line {lineNumber}: header key '{key}' appears twice.");
                        }
                        continue;
                    }
                    CheckHeader(header, fileName, formatKey, formatVersion, requiredKeys);
                    if (line != string.Join('\t', columns))
                    {
                        throw new InvalidDataException(
                            $"{fileName} line {lineNumber}: expected the column header '{string.Join(", ", columns)}', found '{line.Replace('\t', ',')}'.");
                    }
                    columnsRead = true;
                    continue;
                }
                if (line.Length == 0)
                {
                    continue;
                }
                string[] cells = line.Split('\t');
                if (cells.Length != columns.Count)
                {
                    throw new InvalidDataException($"{fileName} line {lineNumber}: expected {columns.Count} cells, found {cells.Length}.");
                }
                rows.Add(parseRow(cells, fileName, lineNumber));
            }

            if (!columnsRead)
            {
                CheckHeader(header, fileName, formatKey, formatVersion, requiredKeys);
                throw new InvalidDataException($"{fileName}: no column header.");
            }
            return new AnnotationTable<T>(fileName, header, rows);
        }

        private static void CheckHeader(Dictionary<string, string> header, string fileName, string formatKey, int formatVersion,
            IReadOnlyCollection<string> requiredKeys)
        {
            if (!header.TryGetValue(formatKey, out var format) || format != formatVersion.ToString(System.Globalization.CultureInfo.InvariantCulture))
            {
                throw new InvalidDataException($"{fileName}: {formatKey} must be {formatVersion}, found '{format ?? "(missing)"}'.");
            }
            foreach (string required in requiredKeys)
            {
                if (!header.TryGetValue(required, out var value) || value.Length == 0)
                {
                    throw new InvalidDataException($"{fileName}: header key '{required}' is missing.");
                }
            }
        }

        /// <summary>Null for an empty cell.</summary>
        public static string Cell(string value) => value.Length == 0 ? null : value;

        /// <summary>A ';'-joined set; empty for an empty cell.</summary>
        public static IReadOnlyList<string> Set(string value) =>
            value.Length == 0 ? Array.Empty<string>() : value.Split(';');

        /// <summary>"true", "false", or null for an empty cell.</summary>
        public static bool? Bool(string value, string column, string fileName, int lineNumber) => value switch
        {
            "true" => true,
            "false" => false,
            "" => null,
            _ => throw new InvalidDataException($"{fileName} line {lineNumber}: {column} '{value}' is not true or false.")
        };
    }
}
