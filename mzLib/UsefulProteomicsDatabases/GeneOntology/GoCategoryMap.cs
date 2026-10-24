using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// One row of a category map: a GO term that anchors a category, and optionally a subcategory within it.
    /// </summary>
    /// <param name="Category">The category label, e.g. "mitochondrion". Never empty.</param>
    /// <param name="Subcategory">The subcategory label as written in the map, unqualified (e.g.
    /// "inner_membrane"), or null when the row anchors the category itself.</param>
    /// <param name="AnchorGoId">The anchoring term as written in the map; it may be an alternative id.</param>
    public sealed record GoCategoryAnchor(string Category, string Subcategory, string AnchorGoId);

    /// <summary>
    /// A consumer's term-to-category map: which GO terms anchor which categories.
    ///
    /// mzLib ships no map. The rows are the consumer's science -- an organelle vocabulary, a complex list,
    /// anything -- and mzLib owns only the file format and the rule that applies it
    /// (<see cref="GoCategoryResolver"/>). The file is:
    /// <code>
    /// #!category_map_format 1
    /// #!map_name organelle
    /// #!map_version 1
    /// category	subcategory	anchor_go_id
    /// mitochondrion		GO:0005739
    /// mitochondrion	inner_membrane	GO:0005743
    /// </code>
    /// An empty subcategory anchors the category itself. Several rows may share a label, so one category
    /// can be reached by unrelated terms. Labels never contain ':' -- the resolver writes a subcategory as
    /// "category:subcategory", and a label already carrying one would be qualified twice.
    ///
    /// The map records its file name and sha256 as well as its declared name and version, because a
    /// version string is a claim and the hash is what shows two files are the same map.
    /// </summary>
    public sealed class GoCategoryMap
    {
        /// <summary>The only format version this reader understands.</summary>
        public const string FormatVersion = "1";

        private static readonly string[] Columns = { "category", "subcategory", "anchor_go_id" };
        private static readonly Regex GoId = new(@"^GO:\d{7}$", RegexOptions.Compiled);

        private GoCategoryMap(string mapName, string mapVersion, IReadOnlyList<GoCategoryAnchor> anchors,
            string sourceFileName, string sourceSha256)
        {
            MapName = mapName;
            MapVersion = mapVersion;
            Anchors = anchors;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
        }

        /// <summary>The map's declared name, e.g. "organelle".</summary>
        public string MapName { get; }

        /// <summary>The map's declared version. The consumer's to bump; mzLib only records it.</summary>
        public string MapVersion { get; }

        /// <summary>Every row, in file order.</summary>
        public IReadOnlyList<GoCategoryAnchor> Anchors { get; }

        /// <summary>The file name the map was read from.</summary>
        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read.</summary>
        public string SourceSha256 { get; }

        /// <summary>
        /// Reads a category map. Strict throughout: a map that is read wrongly assigns categories wrongly,
        /// and nothing downstream could tell.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">
        /// The format version is missing or not <see cref="FormatVersion"/>; map_name or map_version is
        /// missing; a header key is unknown or repeated; the column header is not exactly
        /// "category, subcategory, anchor_go_id"; or a row has the wrong number of cells, an empty category,
        /// an anchor that is not a GO id, a ':' in a label, or repeats an earlier row.
        /// </exception>
        public static GoCategoryMap Load(string path)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Category map not found.", path);
            }

            string fileName = Path.GetFileName(path);
            string sha256;
            using (var hashStream = File.OpenRead(path))
            {
                sha256 = Convert.ToHexString(SHA256.HashData(hashStream)).ToLowerInvariant();
            }

            var header = new Dictionary<string, string>(StringComparer.Ordinal);
            var anchors = new List<GoCategoryAnchor>();
            var seen = new HashSet<GoCategoryAnchor>();
            bool columnsRead = false;

            using var reader = new StreamReader(path);
            int lineNumber = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                lineNumber++;
                if (line.Length == 0)
                {
                    continue;
                }

                if (!columnsRead)
                {
                    if (line.StartsWith("#!", StringComparison.Ordinal))
                    {
                        ReadHeaderLine(line, header, fileName, lineNumber);
                        continue;
                    }
                    CheckHeader(header, fileName);
                    if (line != string.Join('\t', Columns))
                    {
                        throw new InvalidDataException(
                            $"{fileName} line {lineNumber}: expected the column header '{string.Join(", ", Columns)}', found '{line.Replace('\t', ',')}'.");
                    }
                    columnsRead = true;
                    continue;
                }

                var anchor = ReadRow(line, fileName, lineNumber);
                if (!seen.Add(anchor))
                {
                    throw new InvalidDataException($"{fileName} line {lineNumber}: repeats an earlier row.");
                }
                anchors.Add(anchor);
            }

            if (!columnsRead)
            {
                CheckHeader(header, fileName);
                throw new InvalidDataException($"{fileName}: no column header.");
            }

            return new GoCategoryMap(header["map_name"], header["map_version"], anchors, fileName, sha256);
        }

        private static void ReadHeaderLine(string line, Dictionary<string, string> header, string fileName, int lineNumber)
        {
            string[] parts = line.Substring(2).Split(' ', 2);
            string key = parts[0];
            string value = parts.Length == 2 ? parts[1].Trim() : "";
            if (key != "category_map_format" && key != "map_name" && key != "map_version")
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: unknown header key '{key}'.");
            }
            if (value.Length == 0)
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: header key '{key}' has no value.");
            }
            // The category table's header writes "#!category_map <name> <version> <sha256>", split on spaces.
            if (value.Any(char.IsWhiteSpace))
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: header key '{key}' may not contain whitespace.");
            }
            if (!header.TryAdd(key, value))
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: header key '{key}' appears twice.");
            }
        }

        private static void CheckHeader(Dictionary<string, string> header, string fileName)
        {
            if (!header.TryGetValue("category_map_format", out var format) || format != FormatVersion)
            {
                throw new InvalidDataException(
                    $"{fileName}: category_map_format must be {FormatVersion}, found '{format ?? "(missing)"}'.");
            }
            foreach (string required in new[] { "map_name", "map_version" })
            {
                if (!header.ContainsKey(required))
                {
                    throw new InvalidDataException($"{fileName}: header key '{required}' is missing.");
                }
            }
        }

        private static GoCategoryAnchor ReadRow(string line, string fileName, int lineNumber)
        {
            string[] cells = line.Split('\t');
            if (cells.Length != Columns.Length)
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: expected {Columns.Length} cells, found {cells.Length}.");
            }
            string category = cells[0].Trim();
            string subcategory = cells[1].Trim();
            string anchor = cells[2].Trim();

            if (category.Length == 0)
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: category is empty.");
            }
            if (category.Contains(':') || subcategory.Contains(':'))
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: labels may not contain ':'.");
            }
            if (!GoId.IsMatch(anchor))
            {
                throw new InvalidDataException($"{fileName} line {lineNumber}: anchor '{anchor}' is not a GO id.");
            }
            return new GoCategoryAnchor(category, subcategory.Length == 0 ? null : subcategory, anchor);
        }
    }
}
