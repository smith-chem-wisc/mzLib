using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>What kind of homology Compara asserts, derived from its homology_type.</summary>
    public enum HomologyClass
    {
        /// <summary>ortholog_one2one, ortholog_one2many or ortholog_many2many: between two species.</summary>
        Ortholog,

        /// <summary>within_species_paralog, other_paralog or gene_split: within one species.</summary>
        Paralog,

        /// <summary>homoeolog_*: between subgenomes of one polyploid species.</summary>
        Homoeolog
    }

    /// <summary>
    /// One relationship as Compara states it, values verbatim. A NULL in the dump is null here, never 0:
    /// a source declining to say is not a source saying "no". Side A is whichever gene the row lists
    /// first; <see cref="OrthologySnapshot"/> orients cross-species rows so that A is the species that
    /// sorts first.
    /// </summary>
    /// <param name="HomologyId">Compara's id for the relationship; unique within a release.</param>
    /// <param name="HomologyType">Compara's homology_type, e.g. "ortholog_one2many".</param>
    /// <param name="IdentityA">Percentage of A's protein identical to B's in the alignment.</param>
    /// <param name="IdentityB">Percentage of B's protein identical to A's.</param>
    /// <param name="GocScore">Gene order conservation score, 0 to 100.</param>
    /// <param name="WgaCoverage">Whole-genome alignment coverage, 0 to 100.</param>
    /// <param name="IsHighConfidence">Compara's high-confidence flag; null when it gives none.</param>
    /// <param name="SourceFile">The dump the row was read from. Which of a pair's two dumps holds a
    /// relationship is arbitrary, so this is provenance, not meaning.</param>
    public sealed record ComparaHomology(
        string HomologyId,
        string HomologyType,
        HomologyClass Class,
        string SpeciesA, string GeneA, string ProteinA, double? IdentityA,
        string SpeciesB, string GeneB, string ProteinB, double? IdentityB,
        double? Dn, double? Ds, int? GocScore, double? WgaCoverage, bool? IsHighConfidence,
        string SourceFile)
    {
        /// <summary>The same relationship with its sides exchanged; each identity moves with its gene.</summary>
        public ComparaHomology Swapped() => this with
        {
            SpeciesA = SpeciesB, GeneA = GeneB, ProteinA = ProteinB, IdentityA = IdentityB,
            SpeciesB = SpeciesA, GeneB = GeneA, ProteinB = ProteinA, IdentityB = IdentityA
        };

        /// <summary>True when the two rows state the same relationship, whichever dump each came from.</summary>
        public bool SameRelationship(ComparaHomology other) =>
            other != null && this with { SourceFile = null } == other with { SourceFile = null };
    }

    /// <summary>
    /// One of Ensembl Compara's genome-specific homology dumps
    /// (homologies/&lt;species&gt;/Compara.&lt;release&gt;.protein_default.homologies.tsv.gz), holding the rows
    /// whose two species are both in a requested set.
    ///
    /// Compara's README warns that each genome's dump holds an ARBITRARY subset of that genome's
    /// relationships, and in practice they partition: in release 116 every human-mouse orthology is in
    /// the mouse dump and none is in the human one. Reading the dump named after a species therefore
    /// gives a complete-looking, badly incomplete answer. <see cref="OrthologySnapshot.Build"/> requires
    /// the dump of every species it is given for exactly this reason.
    ///
    /// The reader refuses rather than guesses. A different column layout, a short row, a value that is
    /// not a number, a homology_type it does not know, a paralog between two species, an ortholog within
    /// one, or a row whose first species is not the dump's own genome each throw.
    /// </summary>
    public sealed class ComparaHomologyDump
    {
        /// <summary>The column layout this reader accepts. Anything else is refused rather than guessed at.</summary>
        public static readonly IReadOnlyList<string> Columns = new[]
        {
            "gene_stable_id", "protein_stable_id", "species", "identity", "homology_type",
            "homology_gene_stable_id", "homology_protein_stable_id", "homology_species", "homology_identity",
            "dn", "ds", "goc_score", "wga_coverage", "is_high_confidence", "homology_id"
        };

        private static readonly Dictionary<string, HomologyClass> Classes = new(StringComparer.Ordinal)
        {
            ["ortholog_one2one"] = HomologyClass.Ortholog,
            ["ortholog_one2many"] = HomologyClass.Ortholog,
            ["ortholog_many2many"] = HomologyClass.Ortholog,
            ["within_species_paralog"] = HomologyClass.Paralog,
            ["other_paralog"] = HomologyClass.Paralog,
            ["gene_split"] = HomologyClass.Paralog,
            ["homoeolog_one2one"] = HomologyClass.Homoeolog,
            ["homoeolog_one2many"] = HomologyClass.Homoeolog,
            ["homoeolog_many2many"] = HomologyClass.Homoeolog,
        };

        /// <summary>Ensembl names the dumps Compara.&lt;release&gt;.&lt;collection&gt;.homologies.tsv[.gz].</summary>
        private static readonly Regex ReleaseInFileName =
            new(@"^Compara\.(\d+)\.", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private ComparaHomologyDump(string genome, IReadOnlyList<ComparaHomology> rows, IReadOnlyList<string> species,
            string sourceFileName, string sourceSha256, string release, int rowsRead)
        {
            Genome = genome;
            Rows = rows;
            Species = species;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
            RowsRead = rowsRead;
        }

        /// <summary>The genome the dump belongs to: the species column of every row. Null for an empty dump.</summary>
        public string Genome { get; }

        /// <summary>The rows kept: both species in <see cref="Species"/>. In file order, sides as written.</summary>
        public IReadOnlyList<ComparaHomology> Rows { get; }

        /// <summary>The species the rows were restricted to, in ordinal order.</summary>
        public IReadOnlyList<string> Species { get; }

        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a .gz).</summary>
        public string SourceSha256 { get; }

        /// <summary>The Ensembl release parsed from the file name, or null when it carries none.</summary>
        public string Release { get; }

        /// <summary>Every data row read and checked, kept or not.</summary>
        public int RowsRead { get; }

        /// <summary>The homology_type values this reader knows. Any other value is refused.</summary>
        public static IReadOnlyCollection<string> KnownHomologyTypes => Classes.Keys;

        /// <summary>
        /// Reads the dump, keeping the rows whose two species are both in <paramref name="species"/>. Every
        /// row is checked, kept or not, so a malformed file is refused whichever species are asked for.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file breaks one of the rules in the class summary.</exception>
        public static ComparaHomologyDump Load(string path, IEnumerable<string> species)
        {
            ArgumentNullException.ThrowIfNull(species);
            var keep = new HashSet<string>(species, StringComparer.Ordinal);
            if (keep.Count == 0)
            {
                throw new ArgumentException("At least one species is required.", nameof(species));
            }

            string sha256 = EnsemblFile.Sha256(path, "Compara homology dump");
            string name = Path.GetFileName(path);
            using var reader = EnsemblFile.OpenText(path);

            string header = reader.ReadLine() ?? "";
            if (!header.Split('\t').SequenceEqual(Columns))
            {
                throw new InvalidDataException(
                    $"{name}: unexpected columns. Got [{header.Replace('\t', ',')}], expected [{string.Join(",", Columns)}].");
            }

            var rows = new List<ComparaHomology>();
            string genome = null;
            int lineNumber = 1;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                lineNumber++;
                string[] c = line.Split('\t');
                if (c.Length != Columns.Count)
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: {c.Length} cells, expected {Columns.Count}.");
                }

                string speciesA = c[2], speciesB = c[7], type = c[4], id = c[14];
                genome ??= speciesA;
                if (speciesA != genome)
                {
                    throw new InvalidDataException(
                        $"{name} line {lineNumber}: species {speciesA} in a dump of {genome}.");
                }

                if (!Classes.TryGetValue(type, out var cls))
                {
                    throw new InvalidDataException(
                        $"{name} line {lineNumber}: unknown homology_type '{type}'. Read the release notes before adding it.");
                }

                if (cls == HomologyClass.Paralog && speciesA != speciesB)
                {
                    throw new InvalidDataException(
                        $"{name} line {lineNumber}: homology {id} is a {type} between {speciesA} and {speciesB}.");
                }

                if (cls == HomologyClass.Ortholog && speciesA == speciesB)
                {
                    throw new InvalidDataException(
                        $"{name} line {lineNumber}: homology {id} is a {type} within {speciesA}.");
                }

                if (!keep.Contains(speciesA) || !keep.Contains(speciesB))
                {
                    continue;
                }

                rows.Add(new ComparaHomology(
                    id, type, cls,
                    speciesA, c[0], c[1], Real(c[3], "identity", name, lineNumber),
                    speciesB, c[5], c[6], Real(c[8], "homology_identity", name, lineNumber),
                    Real(c[9], "dn", name, lineNumber), Real(c[10], "ds", name, lineNumber),
                    Whole(c[11], "goc_score", name, lineNumber), Real(c[12], "wga_coverage", name, lineNumber),
                    Flag(c[13], name, lineNumber), name));
            }

            var release = ReleaseInFileName.Match(name);
            return new ComparaHomologyDump(genome, rows, keep.OrderBy(s => s, StringComparer.Ordinal).ToList(),
                name, sha256, release.Success ? release.Groups[1].Value : null, lineNumber - 1);
        }

        /// <summary>The class of a homology_type, or an exception naming the value.</summary>
        public static HomologyClass ClassOf(string homologyType) =>
            homologyType != null && Classes.TryGetValue(homologyType, out var cls)
                ? cls
                : throw new ArgumentOutOfRangeException(nameof(homologyType), homologyType, "Not a Compara homology_type.");

        /// <summary>The class as written in a table, e.g. "paralog".</summary>
        public static string ClassName(HomologyClass cls) => cls switch
        {
            HomologyClass.Ortholog => "ortholog",
            HomologyClass.Paralog => "paralog",
            HomologyClass.Homoeolog => "homoeolog",
            _ => throw new ArgumentOutOfRangeException(nameof(cls), cls, null)
        };

        private static bool IsNull(string cell) => cell.Length == 0 || cell == "NULL" || cell == "\\N";

        private static double? Real(string cell, string column, string name, int lineNumber)
        {
            if (IsNull(cell))
            {
                return null;
            }

            return double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double v)
                ? v
                : throw new InvalidDataException($"{name} line {lineNumber}: {column} '{cell}' is not a number.");
        }

        private static int? Whole(string cell, string column, string name, int lineNumber)
        {
            if (IsNull(cell))
            {
                return null;
            }

            return int.TryParse(cell, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int v)
                ? v
                : throw new InvalidDataException($"{name} line {lineNumber}: {column} '{cell}' is not a whole number.");
        }

        private static bool? Flag(string cell, string name, int lineNumber) => cell switch
        {
            "1" => true,
            "0" => false,
            _ when IsNull(cell) => null,
            _ => throw new InvalidDataException(
                $"{name} line {lineNumber}: is_high_confidence '{cell}' is neither 0, 1 nor NULL.")
        };
    }
}
