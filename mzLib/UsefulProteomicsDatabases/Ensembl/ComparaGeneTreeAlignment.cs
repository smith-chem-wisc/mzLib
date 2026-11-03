using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>A protein's row in one of Compara's gene-tree alignments.</summary>
    /// <param name="ProteinId">The Ensembl protein stable id, as the alignment's header gives it.</param>
    /// <param name="AlignmentIndex">Which alignment holds it: 0 for the first in the file, and so on.</param>
    /// <param name="TreeId">The gene tree its gene belongs to, from the <see cref="ComparaGeneTreeContent"/> the alignment was restricted to.</param>
    /// <param name="Sequence">The protein's residues, with the alignment's gaps removed.</param>
    public sealed record ComparaAlignedProtein(string ProteinId, int AlignmentIndex, string TreeId, string Sequence);

    /// <summary>What one alignment column says about a residue in another protein.</summary>
    public enum ComparaColumnOutcome
    {
        /// <summary>The target protein has a residue in the source residue's column.</summary>
        Aligned,

        /// <summary>The target protein has a gap in that column. No residue corresponds.</summary>
        GapInTarget,

        /// <summary>The two proteins are in different alignments (different gene trees). The alignment says nothing about them.</summary>
        DifferentAlignments,

        /// <summary>The source protein was not kept from the file: not in it, or outside the restriction.</summary>
        SourceNotInAlignment,

        /// <summary>The target protein was not kept from the file.</summary>
        TargetNotInAlignment,
    }

    /// <summary>The answer to <see cref="ComparaGeneTreeAlignment.MapResidue"/>.</summary>
    /// <param name="Outcome">Exactly one outcome.</param>
    /// <param name="Column">The source residue's 1-based alignment column; null unless both proteins share an alignment.</param>
    /// <param name="TargetPosition">The 1-based position in the target protein; null unless <see cref="ComparaColumnOutcome.Aligned"/>.</param>
    /// <param name="TargetResidue">The target's residue there; null unless <see cref="ComparaColumnOutcome.Aligned"/>.</param>
    public sealed record ComparaColumnMapping(ComparaColumnOutcome Outcome, int? Column, int? TargetPosition, char? TargetResidue);

    /// <summary>
    /// Compara's peptide alignments of its gene trees
    /// (emf/ensembl-compara/homologies/Compara.&lt;release&gt;.protein_&lt;collection&gt;.aa.fasta.gz): one gapped
    /// multiple alignment per tree, as FASTA, the alignments separated by a "//" line. Each header is a bare
    /// Ensembl protein id.
    ///
    /// It answers which residue of one protein sits in the same column as a residue of another protein of
    /// the same tree. That is the step that crosses a homology edge at residue grain; reaching it from a
    /// UniProt position needs a pairwise alignment to the Ensembl protein first.
    ///
    /// The file is large (908 MB compressed in release 116, every species in the collection), so
    /// <see cref="Load"/> keeps only the proteins it is asked for: the canonical proteins of a
    /// <see cref="ComparaGeneTreeContent"/>. Every row is checked whether it is kept or not.
    ///
    /// The reader refuses rather than guesses. Each of these throws:
    /// a sequence line before any header; an empty header; a character that is not a letter, '*' or '-';
    /// an alignment whose rows differ in width; an empty alignment; a kept protein that appears twice; and
    /// an alignment whose kept proteins belong to more than one gene tree.
    /// </summary>
    public sealed class ComparaGeneTreeAlignment
    {
        private static readonly Regex ReleaseInFileName =
            new(@"^Compara\.(\d+)\.", RegexOptions.Compiled);

        private static readonly Regex CollectionInFileName =
            new(@"^Compara\.\d+\.([A-Za-z_]+)\.aa\.fasta", RegexOptions.Compiled);

        private readonly Dictionary<string, (ComparaAlignedProtein Protein, int[] Columns)> _proteins;

        private ComparaGeneTreeAlignment(Dictionary<string, (ComparaAlignedProtein, int[])> proteins,
            string sourceFileName, string sourceSha256, string release, string collection, int alignmentCount,
            string restrictedTo)
        {
            _proteins = proteins;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
            Collection = collection;
            AlignmentCount = alignmentCount;
            RestrictedToGeneTreeContent = restrictedTo;
        }

        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read (the compressed bytes for a .gz).</summary>
        public string SourceSha256 { get; }

        /// <summary>The Ensembl release parsed from the file name, or null when it carries none.</summary>
        public string Release { get; }

        /// <summary>The collection parsed from the file name (e.g. "protein_default"), or null.</summary>
        public string Collection { get; }

        /// <summary>The number of alignments in the file, kept proteins or not.</summary>
        public int AlignmentCount { get; }

        /// <summary>The sha256 of the gene-tree content whose canonical proteins were kept.</summary>
        public string RestrictedToGeneTreeContent { get; }

        /// <summary>The number of proteins kept.</summary>
        public int Count => _proteins.Count;

        /// <summary>Every kept protein, in ordinal order of protein id.</summary>
        public IEnumerable<ComparaAlignedProtein> Proteins =>
            _proteins.Values.Select(p => p.Protein).OrderBy(p => p.ProteinId, StringComparer.Ordinal);

        public bool TryGetProtein(string proteinId, out ComparaAlignedProtein protein)
        {
            protein = null;
            if (proteinId == null || !_proteins.TryGetValue(proteinId, out var entry))
            {
                return false;
            }
            protein = entry.Protein;
            return true;
        }

        /// <summary>
        /// The residue of <paramref name="targetProteinId"/> in the same alignment column as residue
        /// <paramref name="oneBasedPosition"/> of <paramref name="sourceProteinId"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">A protein id is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The position is outside the source protein.</exception>
        public ComparaColumnMapping MapResidue(string sourceProteinId, int oneBasedPosition, string targetProteinId)
        {
            ArgumentNullException.ThrowIfNull(sourceProteinId);
            ArgumentNullException.ThrowIfNull(targetProteinId);

            if (!_proteins.TryGetValue(sourceProteinId, out var source))
            {
                return new ComparaColumnMapping(ComparaColumnOutcome.SourceNotInAlignment, null, null, null);
            }
            if (oneBasedPosition < 1 || oneBasedPosition > source.Columns.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(oneBasedPosition), oneBasedPosition,
                    $"{sourceProteinId} has {source.Columns.Length} residues.");
            }
            if (!_proteins.TryGetValue(targetProteinId, out var target))
            {
                return new ComparaColumnMapping(ComparaColumnOutcome.TargetNotInAlignment, null, null, null);
            }
            if (source.Protein.AlignmentIndex != target.Protein.AlignmentIndex)
            {
                return new ComparaColumnMapping(ComparaColumnOutcome.DifferentAlignments, null, null, null);
            }

            int column = source.Columns[oneBasedPosition - 1];
            int index = Array.BinarySearch(target.Columns, column);
            return index < 0
                ? new ComparaColumnMapping(ComparaColumnOutcome.GapInTarget, column, null, null)
                : new ComparaColumnMapping(ComparaColumnOutcome.Aligned, column, index + 1, target.Protein.Sequence[index]);
        }

        /// <summary>
        /// Reads the alignments, keeping the canonical protein of every gene in <paramref name="restrictTo"/>.
        /// Every row's characters and width are checked, kept or not.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="restrictTo"/> is null.</exception>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file breaks one of the rules in the class summary.</exception>
        public static ComparaGeneTreeAlignment Load(string path, ComparaGeneTreeContent restrictTo)
        {
            ArgumentNullException.ThrowIfNull(restrictTo);

            string sha256 = EnsemblFile.Sha256(path, "Compara gene-tree alignment");
            string name = Path.GetFileName(path);
            var treeByProtein = restrictTo.Members.ToDictionary(m => m.CanonicalProteinId, m => m.TreeId, StringComparer.Ordinal);
            var kept = new Dictionary<string, (ComparaAlignedProtein, int[])>(StringComparer.Ordinal);

            int alignmentIndex = 0, rowsInAlignment = 0, width = -1;
            string treeOfAlignment = null, keptInAlignment = null;
            string protein = null;
            bool keep = false;
            int rowWidth = 0;
            StringBuilder residues = new();
            List<int> columns = new();

            using (var reader = EnsemblFile.OpenText(path))
            {
                int lineNumber = 0;
                string line;

                void EndRow()
                {
                    if (protein == null)
                    {
                        return;
                    }
                    if (width < 0)
                    {
                        width = rowWidth;
                    }
                    else if (rowWidth != width)
                    {
                        throw new InvalidDataException(
                            $"{name} line {lineNumber}: {protein} is {rowWidth} columns wide, the alignment {width}.");
                    }
                    if (keep)
                    {
                        string tree = treeByProtein[protein];
                        if (treeOfAlignment != null && treeOfAlignment != tree)
                        {
                            throw new InvalidDataException(
                                $"{name} line {lineNumber}: alignment {alignmentIndex} holds {keptInAlignment} of {treeOfAlignment} and {protein} of {tree}.");
                        }
                        treeOfAlignment = tree;
                        keptInAlignment = protein;
                        var row = new ComparaAlignedProtein(protein, alignmentIndex, tree, residues.ToString());
                        if (!kept.TryAdd(protein, (row, columns.ToArray())))
                        {
                            throw new InvalidDataException($"{name} line {lineNumber}: {protein} appears twice.");
                        }
                    }
                    rowsInAlignment++;
                    protein = null;
                }

                void EndAlignment()
                {
                    EndRow();
                    if (rowsInAlignment == 0)
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: an alignment with no rows.");
                    }
                    alignmentIndex++;
                    rowsInAlignment = 0;
                    width = -1;
                    treeOfAlignment = null;
                    keptInAlignment = null;
                }

                while ((line = reader.ReadLine()) != null)
                {
                    lineNumber++;
                    string text = line.Trim();
                    if (text.Length == 0)
                    {
                        continue;
                    }
                    if (text == "//")
                    {
                        EndAlignment();
                        continue;
                    }
                    if (text[0] == '>')
                    {
                        EndRow();
                        string[] header = text.Substring(1).Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                        if (header.Length == 0)
                        {
                            throw new InvalidDataException($"{name} line {lineNumber}: a header with no protein id.");
                        }
                        protein = header[0];
                        keep = treeByProtein.ContainsKey(protein);
                        rowWidth = 0;
                        residues.Clear();
                        columns.Clear();
                        continue;
                    }
                    if (protein == null)
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: a sequence line before any header.");
                    }
                    foreach (char c in text)
                    {
                        rowWidth++;
                        if (c == '-')
                        {
                            continue;
                        }
                        if (!char.IsAsciiLetter(c) && c != '*')
                        {
                            throw new InvalidDataException($"{name} line {lineNumber}: character '{c}' in {protein} is not a residue or '-'.");
                        }
                        if (keep)
                        {
                            residues.Append(c);
                            columns.Add(rowWidth);
                        }
                    }
                }

                // The README separates alignments with "//"; tolerate a file that does not end with one.
                if (protein != null || rowsInAlignment > 0)
                {
                    EndAlignment();
                }
            }

            var release = ReleaseInFileName.Match(name);
            var collection = CollectionInFileName.Match(name);
            return new ComparaGeneTreeAlignment(kept, name, sha256,
                release.Success ? release.Groups[1].Value : null,
                collection.Success ? collection.Groups[1].Value : null,
                alignmentIndex, restrictTo.SourceSha256);
        }
    }
}
