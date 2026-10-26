using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace Omics.SequenceAlignment;

/// <summary>
/// A residue substitution matrix in NCBI's text format (the format of
/// <c>ftp.ncbi.nlm.nih.gov/blast/matrices/</c>): <c>#</c> comment lines, one header line of
/// single-letter residue codes, then one row per code.
/// </summary>
/// <remarks>
/// <para>
/// A residue the matrix does not name (selenocysteine <c>U</c>, pyrrolysine <c>O</c>, <c>J</c>,
/// or any other letter) is scored as <c>X</c>, the matrix's unknown residue. Letters are compared
/// case-insensitively. A matrix with no <c>X</c> refuses an unnamed residue instead.
/// </para>
/// <para>
/// The reader refuses rather than guesses: a row whose length differs from the header, a row
/// label that is not the next header code, a code that appears twice, a score that is not an
/// integer, or a matrix that is not symmetric each throw <see cref="InvalidDataException"/>.
/// </para>
/// </remarks>
public sealed class SubstitutionMatrix
{
    private static readonly Lazy<SubstitutionMatrix> LazyBlosum62 = new(() =>
    {
        var assembly = Assembly.GetExecutingAssembly();
        string resource = $"{assembly.GetName().Name}.Resources.BLOSUM62.txt";
        using var stream = assembly.GetManifestResourceStream(resource)
            ?? throw new InvalidOperationException($"Embedded resource {resource} is missing.");
        return Parse("BLOSUM62", new StreamReader(stream));
    });

    /// <summary>
    /// NCBI's BLOSUM62, verbatim. It is the matrix BLASTP uses by default, with a gap open of 11
    /// and a gap extension of 1.
    /// </summary>
    public static SubstitutionMatrix Blosum62 => LazyBlosum62.Value;

    // Indexed by upper-case ASCII letter or '*'; -1 means the matrix does not name it.
    private readonly int[] _indexByChar = Enumerable.Repeat(-1, 128).ToArray();
    private readonly int[,] _scores;
    private readonly int _unknownIndex;

    private SubstitutionMatrix(string name, IReadOnlyList<char> codes, int[,] scores)
    {
        Name = name;
        Codes = codes;
        _scores = scores;
        for (int k = 0; k < codes.Count; k++)
        {
            _indexByChar[codes[k]] = k;
        }
        _unknownIndex = _indexByChar['X'];
        IdenticalSequencesAlignToThemselves = ComputeIdentityIsOptimal();
    }

    /// <summary>The matrix's name, as it appears in an aligner id (<see cref="PairwiseAligner.Id"/>).</summary>
    public string Name { get; }

    /// <summary>The residue codes the matrix names, in the header's order.</summary>
    public IReadOnlyList<char> Codes { get; }

    /// <summary>
    /// True when every residue scores at least 0 against itself and, for every pair,
    /// <c>Score(a, b) &lt;= (Score(a, a) + Score(b, b)) / 2</c>, ignoring <c>X</c> and <c>*</c>.
    /// When it holds, two identical sequences free of <c>X</c> and <c>*</c> score highest aligned
    /// residue for residue, whatever the gap penalties: every other alignment pairs each residue at
    /// most once, at no more than the mean of the two self-scores, and leaves the rest unpaired.
    /// <see cref="PairwiseAligner"/> skips the dynamic programme only in that case.
    /// </summary>
    public bool IdenticalSequencesAlignToThemselves { get; }

    /// <summary>The score of aligning residue <paramref name="a"/> with residue <paramref name="b"/>.</summary>
    /// <exception cref="ArgumentException">A residue the matrix does not name, when it has no <c>X</c>.</exception>
    public int Score(char a, char b) => _scores[IndexOf(a), IndexOf(b)];

    internal int IndexOf(char c)
    {
        char upper = char.ToUpperInvariant(c);
        int index = upper < 128 ? _indexByChar[upper] : -1;
        if (index >= 0)
        {
            return index;
        }
        if (_unknownIndex >= 0)
        {
            return _unknownIndex;
        }
        throw new ArgumentException($"{Name} does not score residue '{c}' and has no X to score it as.");
    }

    internal int ScoreByIndex(int a, int b) => _scores[a, b];

    /// <summary>Reads a matrix in NCBI's text format.</summary>
    /// <param name="name">The name the matrix is known by, carried into <see cref="PairwiseAligner.Id"/>.</param>
    /// <param name="reader">The text. It is read to the end and not disposed.</param>
    /// <exception cref="InvalidDataException">The text is not a square, symmetric matrix of integers.</exception>
    public static SubstitutionMatrix Parse(string name, TextReader reader)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentNullException.ThrowIfNull(reader);

        List<char>? codes = null;
        var rows = new List<int[]>();
        string? line;
        int lineNumber = 0;
        while ((line = reader.ReadLine()) != null)
        {
            lineNumber++;
            string trimmed = line.Trim();
            if (trimmed.Length == 0 || trimmed.StartsWith('#'))
            {
                continue;
            }
            string[] cells = trimmed.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (codes == null)
            {
                codes = new List<char>();
                foreach (string cell in cells)
                {
                    char code = char.ToUpperInvariant(cell[0]);
                    if (cell.Length != 1 || code >= 128 || codes.Contains(code))
                    {
                        throw new InvalidDataException($"{name} line {lineNumber}: residue code '{cell}' is not a single new letter.");
                    }
                    codes.Add(code);
                }
                continue;
            }
            if (rows.Count == codes.Count)
            {
                throw new InvalidDataException($"{name} line {lineNumber}: more rows than the {codes.Count} codes in the header.");
            }
            if (cells.Length != codes.Count + 1)
            {
                throw new InvalidDataException($"{name} line {lineNumber}: {cells.Length - 1} scores, expected {codes.Count}.");
            }
            if (cells[0].Length != 1 || char.ToUpperInvariant(cells[0][0]) != codes[rows.Count])
            {
                throw new InvalidDataException($"{name} line {lineNumber}: row '{cells[0]}', expected '{codes[rows.Count]}'.");
            }
            var row = new int[codes.Count];
            for (int k = 0; k < codes.Count; k++)
            {
                if (!int.TryParse(cells[k + 1], out row[k]))
                {
                    throw new InvalidDataException($"{name} line {lineNumber}: score '{cells[k + 1]}' is not an integer.");
                }
            }
            rows.Add(row);
        }

        if (codes == null || rows.Count != codes.Count)
        {
            throw new InvalidDataException($"{name}: {rows.Count} rows, expected {codes?.Count ?? 0}.");
        }
        var scores = new int[codes.Count, codes.Count];
        for (int r = 0; r < codes.Count; r++)
        {
            for (int c = 0; c < codes.Count; c++)
            {
                scores[r, c] = rows[r][c];
            }
        }
        for (int r = 0; r < codes.Count; r++)
        {
            for (int c = r + 1; c < codes.Count; c++)
            {
                if (scores[r, c] != scores[c, r])
                {
                    throw new InvalidDataException($"{name}: not symmetric, {codes[r]}/{codes[c]} is {scores[r, c]} but {codes[c]}/{codes[r]} is {scores[c, r]}.");
                }
            }
        }
        return new SubstitutionMatrix(name, codes.AsReadOnly(), scores);
    }

    private bool ComputeIdentityIsOptimal()
    {
        var residues = Enumerable.Range(0, Codes.Count).Where(k => Codes[k] != 'X' && Codes[k] != '*').ToList();
        foreach (int a in residues)
        {
            if (_scores[a, a] < 0)
            {
                return false;
            }
            foreach (int b in residues)
            {
                if (2 * _scores[a, b] > _scores[a, a] + _scores[b, b])
                {
                    return false;
                }
            }
        }
        return true;
    }

    internal bool IsUnknownOrStop(char c)
    {
        char upper = char.ToUpperInvariant(c);
        return upper == 'X' || upper == '*' || upper >= 128 || _indexByChar[upper] < 0;
    }
}
