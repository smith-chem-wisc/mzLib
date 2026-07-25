#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;

namespace TopDownSimulator.Ms2;

/// <summary>
/// Amino-acid-pair peptide-bond cleavage propensities from Haverland et al., "Defining Gas-Phase
/// Fragmentation Propensities of Intact Proteins During Native Top-Down Mass Spectrometry"
/// (J Am Soc Mass Spectrom), Supporting Information worksheet "Data".
/// </summary>
/// <remarks>
/// <para>
/// The source worksheet is long-format: 400 rows (a complete ordered 20 x 20 residue grid) by 15
/// columns. Only the denatured top-down column, "Propensity (dTD)", is kept: all 400 pairs are
/// defined and non-zero there, so a cascade over any sequence is well behaved. The paper's native
/// top-down column is not carried — it has 83 exact zeros and one unsampled pair.
/// </para>
/// <para>
/// The value is already <c>matched fragmentation events / possible fragmentation events</c>, i.e. a
/// dimensionless fraction in [0, 1] — so <b>no normalization is applied on load</b>. The published
/// values are carried verbatim into <c>Ms2/ResiduePairCleavagePropensities.tsv</c>, which is
/// embedded in this assembly and parsed on first use.
/// </para>
/// <para>
/// Any rescaling you want is an explicit knob on <see cref="ResiduePairCleavagePropensityModel"/>
/// (<c>PropensityScale</c>), not a hidden constant in here.
/// </para>
/// </remarks>
public sealed class ResiduePairCleavagePropensityTable
{
    /// <summary>Resource name of the embedded copy of the table.</summary>
    public const string EmbeddedResourceName = "TopDownSimulator.Ms2.ResiduePairCleavagePropensities.tsv";

    private static readonly Lazy<ResiduePairCleavagePropensityTable> LazyDefault =
        new(LoadEmbedded, isThreadSafe: true);

    private readonly Dictionary<(char N, char C), double> _propensities;

    private ResiduePairCleavagePropensityTable(Dictionary<(char, char), double> propensities)
    {
        _propensities = propensities;
        Residues = propensities.Keys
            .SelectMany(k => new[] { k.Item1, k.Item2 })
            .Distinct()
            .OrderBy(c => c)
            .ToArray();
    }

    /// <summary>The embedded published table. Parsed once, then cached.</summary>
    public static ResiduePairCleavagePropensityTable Default => LazyDefault.Value;

    /// <summary>Number of residue pairs in the table. The published table has 400.</summary>
    public int PairCount => _propensities.Count;

    /// <summary>Distinct residue letters appearing in the table, ascending.</summary>
    public IReadOnlyList<char> Residues { get; }

    /// <summary>Reads the embedded copy of the published table.</summary>
    public static ResiduePairCleavagePropensityTable LoadEmbedded()
    {
        var assembly = typeof(ResiduePairCleavagePropensityTable).GetTypeInfo().Assembly;
        using var stream = assembly.GetManifestResourceStream(EmbeddedResourceName)
            ?? throw new InvalidOperationException(
                $"Embedded resource '{EmbeddedResourceName}' was not found in {assembly.GetName().Name}.");
        using var reader = new StreamReader(stream);
        return Parse(ReadLines(reader));
    }

    /// <summary>Reads a table from a TSV file on disk, for swapping in a different propensity set.</summary>
    public static ResiduePairCleavagePropensityTable Load(string path)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("A path is required.", nameof(path));

        return Parse(File.ReadLines(path));
    }

    /// <summary>
    /// Parses TSV lines. Blank lines and lines starting with '#' are skipped; the first surviving
    /// line is the header. Required columns are <c>NTerminalResidue</c>,
    /// <c>CTerminalResidue</c> and <c>DenaturedTopDownPropensity</c>. A value of <c>NaN</c> means
    /// "not measured".
    /// </summary>
    public static ResiduePairCleavagePropensityTable Parse(IEnumerable<string> lines)
    {
        if (lines is null)
            throw new ArgumentNullException(nameof(lines));

        var table = new Dictionary<(char, char), double>();
        string[]? header = null;
        int nCol = -1, cCol = -1, propensityCol = -1;

        foreach (string raw in lines)
        {
            if (string.IsNullOrWhiteSpace(raw) || raw[0] == '#')
                continue;

            var fields = raw.Split('\t');

            if (header is null)
            {
                header = fields;
                nCol = Array.IndexOf(header, "NTerminalResidue");
                cCol = Array.IndexOf(header, "CTerminalResidue");
                propensityCol = Array.IndexOf(header, "DenaturedTopDownPropensity");

                if (nCol < 0 || cCol < 0 || propensityCol < 0)
                    throw new FormatException(
                        "Propensity table header must contain NTerminalResidue, CTerminalResidue " +
                        "and DenaturedTopDownPropensity.");
                continue;
            }

            if (fields.Length <= Math.Max(Math.Max(nCol, cCol), propensityCol))
                throw new FormatException($"Propensity table row has too few columns: '{raw}'.");

            char nResidue = SingleResidue(fields[nCol]);
            char cResidue = SingleResidue(fields[cCol]);

            table[(nResidue, cResidue)] = ParsePropensity(fields[propensityCol]);
        }

        if (header is null || table.Count == 0)
            throw new FormatException("Propensity table contained no data rows.");

        return new ResiduePairCleavagePropensityTable(table);
    }

    /// <summary>
    /// Looks up the propensity for the bond whose N-terminal side is
    /// <paramref name="nTerminalResidue"/> and C-terminal side is <paramref name="cTerminalResidue"/>.
    /// Returns false when the pair is absent from the table or its value was not measured (NaN).
    /// </summary>
    public bool TryGetPropensity(char nTerminalResidue, char cTerminalResidue, out double propensity)
    {
        propensity = double.NaN;
        if (!_propensities.TryGetValue((nTerminalResidue, cTerminalResidue), out double value))
            return false;

        if (double.IsNaN(value))
            return false;

        propensity = value;
        return true;
    }

    /// <summary>Looks up a propensity, throwing when the pair is absent or unmeasured.</summary>
    public double GetPropensity(char nTerminalResidue, char cTerminalResidue)
    {
        if (!TryGetPropensity(nTerminalResidue, cTerminalResidue, out double propensity))
            throw new KeyNotFoundException(
                $"No cleavage propensity for residue pair {nTerminalResidue}|{cTerminalResidue}.");

        return propensity;
    }

    /// <summary>True when the table has a measured value for the pair.</summary>
    public bool Contains(char nTerminalResidue, char cTerminalResidue) =>
        TryGetPropensity(nTerminalResidue, cTerminalResidue, out _);

    private static IEnumerable<string> ReadLines(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) is not null)
            yield return line;
    }

    private static char SingleResidue(string field)
    {
        string trimmed = field.Trim();
        if (trimmed.Length != 1)
            throw new FormatException($"Expected a single residue letter but found '{field}'.");

        return trimmed[0];
    }

    private static double ParsePropensity(string field)
    {
        string trimmed = field.Trim();
        if (trimmed.Length == 0 ||
            trimmed.Equals("NA", StringComparison.OrdinalIgnoreCase) ||
            trimmed.Equals("NaN", StringComparison.OrdinalIgnoreCase))
        {
            return double.NaN;
        }

        if (!double.TryParse(trimmed, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
            throw new FormatException($"Could not parse propensity '{field}'.");

        if (value < 0 || value > 1)
            throw new FormatException(
                $"Propensity '{field}' is outside [0, 1]; the table is expected to hold fractions, " +
                "not percentages or arbitrary scores.");

        return value;
    }
}
