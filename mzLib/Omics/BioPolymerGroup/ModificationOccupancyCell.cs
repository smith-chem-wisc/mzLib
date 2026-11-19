using System.Globalization;
using System.Text.RegularExpressions;

namespace Omics.BioPolymerGroup;

/// <summary>
/// One site read back from an occupancy cell: the inverse of
/// <see cref="SiteSpecificModificationOccupancy.ToModInfoString"/>.
/// </summary>
/// <param name="Position">
/// The position as written: zero-based residue numbering in which 0 is the N-terminus, 1 to length
/// are side chains, and length + 1 is the C-terminus. The C-terminus cannot be recognised without the
/// sequence length, so only <see cref="IsNTerminus"/> is offered.
/// </param>
/// <param name="ModificationIdWithMotif">The modification, e.g. "Oxidation on M". May contain commas,
/// brackets and parentheses ("N6,N6,N6-trimethyllysine on K").</param>
/// <param name="Fraction">The fraction as printed: two decimals for counts, four for intensities.</param>
/// <param name="Numerator">Modified count, or modified intensity to four significant digits.</param>
/// <param name="Denominator">Total count, or total intensity to four significant digits.</param>
public sealed record OccupancySite(int Position, string ModificationIdWithMotif, double Fraction,
    double Numerator, double Denominator)
{
    /// <summary>Position 0 is the N-terminus, not a residue.</summary>
    public bool IsNTerminus => Position == 0;
}

/// <summary>
/// A <c>CountOccupancy_</c> or <c>IntensityOccupancy_</c> cell of a MetaMorpheus protein-group or
/// peptide table, parsed. The format is written by
/// <see cref="SampleGroupResult.FormatOccupancy"/> and <see cref="SiteSpecificModificationOccupancy.ToModInfoString"/>,
/// and this type is its inverse: sites joined by <c>;</c> within an entity, entities joined by <c>|</c>.
/// </summary>
/// <remarks>
/// <para>An entity is a member of the group (protein level) or a base sequence (peptide level).
/// The cell does not name it, and an entity with no sites is dropped, so the entities cannot be
/// zipped with the accession list by position.</para>
/// <para>Count and intensity cells round differently, so trust the pair in a count cell and the
/// fraction in an intensity cell. Exact intensities are in the <c>Intensity_</c> columns.</para>
/// </remarks>
public sealed class ModificationOccupancyCell
{
    /// <summary>What MetaMorpheus 1.1.x writes in place of any cell over Excel's limit.</summary>
    public const string ExcelTruncationText = "Output too long for Excel";

    // Anchored on ",info:fraction=<number>(" so a modification name may contain commas, brackets and
    // parentheses. The name may not itself span that anchor, so two sites run together without a
    // separator are refused rather than read as one site with a very long name. Nor may it contain a "]"
    // followed by a separator, so a malformed site cannot swallow the one after it.
    private static readonly Regex Site = new(
        @"\Gpos(?<pos>\d+)\[(?<name>(?:(?!,info:fraction=|\][;|]).)+),info:fraction=(?<f>[-+0-9.eE]+|NaN)\((?<n>[^/()]+)/(?<d>[^()]+)\)\](?<sep>[;|]|$)",
        RegexOptions.Compiled | RegexOptions.CultureInvariant);

    // What a cut leaves: the start of one site running to the end of the cell, each part of Site in turn,
    // each part optional, down to a lone "p". The name may not contain a "]" that closes a site (one followed by a
    // separator or the end), so a finished site with a misspelt field is refused rather than read as cut.
    private static readonly Regex CutSite = new(
        @"\G(?:po?|pos(?:\d+(?:\[(?:(?!,info:fraction=|\][;|]|\]$).)*(?:,info:fraction=(?:[-+0-9.eE]*|N|Na|NaN)(?:\([^/()]*(?:/[^()]*\)?)?)?)?)?)?)$",
        RegexOptions.Compiled | RegexOptions.CultureInvariant);

    private ModificationOccupancyCell(IReadOnlyList<IReadOnlyList<OccupancySite>> entities, bool truncated)
    {
        Entities = entities;
        IsTruncated = truncated;
    }

    /// <summary>An empty cell: no modified sites were reported.</summary>
    public static ModificationOccupancyCell Empty { get; } = new([], false);

    /// <summary>Sites grouped by entity, in the order written. Empty entities are not represented.</summary>
    public IReadOnlyList<IReadOnlyList<OccupancySite>> Entities { get; }

    /// <summary>Every site, across entities.</summary>
    public IEnumerable<OccupancySite> Sites => Entities.SelectMany(e => e);

    /// <summary>
    /// The cell was cut short: either replaced wholesale by <see cref="ExcelTruncationText"/>, or cut
    /// mid-site (mzLib's writer cuts at <see cref="BioPolymerGroupTsvSchema.MaxStringLength"/>, whatever
    /// it was set to when the file was written). Complete sites before the cut are kept. A truncated
    /// cell is not an empty one.
    /// </summary>
    public bool IsTruncated { get; }

    /// <summary>
    /// Parses one cell. Blank means no sites. Throws <see cref="FormatException"/> for text that is not
    /// an occupancy cell at all, rather than returning a partial reading of it.
    /// </summary>
    public static ModificationOccupancyCell Parse(string? cell)
    {
        if (string.IsNullOrWhiteSpace(cell))
            return Empty;
        if (cell == ExcelTruncationText)
            return new([], true);

        var entities = new List<IReadOnlyList<OccupancySite>>();
        var current = new List<OccupancySite>();
        int at = 0;
        while (at < cell.Length)
        {
            var m = Site.Match(cell, at);
            if (!m.Success)
            {
                // A cut leaves an unfinished site at the end. It is recognised by its shape, not by the
                // cell's length: the limit it was cut at is a writer setting this reader cannot know.
                // Anything else is not this format, and guessing at it would be the silent error this
                // type exists to stop.
                if (CutSite.IsMatch(cell, at))
                {
                    if (current.Count > 0) entities.Add(current);
                    return new(entities, true);
                }
                throw new FormatException($"Not an occupancy cell at character {at}: \"{Excerpt(cell, at)}\"");
            }
            current.Add(new OccupancySite(
                int.Parse(m.Groups["pos"].Value, CultureInfo.InvariantCulture),
                m.Groups["name"].Value,
                ParseNumber(m.Groups["f"].Value),
                ParseNumber(m.Groups["n"].Value),
                ParseNumber(m.Groups["d"].Value)));
            string sep = m.Groups["sep"].Value;
            if (sep != ";")
            {
                entities.Add(current);
                current = [];
            }
            at = m.Index + m.Length;
        }
        return new(entities, false);
    }

    private static double ParseNumber(string s) =>
        double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static string Excerpt(string cell, int at) =>
        cell.Substring(at, Math.Min(40, cell.Length - at));
}
