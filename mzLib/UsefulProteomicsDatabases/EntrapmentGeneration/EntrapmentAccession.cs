#nullable enable
using System;
using System.Globalization;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>
/// The accession an entrapment protein wears, and how to read the target and fold back out of it.
/// </summary>
/// <remarks>
/// <para>Two existing conventions are followed rather than a third invented. The entrapment marker
/// goes <b>first</b>, because <see cref="ProteinDbLoader"/> decides a protein is entrapment by
/// looking for the identifier anywhere in its accession, and prepends it when absent -- so an
/// accession shaped this way reloads as entrapment with no extra argument and is never
/// double-prefixed. The fold goes <b>last</b>, matching
/// <c>VariantApplication.GetAccession</c>, which appends derived detail as a suffix.</para>
/// <para>Everything needed to pair a discovery back to its target therefore travels in the
/// accession, which a search already reports. Nothing has to be carried alongside the database, and
/// nothing leaks into the sequence -- which matters, because target and partner are deliberately
/// identical in mass and composition and must stay that way.</para>
/// </remarks>
public static class EntrapmentAccession
{
    private const string FoldMarker = "_f";

    /// <summary>
    /// Marks an entry taken from a foreign proteome rather than rearranged from a target.
    /// </summary>
    /// <remarks>
    /// It sits immediately after the identifier, so a foreign entry reads
    /// <c>Random_foreign_P12345</c>. Two things had to be true at once: a consumer scanning for the
    /// entrapment identifier must still find it -- miss one and a foreign entry is counted as a
    /// <i>target</i>, which corrupts the estimate in the worst direction -- and
    /// <see cref="TryParse"/> must refuse it, because there is no target to name.
    /// </remarks>
    private const string ForeignMarker = "foreign_";

    /// <summary>The accession for one fold of one target's entrapment partner.</summary>
    /// <exception cref="MzLibUtil.MzLibException">The target accession is missing, or the fold is negative.</exception>
    public static string Format(string targetAccession, int fold,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (string.IsNullOrEmpty(targetAccession))
        {
            throw new MzLibUtil.MzLibException("An entrapment accession needs a target accession.");
        }
        if (fold < 0)
        {
            throw new MzLibUtil.MzLibException($"Fold must not be negative, but was {fold}.");
        }

        return $"{entrapmentIdentifier}_{targetAccession}{FoldMarker}{fold.ToString(CultureInfo.InvariantCulture)}";
    }

    /// <summary>The accession for an entry taken from a foreign proteome.</summary>
    /// <param name="foreignAccession">The protein's accession in its own database.</param>
    /// <exception cref="MzLibUtil.MzLibException">The foreign accession is missing.</exception>
    public static string FormatForeign(string foreignAccession,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        if (string.IsNullOrEmpty(foreignAccession))
        {
            throw new MzLibUtil.MzLibException("A foreign entrapment accession needs an accession.");
        }

        return $"{entrapmentIdentifier}_{ForeignMarker}{foreignAccession}";
    }

    /// <summary>Reads a foreign entry's own accession back out.</summary>
    /// <returns>False when this is not a foreign entrapment accession.</returns>
    public static bool TryParseForeign(string? accession, out string foreignAccession,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        foreignAccession = string.Empty;
        string prefix = entrapmentIdentifier + "_" + ForeignMarker;
        if (string.IsNullOrEmpty(accession)
            || !accession.StartsWith(prefix, StringComparison.OrdinalIgnoreCase))
        {
            return false;
        }

        foreignAccession = accession.Substring(prefix.Length);
        return foreignAccession.Length > 0;
    }

    /// <summary>
    /// Reads the target accession and fold back out, or reports that this is not one of ours.
    /// </summary>
    /// <remarks>
    /// The fold is taken from the <em>last</em> "_f&lt;digits&gt;" rather than by splitting on
    /// underscores: UniProt entry names contain them, so a naive split loses the target.
    /// </remarks>
    public static bool TryParse(string? accession, out string targetAccession, out int fold,
        string entrapmentIdentifier = ProteinDbLoader.DefaultEntrapmentIdentifier)
    {
        targetAccession = string.Empty;
        fold = 0;

        string prefix = entrapmentIdentifier + "_";
        if (string.IsNullOrEmpty(accession) || !accession.StartsWith(prefix, StringComparison.OrdinalIgnoreCase))
        {
            return false;
        }

        // A foreign entry has no target, and its own accession may well end in something that looks
        // like a fold marker -- "Random_foreign_ABC_f1" would otherwise parse as target "foreign_ABC".
        // Refusing on the marker rather than on the shape keeps that from ever being a guess.
        if (accession.StartsWith(prefix + ForeignMarker, StringComparison.OrdinalIgnoreCase))
        {
            return false;
        }

        int marker = accession.LastIndexOf(FoldMarker, StringComparison.Ordinal);
        if (marker <= prefix.Length - 1)
        {
            return false;
        }

        string digits = accession.Substring(marker + FoldMarker.Length);
        if (digits.Length == 0
            || !int.TryParse(digits, NumberStyles.None, CultureInfo.InvariantCulture, out int parsedFold))
        {
            return false;
        }

        targetAccession = accession.Substring(prefix.Length, marker - prefix.Length);
        if (targetAccession.Length == 0)
        {
            return false;
        }

        fold = parsedFold;
        return true;
    }
}
