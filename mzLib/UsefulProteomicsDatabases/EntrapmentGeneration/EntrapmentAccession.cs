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
