using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Security.Cryptography;
using System.Text;

namespace Quantification.Differential;

/// <summary>
/// The pieces of the observation tables' canonical text (STAT1 M5b ruling 10): one text per table that is the same on any
/// machine and culture and differs whenever anything the analysis reads differs, so its sha256 can stand for the table in
/// the analysis id. Fields are joined by tabs and lines by <c>\n</c>; text fields are escaped so no field can hold either.
/// </summary>
internal static class CanonicalText
{
    /// <summary>A number in round-trip invariant form; NaN as <c>NaN</c>.</summary>
    public static string Number(double value) =>
        double.IsNaN(value) ? "NaN" : value.ToString("R", CultureInfo.InvariantCulture);

    /// <summary>An integer in invariant form.</summary>
    public static string Number(int value) => value.ToString(CultureInfo.InvariantCulture);

    /// <summary>A text field with <c>\</c>, tab, line feed and carriage return escaped; null as <c>\N</c>, which no text gives.</summary>
    public static string Escape(string? text) =>
        text is null ? @"\N" : text.Replace(@"\", @"\\").Replace("\t", @"\t").Replace("\n", @"\n").Replace("\r", @"\r");

    /// <summary>The lower-case hex sha256 of the text's UTF-8 bytes.</summary>
    public static string Sha256(string text) =>
        Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(text))).ToLowerInvariant();

    /// <summary>One line per described protein group, by name (ordinal): name, q-value, decoy, contaminant, entrapment, genes.</summary>
    public static void AppendGroups(StringBuilder text, IReadOnlyDictionary<string, ProteinGroupInfo> groups)
    {
        foreach (var name in groups.Keys.Order(StringComparer.Ordinal))
        {
            var g = groups[name];
            text.Append("group\t").Append(Escape(name)).Append('\t').Append(Number(g.QValue))
                .Append('\t').Append(g.IsDecoy ? '1' : '0').Append('\t').Append(g.IsContaminant ? '1' : '0')
                .Append('\t').Append(g.IsEntrapment ? '1' : '0').Append('\t').Append(Escape(g.Genes)).Append('\n');
        }
    }

    /// <summary>A peptide's identity: full sequence, base sequence and its groups (in the peptide's order, <c>;</c>-joined).</summary>
    public static void AppendPeptide(StringBuilder text, ObservedPeptide peptide) =>
        text.Append(Escape(peptide.FullSequence)).Append('\t').Append(Escape(peptide.BaseSequence)).Append('\t')
            .Append(Escape(string.Join(';', peptide.ProteinGroups)));

    /// <summary>One line per warning, in order.</summary>
    public static void AppendWarnings(StringBuilder text, IReadOnlyList<string> warnings)
    {
        foreach (var w in warnings) text.Append("warning\t").Append(Escape(w)).Append('\n');
    }
}
