namespace Omics.BioPolymerGroup;

/// <summary>
/// One column of tab-separated output: its header text paired with an accessor that reads the
/// value out of an item.
/// </summary>
/// <typeparam name="T">The record type being written, e.g. a biopolymer group.</typeparam>
public sealed class TsvColumn<T>
{
    /// <summary>Column name written to the header line.</summary>
    public string Header { get; }

    /// <summary>Reads this column's value from one item.</summary>
    public Func<T, string> GetValue { get; }

    public TsvColumn(string header, Func<T, string> getValue)
    {
        Header = header ?? throw new ArgumentNullException(nameof(header));
        GetValue = getValue ?? throw new ArgumentNullException(nameof(getValue));
    }
}

/// <summary>
/// Writes records as a tab-separated file from a schema — an ordered list of <see cref="TsvColumn{T}"/>.
///
/// The schema is built once for a dataset and then applied to every record, so the header and all
/// rows are necessarily the same width. Records do not format themselves: a record that lacks a
/// value for some column contributes an empty field rather than omitting the column, which is what
/// keeps a ragged dataset from shifting every subsequent field on the affected rows.
/// </summary>
public static class TsvWriter
{
    /// <summary>
    /// Maximum length for a single field. Longer values are truncated by <see cref="Truncate"/>.
    ///
    /// Default is 32,000 characters, slightly below Excel's cell limit of 32,767, so output files
    /// open in Excel without truncation or corruption. Set to 0 or negative to disable truncation
    /// (useful for programmatic processing where Excel compatibility does not matter).
    /// </summary>
    /// <remarks>
    /// Excel specification: a cell can contain up to 32,767 characters.
    /// See: https://support.microsoft.com/en-us/office/excel-specifications-and-limits
    /// </remarks>
    public static int MaxFieldLength { get; set; } = 32000;

    /// <summary>
    /// Caps a field at <see cref="MaxFieldLength"/>. Null and empty both render as an empty field.
    /// </summary>
    public static string Truncate(string? value)
    {
        if (string.IsNullOrEmpty(value))
            return string.Empty;

        if (MaxFieldLength <= 0 || value.Length <= MaxFieldLength)
            return value;

        return value.Substring(0, MaxFieldLength);
    }

    /// <summary>Renders the header line for a schema.</summary>
    public static string HeaderLine<T>(IReadOnlyList<TsvColumn<T>> schema)
        => string.Join('\t', schema.Select(c => c.Header));

    /// <summary>
    /// Renders one record. A column accessor returning null yields an empty field so the row keeps
    /// its alignment with the header.
    /// </summary>
    public static string RowLine<T>(IReadOnlyList<TsvColumn<T>> schema, T item)
        => string.Join('\t', schema.Select(c => c.GetValue(item) ?? string.Empty));

    /// <summary>
    /// Writes the header followed by one line per item.
    /// </summary>
    public static void Write<T>(TextWriter output, IReadOnlyList<TsvColumn<T>> schema, IEnumerable<T> items)
    {
        output.WriteLine(HeaderLine(schema));
        foreach (var item in items)
        {
            output.WriteLine(RowLine(schema, item));
        }
    }
}
