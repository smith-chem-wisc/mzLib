using System.Reflection;
using CsvHelper.Configuration.Attributes;
using MzLibUtil;
using Parquet;
using Parquet.Schema;

namespace Readers;

/// <summary>
/// DIA-NN's main report as written by DIA-NN 2.x: a report.parquet with one row per precursor per run.
/// It yields the same <see cref="DiaNnPrecursor"/> records as the 1.8/1.9 TSV read by
/// <see cref="DiaNnReportFile"/>, so FlashLFQ and every other consumer of that type read it unchanged.
/// <para>
/// DIA-NN 2.x changed the columns it writes:
/// <list type="bullet">
/// <item>Removed: File.Name, PG.Quantity, PG.Normalised, Genes.Quantity, Genes.Normalised and MS2.Scan.
/// Their properties are left "not reported": null for strings and nullable numbers, NaN for
/// other doubles, and 0 for <see cref="DiaNnPrecursor.Ms2ScanNumber"/>, the one integer.</item>
/// <item>Renamed: Lib.Index is now Precursor.Lib.Index.</item>
/// <item>Added: a Decoy column, written when DIA-NN runs with --report-decoys. Decoy rows are
/// skipped, because <see cref="DiaNnPrecursor.IsDecoy"/> is always false and a loaded decoy would be
/// quantified as a target.</item>
/// </list>
/// </para>
/// </summary>
public class DiaNnParquetReportFile : DiaNnReportFile
{
    /// <summary>
    /// The columns that identify a DIA-NN main report. Run is required in place of the TSV's
    /// File.Name, which 2.x no longer writes.
    /// </summary>
    private static readonly string[] IdentifyingColumns = ["Precursor.Id", "Stripped.Sequence", "Run"];

    /// <summary>
    /// 1.8/1.9 column name to the name DIA-NN 2.x writes instead.
    /// </summary>
    private static readonly Dictionary<string, string> RenamedInDiaNn2 = new()
    {
        ["Lib.Index"] = "Precursor.Lib.Index",
    };

    private const string DecoyColumn = "Decoy";

    /// <summary>
    /// Every <see cref="DiaNnPrecursor"/> property that maps to a report column, with the name that
    /// column has in a DIA-NN 2.x report. The CsvHelper [Name] attributes stay the single source of
    /// column names for both readers.
    /// </summary>
    private static readonly (PropertyInfo Property, string Column)[] MappedProperties = typeof(DiaNnPrecursor)
        .GetProperties(BindingFlags.Public | BindingFlags.Instance)
        .Where(property => property.CanWrite)
        .Select(property => (Property: property, Name: property.GetCustomAttribute<NameAttribute>()?.Names.FirstOrDefault()))
        .Where(mapping => mapping.Name is not null)
        .Select(mapping => (mapping.Property, RenamedInDiaNn2.GetValueOrDefault(mapping.Name!, mapping.Name!)))
        .ToArray();

    public override SupportedFileType FileType => SupportedFileType.DiaNnReportParquet;

    public DiaNnParquetReportFile(string filePath) : base(filePath) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public DiaNnParquetReportFile() : base()
    {
        Software = Software.DiaNn;
    }

    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    /// <exception cref="MzLibException">The file is not a readable parquet, naming the file.</exception>
    public override void LoadResults()
    {
        if (!File.Exists(FilePath))
            throw new FileNotFoundException($"DIA-NN parquet report not found: '{FilePath}'", FilePath);

        try
        {
            Results = ReadTargets(FilePath);
        }
        catch (Exception e) when (e is not MzLibException)
        {
            throw new MzLibException($"Could not read DIA-NN parquet report '{FilePath}': {e.Message}", e);
        }
    }

    /// <exception cref="NotSupportedException">Always.</exception>
    public override void WriteResults(string outputPath) =>
        throw new NotSupportedException("Writing a DIA-NN parquet report is not supported; DiaNnParquetReportFile only reads it.");

    /// <summary>
    /// True when <paramref name="filePath"/> is a parquet file carrying the columns of a DIA-NN main
    /// report. A file that exists but is not a readable parquet gives false rather than throwing.
    /// </summary>
    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    /// <exception cref="FileNotFoundException">The file does not exist.</exception>
    internal static bool HasDiaNnReportColumns(string filePath)
    {
        if (!File.Exists(filePath))
            throw new FileNotFoundException($"File not found: '{filePath}'", filePath);

        try
        {
            ParquetReader reader = ParquetReader.CreateAsync(filePath).GetAwaiter().GetResult();
            try
            {
                var columns = reader.Schema.GetDataFields().Select(field => field.Name).ToHashSet();
                return IdentifyingColumns.All(columns.Contains);
            }
            finally
            {
                reader.DisposeAsync().AsTask().GetAwaiter().GetResult();
            }
        }
        catch (Exception e) when (e is not OutOfMemoryException)
        {
            // Parquet.Net reports a malformed file through several exception types; to a caller asking
            // "is this a DIA-NN report?" they all mean no.
            return false;
        }
    }

    private static List<DiaNnPrecursor> ReadTargets(string filePath)
    {
        ParquetReader reader = ParquetReader.CreateAsync(filePath).GetAwaiter().GetResult();
        try
        {
            return ReadTargets(reader, filePath);
        }
        finally
        {
            reader.DisposeAsync().AsTask().GetAwaiter().GetResult();
        }
    }

    private static List<DiaNnPrecursor> ReadTargets(ParquetReader reader, string filePath)
    {
        Dictionary<string, DataField> fields = reader.Schema.GetDataFields().ToDictionary(field => field.Name);

        var missing = IdentifyingColumns.Where(column => !fields.ContainsKey(column)).ToList();
        if (missing.Count > 0)
            throw new MzLibException($"'{filePath}' is not a DIA-NN report: missing column(s) {string.Join(", ", missing)}");

        // Strings repeat heavily across rows (run names, proteins, one peptide per run), so share
        // one instance of each, as the TSV reader does.
        var pool = new Dictionary<string, string>(StringComparer.Ordinal);
        var results = new List<DiaNnPrecursor>();

        for (int rowGroup = 0; rowGroup < reader.RowGroupCount; rowGroup++)
        {
            using ParquetRowGroupReader rowGroupReader = reader.OpenRowGroupReader(rowGroup);
            int rowCount = checked((int)rowGroupReader.RowCount);

            bool[] isDecoy = fields.TryGetValue(DecoyColumn, out DataField? decoyField)
                ? ReadColumn(rowGroupReader, decoyField, rowCount).Select(value => value is not null && System.Convert.ToInt64(value) != 0).ToArray()
                : new bool[rowCount];

            var precursors = new DiaNnPrecursor[rowCount];
            for (int row = 0; row < rowCount; row++)
                precursors[row] = NewWithUnreportedDefaults();

            foreach (var (property, column) in MappedProperties)
            {
                if (!fields.TryGetValue(column, out DataField? field))
                    continue;

                object?[] values = ReadColumn(rowGroupReader, field, rowCount);
                for (int row = 0; row < rowCount; row++)
                    if (!isDecoy[row])
                        property.SetValue(precursors[row], ConvertValue(values[row], property.PropertyType, pool));
            }

            for (int row = 0; row < rowCount; row++)
                if (!isDecoy[row])
                    results.Add(precursors[row]);
        }

        return results;
    }

    private static readonly MethodInfo ReadStructColumnMethod =
        typeof(DiaNnParquetReportFile).GetMethod(nameof(ReadStructColumn), BindingFlags.NonPublic | BindingFlags.Static)!;

    /// <summary>
    /// Reads one column of a row group as boxed values, null where the file holds a null. Parquet.Net
    /// reads are typed by the column's CLR type, so dispatch on it.
    /// </summary>
    private static object?[] ReadColumn(ParquetRowGroupReader rowGroupReader, DataField field, int rowCount)
    {
        Type clrType = Nullable.GetUnderlyingType(field.ClrType) ?? field.ClrType;
        // Parquet.Net 6 describes a string column as ReadOnlyMemory<char>; its string[] overload decodes it.
        if (clrType == typeof(string) || clrType == typeof(ReadOnlyMemory<char>))
        {
            var strings = new string?[rowCount];
            rowGroupReader.ReadAsync(field, strings).AsTask().GetAwaiter().GetResult();
            return strings;
        }

        return (object?[])ReadStructColumnMethod.MakeGenericMethod(clrType).Invoke(null, [rowGroupReader, field, rowCount])!;
    }

    private static object?[] ReadStructColumn<T>(ParquetRowGroupReader rowGroupReader, DataField field, int rowCount) where T : struct
    {
        if (field.IsNullable)
        {
            var nullable = new T?[rowCount];
            rowGroupReader.ReadAsync(field, nullable.AsMemory()).AsTask().GetAwaiter().GetResult();
            return nullable.Select(value => value.HasValue ? (object?)value.Value : null).ToArray();
        }

        var values = new T[rowCount];
        rowGroupReader.ReadAsync(field, values.AsMemory()).AsTask().GetAwaiter().GetResult();
        return values.Select(value => (object?)value).ToArray();
    }

    /// <summary>
    /// A record whose non-nullable doubles start as NaN, so a column DIA-NN 2.x does not write reads
    /// as "not reported" rather than as zero. Reference types and nullable numbers already default
    /// to null.
    /// </summary>
    private static DiaNnPrecursor NewWithUnreportedDefaults()
    {
        var precursor = new DiaNnPrecursor();
        foreach (var (property, _) in MappedProperties)
            if (property.PropertyType == typeof(double))
                property.SetValue(precursor, double.NaN);
        return precursor;
    }

    private static object? ConvertValue(object? value, Type target, Dictionary<string, string> pool)
    {
        if (value is null)
            return target == typeof(double) ? double.NaN : null;

        if (target == typeof(string))
        {
            string text = (string)value;
            if (pool.TryGetValue(text, out string? shared))
                return shared;
            pool[text] = text;
            return text;
        }

        Type underlying = Nullable.GetUnderlyingType(target) ?? target;
        if (underlying == typeof(bool))
            return System.Convert.ToInt64(value) != 0;
        return System.Convert.ChangeType(value, underlying, System.Globalization.CultureInfo.InvariantCulture);
    }
}
