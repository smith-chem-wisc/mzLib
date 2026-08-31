using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.TypeConversion;
using MzLibUtil;
using System.Text;

namespace Readers
{
    /// <summary>
    /// Converts a list of doubles delimited by semicolons to a list of doubles
    /// To be used with CsvHelper
    /// </summary>
    internal class SemicolonDelimitedToDoubleListConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            var splits = text.Split(';');
            return splits.Select(p => double.Parse(p)).ToList();
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var list = value as IEnumerable<double> ?? throw new MzLibException("Cannot convert input to IEnumerable<double>");
            return string.Join(';', list);
        }
    }

    internal class DashToNullOrDoubleConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            if (text == "-")
                return null;
            return double.TryParse(text, out var result) ? result : 0.0;
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            return value as double? == null ? "-" : value.ToString();
        }
    }

    public class DashToNullOrIntegerConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            if (text == "-")
                return null;
            return int.TryParse(text, out var result) ? result : 0;
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            return value as int? == null ? "-" : value.ToString();
        }
    }

    public class CommaDelimitedToIntegerArrayTypeConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            var splits = text.Split(',');
            var toReturn = splits.Where(p => p != "");
            return toReturn.Select(int.Parse).ToArray();
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var list = value as IEnumerable<int> ?? throw new MzLibException("Cannot convert input to IEnumerable<double>");
            return string.Join(',', list);
        }
    }

    public class CommaDelimitedToDoubleArrayTypeConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            var splits = text.Split(',');
            var toReturn = splits.Where(p => p != "");
            return toReturn.Select(double.Parse).ToArray();
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var list = value as IEnumerable<double> ?? throw new MzLibException("Cannot convert input to IEnumerable<double>");
            return string.Join(',', list);
        }
    }

    /// <summary>
    /// Reads a boolean written as 0 or 1 and writes it back the same way. CsvHelper's default
    /// converter reads 0/1 but writes "False"/"True", which would stop the file being readable by
    /// the tools that produced it.
    /// </summary>
    internal class IntegerBooleanConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            if (string.IsNullOrWhiteSpace(text))
                return false;
            if (int.TryParse(text, out var number))
                return number != 0;
            return bool.TryParse(text, out var parsed) && parsed;
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            return value is bool flag && flag ? "1" : "0";
        }
    }

    internal class CommaDelimitedToStringArrayTypeConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            return text.Split(',')
                .Where(p => p != "")
                .Select(p => p.Trim())
                .ToArray();
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var list = value as IEnumerable<string> ?? throw new MzLibException("Cannot convert input to IEnumerable<string>");
            return string.Join(',', list);
        }
    }
}
