using CsvHelper;
using CsvHelper.Configuration;
using CsvHelper.TypeConversion;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Converts a list of doubles delimited by semicolons to a list of doubles
    /// To be used with CsvHelper
    /// </summary>
    public class SemicolonDelimitedToDoubleListConverter : DefaultTypeConverter
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

    public class DashToNullOrDoubleConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            return text == "-" ? null : double.Parse(text);
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            return value as double? == null ? "-" : value.ToString();
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

    public class CommaDelimitedToStringArrayTypeConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            return text.Split(',').Where(p => p != "").ToArray();
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var list = value as IEnumerable<string> ?? throw new MzLibException("Cannot convert input to IEnumerable<string>");
            return string.Join(',', list);
        }
    }
}
