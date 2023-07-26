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
}
