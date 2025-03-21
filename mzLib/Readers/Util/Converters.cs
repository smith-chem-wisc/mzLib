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

    /// <summary>
    /// Converts the chemical formula from MsPathFinderT to MetaMorpheus
    /// MsPathFinderT: "C(460) H(740) N(136) O(146) S(0)"
    /// MetaMorpheus: "C460H740N136O146S"
    /// </summary>
    internal class MsPathFinderTCompositionToChemicalFormulaConverter : DefaultTypeConverter
    {
        public override object ConvertFromString(string text, IReaderRow row, MemberMapData memberMapData)
        {
            var composition = text.Split(' ').Where(p => p != "").ToArray();
            var chemicalFormula = new Chemistry.ChemicalFormula();
            foreach (var element in composition)
            {
                var elementSplit = element.Split('(');
                var elementName = elementSplit[0];
                var elementCount = int.Parse(elementSplit[1].Replace(")", ""));
                chemicalFormula.Add(elementName, elementCount);
            }
            return chemicalFormula;
        }

        public override string ConvertToString(object value, IWriterRow row, MemberMapData memberMapData)
        {
            var chemicalFormula = value as Chemistry.ChemicalFormula ?? throw new Exception("Cannot convert input to ChemicalFormula");
            var sb = new StringBuilder();

            bool onNumber = false;
            foreach (var character in chemicalFormula.Formula)
            {
                if (!char.IsDigit(character)) // if is a letter
                {
                    if (onNumber)
                    {
                        sb.Append(") " + character);
                        onNumber = false;
                    }
                    else
                        sb.Append(character);
                }
                else
                {
                    if (!onNumber)
                    {
                        sb.Append("(" + character);
                        onNumber = true;
                    }
                    else
                        sb.Append(character);
                }
            }

            var stringForm = sb.ToString();
            if (char.IsDigit(stringForm.Last()))
                stringForm += ")";
            else
                stringForm += "(1)";

            return stringForm;
        }
    }
}
