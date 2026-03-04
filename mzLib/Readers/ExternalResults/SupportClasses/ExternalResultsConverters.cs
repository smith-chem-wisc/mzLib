using CsvHelper.Configuration;
using CsvHelper.TypeConversion;
using CsvHelper;
using System.Text;
using MzLibUtil;

namespace Readers;

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
        var chemicalFormula = value as Chemistry.ChemicalFormula ?? throw new MzLibException("Cannot convert input to ChemicalFormula");
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

