using System;
using System.Text;
using System.Text.RegularExpressions;

namespace Omics.Modifications.Conversion;

public static class ChronologerSequenceFormatter
{
    private static readonly (Regex pattern, string replacement)[] ModificationPatterns = new[]
    {
        (new Regex(@"M\[\+15\.99\d*\]", RegexOptions.Compiled), "m"),
        (new Regex(@"C\[\+57\.02\d*\]", RegexOptions.Compiled), "c"),
        (new Regex(@"C\[\+39\.99\d*\]", RegexOptions.Compiled), "d"),
        (new Regex(@"\[\-18\.01\d*\]E", RegexOptions.Compiled), "e"),
        (new Regex(@"E\[\-18\.01\d*\]", RegexOptions.Compiled), "e"),
        (new Regex(@"\[\-17\.02\d*\]Q", RegexOptions.Compiled), "e"),
        (new Regex(@"Q\[\-17\.02\d*\]", RegexOptions.Compiled), "e"),
        (new Regex(@"S\[\+79\.96\d*\]", RegexOptions.Compiled), "s"),
        (new Regex(@"T\[\+79\.96\d*\]", RegexOptions.Compiled), "t"),
        (new Regex(@"Y\[\+79\.96\d*\]", RegexOptions.Compiled), "y"),
        (new Regex(@"K\[\+42\.01\d*\]", RegexOptions.Compiled), "a"),
        (new Regex(@"K\[\+100\.0\d*\]", RegexOptions.Compiled), "b"),
        (new Regex(@"K\[\+114\.0\d*\]", RegexOptions.Compiled), "u"),
        (new Regex(@"K\[\+14\.01\d*\]", RegexOptions.Compiled), "n"),
        (new Regex(@"K\[\+28\.03\d*\]", RegexOptions.Compiled), "o"),
        (new Regex(@"K\[\+42\.04\d*\]", RegexOptions.Compiled), "p"),
        (new Regex(@"R\[\+14\.01\d*\]", RegexOptions.Compiled), "q"),
        (new Regex(@"R\[\+28\.03\d*\]", RegexOptions.Compiled), "r"),
        (new Regex(@"K\[\+224\.1\d*\]", RegexOptions.Compiled), "z"),
        (new Regex(@"K\[\+229\.1\d*\]", RegexOptions.Compiled), "x"),
    };

    private static readonly Dictionary<string, char> NTerminusCodes = new()
    {
        { "+42.01", '^' },
        { "+224.1", '&' },
        { "+229.1", '*' }
    };

    public static bool TryFormatChronologerSequence(
        IBioPolymerWithSetMods bioPolymer,
        string massShiftSequence,
        SequenceConversionHandlingMode handlingMode,
        out string? formattedSequence,
        out SequenceConversionFailureReason? reason)
    {
        ArgumentNullException.ThrowIfNull(bioPolymer);
        return TryFormatChronologerSequence(
            bioPolymer.BaseSequence,
            massShiftSequence,
            handlingMode,
            out formattedSequence,
            out reason,
            bioPolymer.FullSequence);
    }

    public static bool TryFormatChronologerSequence(
        string baseSequence,
        string massShiftSequence,
        SequenceConversionHandlingMode handlingMode,
        out string? formattedSequence,
        out SequenceConversionFailureReason? reason,
        string? originalSequence = null)
    {
        formattedSequence = null;
        reason = null;

        var workingSequence = ApplyModificationPatterns(massShiftSequence);
        workingSequence = AddNTerminusToken(workingSequence);
        workingSequence += "_";

        if (!workingSequence.Contains('[') && !workingSequence.Contains(']'))
        {
            formattedSequence = workingSequence;
            return true;
        }

        switch (handlingMode)
        {
            case SequenceConversionHandlingMode.RemoveIncompatibleMods:
                formattedSequence = StripAnnotations(workingSequence);
                reason = SequenceConversionFailureReason.ModificationsRemoved;
                return true;

            case SequenceConversionHandlingMode.UsePrimarySequence:
                formattedSequence = $"-{baseSequence}_";
                reason = SequenceConversionFailureReason.UsedPrimarySequence;
                return true;

            case SequenceConversionHandlingMode.KeepOriginalAnnotation:
                formattedSequence = workingSequence;
                reason = SequenceConversionFailureReason.IncompatibleModifications;
                return true;

            case SequenceConversionHandlingMode.ReturnNull:
                formattedSequence = null;
                reason = SequenceConversionFailureReason.ReturnedNull;
                return false;

            case SequenceConversionHandlingMode.UseMassShifts:
                formattedSequence = null;
                reason = SequenceConversionFailureReason.IncompatibleModifications;
                return false;

            case SequenceConversionHandlingMode.ThrowException:
            default:
                throw new IncompatibleModificationException(originalSequence ?? baseSequence, workingSequence, "Chronologer");
        }
    }

    private static string ApplyModificationPatterns(string sequence)
    {
        var working = sequence;
        foreach (var (pattern, replacement) in ModificationPatterns)
        {
            working = pattern.Replace(working, replacement);
        }

        return working;
    }

    private static string AddNTerminusToken(string sequence)
    {
        if (sequence.Length == 0)
        {
            return "-";
        }

        if (sequence[0] == 'd')
        {
            return ")" + sequence;
        }

        if (sequence[0] == 'e')
        {
            return "(" + sequence;
        }

        if (sequence[0] == '[')
        {
            var close = sequence.IndexOf(']');
            if (close > 0)
            {
                var length = Math.Min(6, Math.Max(0, close - 1));
                var key = length > 0 ? sequence.Substring(1, length) : string.Empty;
                if (!NTerminusCodes.TryGetValue(key, out var nterm))
                {
                    nterm = '-';
                }

                return nterm + sequence[(close + 1)..];
            }
        }

        return "-" + sequence;
    }

    private static string StripAnnotations(string sequence)
    {
        var builder = new StringBuilder(sequence.Length);
        var i = 0;
        while (i < sequence.Length)
        {
            if (sequence[i] == '[')
            {
                var closeIdx = sequence.IndexOf(']', i);
                if (closeIdx == -1)
                {
                    break;
                }

                i = closeIdx + 1;
                continue;
            }

            builder.Append(sequence[i]);
            i++;
        }

        return builder.ToString();
    }
}
