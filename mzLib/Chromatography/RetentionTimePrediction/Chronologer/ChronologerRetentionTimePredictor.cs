using Chromatography.RetentionTimePrediction.Util;
using System.Text.RegularExpressions;
using TorchSharp;
using static TorchSharp.torch;

namespace Chromatography.RetentionTimePrediction.Chronologer;

/// <summary>
/// Chronologer-based retention time predictor using deep learning.
/// Predicts C18 retention times reported in % ACN.
/// </summary>
public class ChronologerRetentionTimePredictor : RetentionTimePredictor, IDisposable
{
    private readonly Chronologer _model;
    protected override int MaxSequenceLength => 50;
    private int EncodedLength => MaxSequenceLength + 2; // +2 for N/C termini tokens
    public override string PredictorName => "Chronologer";
    public override SeparationType SeparationType => SeparationType.HPLC;

    /// <summary>
    /// Initializes a new Chronologer predictor with custom weights file.
    /// </summary>
    public ChronologerRetentionTimePredictor(
        IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.RemoveIncompatibleMods,
        string? weightsPath = null)
        : base(modHandlingMode)
    {
        _model = weightsPath != null
            ? new Chronologer(weightsPath)
            : new Chronologer();
    }

    protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        var baseSequence = peptide.BaseSequence;
        if (baseSequence.Any(aa => Array.IndexOf(CanonicalAminoAcids, aa) == -1))
        {
            failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
            return false;
        }

        return base.ValidateBasicConstraints(peptide, out failureReason);
    }

    protected override double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null)
    {
        // Get formatted sequence if not provided
        formattedSequence ??= GetFormattedSequence(peptide, out RetentionTimeFailureReason? failureReason);
        if (formattedSequence == null)
            return null;

        // Encode to tensor
        var ids = new long[EncodedLength]; // Zero-padded
        for (int i = 0; i < formattedSequence.Length; i++)
        {
            if (!CodeToInt.TryGetValue(formattedSequence[i], out int v))
                return null; // Invalid character

            ids[i] = v;
        }

        // Output shape: [1, MaxPepLen+2], dtype int64
        Tensor sequenceTensor = tensor(ids, dtype: ScalarType.Int64).reshape(1, EncodedLength);

        // Predict retention time
        try
        {
            var prediction = _model.Predict(sequenceTensor);
            return prediction[0].ToDouble();
        }
        finally
        {
            sequenceTensor?.Dispose();
        }
    }

    public override string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;

        // Get full sequence with mass shifts
        string workingSequence = peptide.FullSequenceWithMassShifts;

        // Replace mass shifts with chronologer dictionary codes
        foreach (var (pattern, replacement) in ModificationPatterns)
        {
            workingSequence = pattern.Replace(workingSequence, replacement);
        }

        // Add N-terminus token
        workingSequence = AddNTerminusToken(workingSequence);

        // Add C-terminus token
        workingSequence += "_";

        // At this point we have replaced everything that is chronologer compatible with its chronologer dictionary representation. 
        // If we have any more [] in the full sequence, it means there are incompatible modifications.
        if (!workingSequence.Contains('[') && !workingSequence.Contains(']')) 
            return workingSequence;

        switch (ModHandlingMode)
        {
            case IncompatibleModHandlingMode.ReturnNull:
                failureReason = RetentionTimeFailureReason.IncompatibleModifications;
                return null;

            case IncompatibleModHandlingMode.RemoveIncompatibleMods:
                // Remove incompatible modification annotations from full sequence
                var sb = new System.Text.StringBuilder();
                int i = 0;
                while (i < workingSequence.Length)
                {
                    if (workingSequence[i] == '[')
                    {
                        // Found modification
                        int closeIdx = workingSequence.IndexOf(']', i);
                        if (closeIdx != -1)
                        {
                            // Skip modification annotation
                            i = closeIdx + 1;
                            continue;
                        }
                    }
                    // Copy character
                    sb.Append(workingSequence[i]);
                    i++;
                }
                return sb.ToString();

            // Use base sequence without modifications with default termini
            case IncompatibleModHandlingMode.UsePrimarySequence:
                return $"-{peptide.BaseSequence}_";

            case IncompatibleModHandlingMode.ThrowException:
            default:
                throw new IncompatibleModificationException(peptide.FullSequence, workingSequence, PredictorName);
        }
    }

    #region Sequence Encoding 

    // Chronologer alphabet
    // 20 canonical (1..20) + 17 modified (21..37) + 7 N/C states (38..44) + 10 user slots (45..54)
    private static readonly char[] Residues = (
        "ACDEFGHIKLMNPQRSTVWY" +     // 1-20: canonical amino acids
        "cmdestyabunopqrxz" +         // 21-37: modified residues
        "-^()&*_" +                   // 38-44: N/C terminus states
        "0123456789"                  // 45-54: user-defined slots
    ).ToCharArray();

    private static readonly Dictionary<char, int> CodeToInt =
        Residues.Select((c, i) => (c, i + 1)).ToDictionary(t => t.c, t => t.Item2);

    // Compiled regex patterns for performance
    private static readonly (Regex pattern, string replacement)[] ModificationPatterns = new[]
    {
        (new Regex(@"M\[\+15\.99\d*\]", RegexOptions.Compiled), "m"),  // Oxidation on M
        (new Regex(@"C\[\+57\.02\d*\]", RegexOptions.Compiled), "c"),  // Carbamidomethyl on C
        (new Regex(@"C\[\+39\.99\d*\]", RegexOptions.Compiled), "d"),  // Alternative C mod
        (new Regex(@"\[\-18\.01\d*\]E", RegexOptions.Compiled), "e"), // PyroGlu from E (prefix)
        (new Regex(@"E\[\-18\.01\d*\]", RegexOptions.Compiled), "e"),  // PyroGlu from E (suffix)
        (new Regex(@"\[\-17\.02\d*\]Q", RegexOptions.Compiled), "e"), // PyroGlu from Q (prefix)
        (new Regex(@"Q\[\-17\.02\d*\]", RegexOptions.Compiled), "e"),  // PyroGlu from Q (suffix)
        (new Regex(@"S\[\+79\.96\d*\]", RegexOptions.Compiled), "s"),  // Phosphorylation on S
        (new Regex(@"T\[\+79\.96\d*\]", RegexOptions.Compiled), "t"),  // Phosphorylation on T
        (new Regex(@"Y\[\+79\.96\d*\]", RegexOptions.Compiled), "y"),  // Phosphorylation on Y
        (new Regex(@"K\[\+42\.01\d*\]", RegexOptions.Compiled), "a"),  // Acetylation on K
        (new Regex(@"K\[\+100\.0\d*\]", RegexOptions.Compiled), "b"),  // Succinylation on K
        (new Regex(@"K\[\+114\.0\d*\]", RegexOptions.Compiled), "u"),  // Ubiquitination on K
        (new Regex(@"K\[\+14\.01\d*\]", RegexOptions.Compiled), "n"),  // Methylation on K
        (new Regex(@"K\[\+28\.03\d*\]", RegexOptions.Compiled), "o"),  // Dimethylation on K
        (new Regex(@"K\[\+42\.04\d*\]", RegexOptions.Compiled), "p"),  // Trimethylation on K
        (new Regex(@"R\[\+14\.01\d*\]", RegexOptions.Compiled), "q"),  // Methylation on R
        (new Regex(@"R\[\+28\.03\d*\]", RegexOptions.Compiled), "r"),  // Dimethylation on R
        (new Regex(@"K\[\+224\.1\d*\]", RegexOptions.Compiled), "z"),  // GlyGly on K
        (new Regex(@"K\[\+229\.1\d*\]", RegexOptions.Compiled), "x"),  // Heavy GlyGly on K
    };

    // N-terminus modification codes
    private static readonly Dictionary<string, char> NTerminusCodes = new()
    {
        { "+42.01", '^' },  // N-term acetylation
        { "+224.1", '&' },  // N-term GlyGly
        { "+229.1", '*' }   // N-term heavy GlyGly
    };

    /// <summary>
    /// Adds appropriate N-terminus token based on modifications or default state.
    /// </summary>
    private static string AddNTerminusToken(string seq)
    {
        // Check for PyroGlu at N-terminus
        if (seq[0] == 'd') // pyroGlu at first position
            return ")" + seq;

        if (seq[0] == 'e') // cyclized CAM-Cys at first
            return "(" + seq;

        // Check for N-terminal mass modification
        if (seq[0] == '[')
        {
            // grab [+xx.xx]
            int close = seq.IndexOf(']');

            string key = seq.Substring(1, 6);
            if (!NTerminusCodes.TryGetValue(key, out char nterm))
                nterm = '-'; // Unknown modification - use default

            return nterm + seq.Substring(close + 1);
        }

        // Free N-terminus
        return "-" + seq;
    }

    #endregion

    public void Dispose()
    {
        _model.Dispose();
    }
}