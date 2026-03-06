using TorchSharp;
using static TorchSharp.torch;
using Omics;
using Omics.Modifications.Conversion;

namespace Chromatography.RetentionTimePrediction.Chronologer;

/// <summary>
/// Chronologer-based retention time predictor using deep learning.
/// Predicts C18 retention times reported in % ACN.
/// </summary>
public class ChronologerRetentionTimePredictor : RetentionTimePredictor, IDisposable
{
    private readonly Chronologer _model;
    private readonly object _modelLock = new(); 
    private bool _disposed;

    protected override int MaxSequenceLength => 50;
    private int EncodedLength => MaxSequenceLength + 2; // +2 for N/C termini tokens
    public override string PredictorName => "Chronologer";
    public override SeparationType SeparationType => SeparationType.HPLC;

    /// <summary>
    /// Initializes a new Chronologer predictor with custom weights file. Uses default pretrained weights if none provided.
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
        if (_disposed)
            throw new ObjectDisposedException(nameof(ChronologerRetentionTimePredictor));

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
        using Tensor sequenceTensor = tensor(ids, dtype: ScalarType.Int64).reshape(1, EncodedLength);


        // Predict retention time - keep both prediction AND disposal inside lock
        lock (_modelLock)
        {
            using Tensor prediction = _model.Predict(sequenceTensor);
            return prediction[0].ToDouble();
        }
    }

    public override string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;
        var handlingMode = ModHandlingMode.ToSequenceConversionHandlingMode();

        if (peptide is IBioPolymerWithSetMods withSetMods)
        {
            if (SequenceTargetConverter.TryConvert(withSetMods, SequenceConversionTarget.Chronologer, handlingMode, out var converted, out _))
            {
                return converted;
            }

            failureReason = RetentionTimeFailureReason.IncompatibleModifications;
            return null;
        }

        var massShiftSequence = peptide.FullSequenceWithMassShifts;
        if (ChronologerSequenceFormatter.TryFormatChronologerSequence(peptide.BaseSequence, massShiftSequence, handlingMode, out var formatted, out _))
        {
            return formatted;
        }

        failureReason = RetentionTimeFailureReason.IncompatibleModifications;
        return null;
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

    #endregion

    public void Dispose()
    {
        if (_disposed)
            return;

        _disposed = true;
        _model?.Dispose();
    }
}
