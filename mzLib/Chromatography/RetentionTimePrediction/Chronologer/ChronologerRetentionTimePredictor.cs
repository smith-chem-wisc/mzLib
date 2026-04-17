using TorchSharp;
using static TorchSharp.torch;
using Omics.SequenceConversion;

namespace Chromatography.RetentionTimePrediction.Chronologer;

/// <summary>
/// Chronologer-based retention time predictor using deep learning.
/// Predicts C18 retention times reported in % ACN.
/// </summary>
public class ChronologerRetentionTimePredictor : RetentionTimePredictor
{
    private static readonly SequenceConversionService ConversionService = SequenceConversionService.Default;
    private static readonly string ChronologerFormatName = ChronologerSequenceFormatSchema.Instance.FormatName;
    private static readonly string MzLibFormatName = MzLibSequenceFormatSchema.Instance.FormatName;
    private static readonly string MassShiftFormatName = MassShiftSequenceFormatSchema.Instance.FormatName;

    private readonly Chronologer _model;
    private readonly object _modelLock = new();
    private bool _disposed;

    protected override int MaxSequenceLength => ChronologerSequenceFormatSchema.MaxSequenceLength;
    public override string PredictorName => "Chronologer";
    public override SeparationType SeparationType => SeparationType.HPLC;

    /// <summary>
    /// Initializes a new Chronologer predictor with custom weights file. Uses default pretrained weights if none provided.
    /// </summary>
    public ChronologerRetentionTimePredictor(
        SequenceConversionHandlingMode sequenceHandlingMode = SequenceConversionHandlingMode.RemoveIncompatibleElements,
        string? weightsPath = null)
        : base(sequenceHandlingMode)
    {
        _model = weightsPath != null
            ? new Chronologer(weightsPath)
            : new Chronologer();
    }

    protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        var baseSequence = peptide.BaseSequence;
        if (baseSequence.Any(aa => !CanonicalAminoAcids.Contains(aa)))
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
        var ids = new long[ChronologerSequenceFormatSchema.EncodedLength]; // Zero-padded
        for (int i = 0; i < formattedSequence.Length; i++)
        {
            if (!CodeToInt.TryGetValue(formattedSequence[i], out int v))
                return null; // Invalid character

            ids[i] = v;
        }

        // Output shape: [1, MaxPepLen+2], dtype int64
        using Tensor sequenceTensor = tensor(ids, dtype: ScalarType.Int64).reshape(1, ChronologerSequenceFormatSchema.EncodedLength);


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

        var candidates = new (string? Sequence, string? ForcedFormat)[]
        {
            (peptide.FullSequence, null),
            (peptide.FullSequenceWithMassShifts, MassShiftFormatName),
            (peptide.BaseSequence, MzLibFormatName)
        };

        try
        {
            foreach (var (sequence, forcedFormat) in candidates)
            {
                if (string.IsNullOrWhiteSpace(sequence))
                    continue;

                var attemptWarnings = new ConversionWarnings();

                var sourceFormat = forcedFormat
                    ?? ConversionService.DetectFormat(sequence!)
                    ?? MzLibFormatName;

                var formatted = ConversionService.Convert(
                    sequence!,
                    sourceFormat,
                    ChronologerFormatName,
                    attemptWarnings,
                    SequenceHandlingMode);

                if (formatted != null)
                {
                    return formatted;
                }

                var mappedReason = MapFailureReason(attemptWarnings.FailureReason)
                    ?? (SequenceHandlingMode == SequenceConversionHandlingMode.ReturnNull
                        ? RetentionTimeFailureReason.IncompatibleModifications
                        : null);

                if (SequenceHandlingMode == SequenceConversionHandlingMode.UsePrimarySequence)
                {
                    failureReason = mappedReason;
                    return FormatPrimarySequence(peptide.BaseSequence);
                }

                if (SequenceHandlingMode == SequenceConversionHandlingMode.ReturnNull)
                {
                    failureReason = mappedReason ?? RetentionTimeFailureReason.IncompatibleModifications;
                    return null;
                }

                failureReason ??= mappedReason;
            }
        }
        catch (SequenceConversionException ex)
        {
            HandleConversionException(peptide, ex, ref failureReason);

            if (SequenceHandlingMode == SequenceConversionHandlingMode.UsePrimarySequence)
                return FormatPrimarySequence(peptide.BaseSequence);

            return null;
        }

        failureReason ??= SequenceHandlingMode == SequenceConversionHandlingMode.ReturnNull
            ? RetentionTimeFailureReason.IncompatibleModifications
            : RetentionTimeFailureReason.PredictionError;

        if (SequenceHandlingMode == SequenceConversionHandlingMode.UsePrimarySequence)
        {
            return FormatPrimarySequence(peptide.BaseSequence);
        }

        return null;
    }

    private void HandleConversionException(IRetentionPredictable peptide, SequenceConversionException exception, ref RetentionTimeFailureReason? failureReason)
    {
        failureReason = MapFailureReason(exception.FailureReason);

        if (SequenceHandlingMode == SequenceConversionHandlingMode.ThrowException &&
            exception.FailureReason == ConversionFailureReason.IncompatibleModifications)
        {
            var workingSequence = peptide.FullSequenceWithMassShifts ?? peptide.BaseSequence;
            throw new IncompatibleModificationException(
                peptide.FullSequence,
                workingSequence,
                PredictorName);
        }
    }

    private static string FormatPrimarySequence(string baseSequence)
    {
        return string.Concat(ChronologerSequenceFormatSchema.FreeNTerminus, baseSequence, ChronologerSequenceFormatSchema.CTerminus);
    }

    private static RetentionTimeFailureReason? MapFailureReason(ConversionFailureReason? reason)
    {
        if (!reason.HasValue)
            return null;

        return reason.Value switch
        {
            ConversionFailureReason.IncompatibleModifications => RetentionTimeFailureReason.IncompatibleModifications,
            ConversionFailureReason.InvalidSequence => RetentionTimeFailureReason.EmptySequence,
            _ => RetentionTimeFailureReason.PredictionError
        };
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

    /// <summary>
    /// Releases the underlying Chronologer TorchSharp model. Overrides the base
    /// <see cref="RetentionTimePredictor.Dispose"/> so that virtual dispatch reaches
    /// this implementation even when the caller holds a base-class or
    /// <see cref="IDisposable"/> reference. Safe to call multiple times.
    /// </summary>
    public override void Dispose()
    {
        if (_disposed)
            return;

        _disposed = true;
        _model?.Dispose();
    }
}