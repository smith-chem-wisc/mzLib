using System.Collections.Generic;
using System.ComponentModel;
using Omics.SequenceConversion;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.AbstractClasses;

public abstract class KoinaModelBase<TModelInput, TModelOutput>
{
    private readonly KoinaSequenceConverter _sequenceConverter = new();

    #region Model Metadata
    public abstract string ModelName { get; }
    public abstract int MaxBatchSize { get; }
    public abstract int MaxNumberOfBatchesPerRequest { get; init; }
    public abstract int ThrottlingDelayInMilliseconds { get; init; }
    public abstract int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds { get; }
    #endregion

    #region Input Sequence Validation Constraints
    public abstract SequenceConversionHandlingMode ModHandlingMode { get; init; }
    public abstract int MaxPeptideLength { get; }
    public abstract int MinPeptideLength { get; }

    public virtual IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>();
    protected virtual UnimodSequenceFormatSchema UnimodSchema => UnimodSequenceFormatSchema.Instance;
    protected virtual string AllowedAminoAcidPattern => "^[ACDEFGHIKLMNPQRSTVWY]+$";
    #endregion

    #region Required Client Methods for Koina API Interaction
    protected abstract List<Dictionary<string, object>> ToBatchedRequests(List<TModelInput> validInputs);
    #endregion

    #region Validation and Modification Handling
    /// <summary>
    /// Validates a peptide sequence against model constraints for modifications and basic sequence requirements.
    /// Handles incompatible modifications according to the specified ModHandlingMode.
    /// </summary>
    protected virtual string? TryCleanSequence(
        string sequence,
        out string? apiSequence,
        out WarningException? warning)
    {
        var result = _sequenceConverter.Convert(
            sequence,
            AllowedUnimodIds,
            UnimodSchema,
            AllowedAminoAcidPattern,
            MinPeptideLength,
            MaxPeptideLength,
            ModHandlingMode);

        if (result == null)
        {
            apiSequence = null;
            warning = null;
            return null;
        }

        apiSequence = result.ApiSequence;
        warning = result.Warning;
        return apiSequence;
    }

    protected virtual bool IsValidBaseSequence(string baseSequence)
    {
        return _sequenceConverter.IsValidBaseSequence(
            baseSequence,
            AllowedAminoAcidPattern,
            MinPeptideLength,
            MaxPeptideLength);
    }
    #endregion
}
