using System.Collections.Generic;
using System.ComponentModel;
using Omics.SequenceConversion;
using PredictionClients.Koina.Util;

namespace PredictionClients.Koina.AbstractClasses;

public abstract class KoinaModelBase<TModelInput, TModelOutput>
{
    private readonly KoinaSequenceConverter _sequenceConverter = new();

    #region Model Metadata
    /// <summary>
    /// Gets the model name as registered in the Koina API.
    /// </summary>
    public abstract string ModelName { get; }

    /// <summary>
    /// Gets the maximum number of sequences allowed per API request batch.
    /// Can dig in the Koina github repo to find these values if needed.
    /// </summary>
    public abstract int MaxBatchSize { get; }

    /// <summary>
    /// Gets the maximum number of batches that can be combined into a single API request.
    /// Used to optimize request throughput while respecting API limitations.
    /// </summary>
    public abstract int MaxNumberOfBatchesPerRequest { get; init; }

    /// <summary>
    /// Gets the delay in milliseconds to wait between consecutive API requests.
    /// Used for rate limiting to prevent overwhelming the Koina API server.
    /// </summary>
    public abstract int ThrottlingDelayInMilliseconds { get; init; }       
    /// <summary>
    /// Gets the benchmarked processing time in milliseconds for one batch at MaxBatchSize.
    /// Used for estimating total request duration and optimizing parallelization strategies.
    /// </summary>
    public abstract int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds { get; }
    #endregion

    #region Input Sequence Validation Constraints
    public abstract SequenceConversionHandlingMode ModHandlingMode { get; init; }
    /// <summary>
    /// Gets the maximum allowed peptide base sequence length.
    /// </summary>
    public abstract int MaxPeptideLength { get; }

    /// <summary>
    /// Gets the minimum allowed peptide base sequence length.
    /// </summary>
    public abstract int MinPeptideLength { get; }

    /// <summary>
    /// Unimod Ids that are allowed to b passed to the model. 
    /// </summary>
    public virtual IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>();
    protected virtual UnimodSequenceFormatSchema UnimodSchema => UnimodSequenceFormatSchema.Instance;

    /// <summary>
    /// Gets the regex pattern for validating amino acid sequences.
    /// </summary>
    protected virtual string AllowedAminoAcidPattern => "^[ACDEFGHIKLMNPQRSTVWY]+$";
    #endregion

    #region Required Client Methods for Koina API Interaction
    /// <summary>
    /// Converts peptide sequences and associated data into batched request payloads for the Koina API.
    /// Implementations should group input sequences into batches respecting the MaxBatchSize constraint
    /// and format them according to the specific model's input requirements.
    /// </summary>
    /// <returns>List of request dictionaries, each containing a batch of sequences and parameters</returns>
    /// <remarks>
    /// Each dictionary in the returned list represents one API request batch and should contain:
    /// - Peptide sequences (formatted according to model requirements)
    /// - Model-specific parameters (e.g., charge states, collision energies, NCE values)
    /// - Any additional metadata required by the specific Koina model
    /// Must ensure that only the validated sequences that meet the model's constraints are included in the batches. 
    /// The total number of sequences across all batches should equal PeptideSequences.Count.
    /// </remarks>
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

    #endregion
}
