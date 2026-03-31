using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace PredictionClients.Koina.AbstractClasses;

public abstract class KoinaModelBase<TModelInput, TModelOutput>
{
    private static readonly Regex BaseStripper = new(@"\[[^\]]+\]", RegexOptions.Compiled);

    protected KoinaModelBase(ISequenceConverter sequenceConverter)
    {
        SequenceConverter = sequenceConverter ?? throw new ArgumentNullException(nameof(sequenceConverter));
    }

    protected ISequenceConverter SequenceConverter { get; }

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
        apiSequence = null;
        warning = null;

        var rawBase = BaseStripper.Replace(sequence, string.Empty);
        if (!Regex.IsMatch(rawBase, AllowedAminoAcidPattern))
        {
            HandleFailure(ModHandlingMode, "Invalid base sequence.");
            return null;
        }

        var conversionWarnings = new ConversionWarnings();
        CanonicalSequence? canonical;
        try
        {
            canonical = SequenceConverter.Parse(sequence, conversionWarnings, ModHandlingMode);
            if (!canonical.HasValue)
            {
                HandleFailure(ModHandlingMode, "Failed to parse sequence.");
                warning = BuildWarning(conversionWarnings, null);
                return null;
            }
        }
        catch (SequenceConversionException ex)
        {
            HandleFailure(ModHandlingMode, $"Failed to parse sequence: {ex.Message}");
            warning = BuildWarning(conversionWarnings, ex.Message);
            return null;
        }

        if (!IsValidBaseSequence(canonical.Value.BaseSequence, AllowedAminoAcidPattern, MinPeptideLength, MaxPeptideLength))
        {
            HandleFailure(ModHandlingMode, "Invalid base sequence.");
            return null;
        }

        var cleaned = canonical.Value;
        if (ModHandlingMode == SequenceConversionHandlingMode.UsePrimarySequence && cleaned.HasModifications)
        {
            cleaned = cleaned.WithModifications(Array.Empty<CanonicalModification>());
            conversionWarnings.AddWarning("Sequence modifications were removed for prediction.");
        }

        string? serialized;
        try
        {
            serialized = SequenceConverter.Serialize(cleaned, conversionWarnings, ModHandlingMode);
        }
        catch (SequenceConversionException ex)
        {
            HandleFailure(ModHandlingMode, ex.Message);
            warning = BuildWarning(conversionWarnings, ex.Message);
            return null;
        }

        if (serialized == null)
        {
            warning = BuildWarning(conversionWarnings, null);
            return null;
        }

        apiSequence = serialized;
        warning = BuildWarning(conversionWarnings, null);
        return apiSequence;
    }

    #endregion

    protected static ISequenceConverter CreateUnimodConverter(
        UnimodSequenceFormatSchema schema,
        IReadOnlySet<int> allowedUnimodIds)
    {
        var lookup = CreateLookup(allowedUnimodIds);
        var serializer = new UnimodSequenceSerializer(schema, lookup);
        return new SequenceConverter(MzLibSequenceParser.Instance, serializer);
    }

    private static IModificationLookup CreateLookup(IReadOnlySet<int> allowedUnimodIds)
    {
        if (allowedUnimodIds.Count == 0)
        {
            return new UnimodModificationLookup(Enumerable.Empty<Modification>());
        }

        var candidates = Mods.UnimodModifications
            .Where(m => TryGetUnimodId(m, out var id) && allowedUnimodIds.Contains(id))
            .ToList();

        return new UnimodModificationLookup(candidates);
    }

    private static bool TryGetUnimodId(Modification modification, out int id)
    {
        if (modification.Accession?.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase) == true
            && int.TryParse(modification.Accession[7..], out id))
        {
            return true;
        }

        if (modification.DatabaseReference != null)
        {
            foreach (var kvp in modification.DatabaseReference)
            {
                if (!kvp.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }

                if (kvp.Value.Count > 0)
                {
                    var reference = kvp.Value[0]
                        .Replace("UNIMOD:", string.Empty, StringComparison.OrdinalIgnoreCase)
                        .Replace(":", string.Empty);

                    if (int.TryParse(reference, out id))
                    {
                        return true;
                    }
                }
            }
        }

        id = -1;
        return false;
    }

    private static bool IsValidBaseSequence(string baseSequence, string allowedPattern, int minLength, int maxLength)
    {
        return Regex.IsMatch(baseSequence, allowedPattern)
               && baseSequence.Length <= maxLength
               && baseSequence.Length >= minLength;
    }

    private static void HandleFailure(SequenceConversionHandlingMode mode, string message)
    {
        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new ArgumentException(message);
        }
    }

    private static WarningException? BuildWarning(ConversionWarnings warnings, string? additionalMessage)
    {
        var messages = new List<string>();

        if (!string.IsNullOrWhiteSpace(additionalMessage))
        {
            messages.Add(additionalMessage);
        }

        if (warnings.HasIncompatibleItems)
        {
            messages.Add($"Sequence contains unsupported modifications: {string.Join(", ", warnings.IncompatibleItems)}");
        }

        if (warnings.HasErrors)
        {
            messages.AddRange(warnings.Errors);
        }

        if (warnings.HasWarnings)
        {
            messages.AddRange(warnings.Warnings);
        }

        return messages.Count > 0 ? new WarningException(string.Join(" ", messages)) : null;
    }
}
