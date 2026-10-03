using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.Client;
using Readers.ProForma;

namespace PredictionClients.Koina.AbstractClasses;

public abstract class KoinaModelBase<TModelInput, TModelOutput>
{
    protected static readonly Regex BaseStripper = new(@"\[[^\]]+\]", RegexOptions.Compiled);

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
    /// Unimod modification IDs accepted by the model when converting sequences.
    /// Used by the modification converter layer (not parameter validation).
    /// empty = no modifications are accepted, UNLESS <see cref="AcceptsAllUnimodModifications"/> is true.
    /// </summary>
    public virtual IReadOnlySet<int> AllowedUnimodIds => new HashSet<int>();

    /// <summary>
    /// True when this model accepts any UNIMOD-identified modification instead of restricting to
    /// <see cref="AllowedUnimodIds"/>. Models built via <see cref="CreateUnimodConverterAcceptAll"/>
    /// must override this to true, since an empty <see cref="AllowedUnimodIds"/> otherwise means
    /// "reject every modification".
    /// </summary>
    public virtual bool AcceptsAllUnimodModifications => false;

    /// <summary>
    /// The modified residues and termini the model's preprocessing has a token for, written as Koina writes them
    /// (e.g. "M[UNIMOD:35]", "[UNIMOD:737]-"). An allowed id anywhere else is not a modification the model accepts.
    /// null = any residue or terminus will do for an allowed id.
    /// </summary>
    public virtual IReadOnlySet<string>? AllowedModificationTokens => null;

    /// <summary>
    /// UNIMOD IDs of the N-terminal modification this model requires.
    /// null = no N-terminal modification is required.
    /// empty = an N-terminal modification IS required, and any one the model allows will do
    /// (unlike <see cref="AllowedUnimodIds"/>, where empty means none).
    /// populated = the N-terminal modification must be one of these IDs, each of which must also be allowed.
    /// </summary>
    public virtual IReadOnlySet<int>? RequiredNTerminalUnimodIds => null;

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
    /// - Peptide sequences (formatted according to model requirements): each input's Koina sequence, from the model
    ///   family's GetKoinaSequence, not its ValidatedFullSequence, which is in the input's own format
    /// - Model-specific parameters (e.g., charge states, collision energies, NCE values)
    /// - Any additional metadata required by the specific Koina model
    /// Must ensure that only the validated sequences that meet the model's constraints are included in the batches. 
    /// The total number of sequences across all batches should equal PeptideSequences.Count.
    /// </remarks>
    protected abstract List<Dictionary<string, object>> ToBatchedRequests(List<TModelInput> validInputs);

    /// <summary>
    /// Creates a single batch request dictionary for the Koina API.
    /// Models call this inside their <see cref="ToBatchedRequests"/> loop to build each batch.
    /// </summary>
    protected static Dictionary<string, object> BuildBatchedRequest(int batchIndex, params InputField[] fields)
    {
        var inputs = new List<object>(fields.Length);
        foreach (var f in fields)
        {
            inputs.Add(new
            {
                name = f.Name,
                shape = new[] { f.Data.Length, 1 },
                datatype = f.Datatype,
                data = f.Data
            });
        }
        return new Dictionary<string, object>
        {
            { "id", $"Batch{batchIndex}" },
            { "inputs", inputs }
        };
    }

    /// <summary>
    /// How long a whole prediction session is allowed to take: every batch of it at twice this
    /// model's own benchmarked speed, plus the throttling delay between chunks, rounded up to whole
    /// minutes and never less than one.
    /// </summary>
    /// <remarks>
    /// Hoisted verbatim from the five model bases, which each carried this expression inline -- and
    /// two of them (FragmentIntensityModel, RetentionTimeModel) carried it WITHOUT the
    /// Math.Max(..., 1) the other three had. That divergence was harmless and is preserved as
    /// harmless rather than fixed as a bug: Math.Ceiling of any positive value is already at least
    /// 1, and no shipped model declares a benchmarked time of zero, so the clamp never changed an
    /// outcome. It is kept because it states the floor, and the floor is now load-bearing -- see
    /// below.
    ///
    /// Stated once, it is testable, which five inline copies were not. That matters more than the
    /// de-duplication: this deadline is the line between "Koina stalled" and "our batching estimate
    /// is wrong", and the live Koina tests skip on the first while still failing on the second (see
    /// Test\KoinaTests\KoinaLiveTestFixture.cs). A regression that made this too short would
    /// otherwise be discovered only as a live test quietly reporting Skipped -- which is exactly the
    /// mistake that got "out of memory" removed from KoinaServiceException.ServiceFaultMarkers.
    /// KoinaModelDiscoveryTests.EveryConcreteModel_SessionDeadlineCoversItsOwnBenchmarkedWork
    /// pins it offline instead, for every model.
    ///
    /// The estimate itself, carried over from those copies: two times the benchmarked per-batch
    /// time gives buffer so a healthy run does not hit the deadline, plus the throttling time
    /// between chunks. The benchmark covers the whole of Predict(), so it already includes overhead
    /// beyond the API call. Large requests do not necessarily scale linearly, so this is a rough
    /// estimate chosen as an aggressive upper bound rather than a tight one.
    /// </remarks>
    public TimeSpan SessionDeadline(int batchCount, int chunkCount)
    {
        int minutes = (int)Math.Ceiling(
            (batchCount * 2 * BenchmarkedTimeForOneMaxBatchSizeInMilliseconds
             + ThrottlingDelayInMilliseconds * chunkCount) / 6e4); // 60000ms/min
        return TimeSpan.FromMinutes(Math.Max(minutes, 1));
    }

    /// <summary>
    /// Sends a single batch request to the Koina API and returns the raw JSON response.
    /// Seam for testing: overriding this lets the batched prediction pipeline run against a
    /// canned transport instead of the network.
    /// </summary>
    protected virtual Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
        => HTTP.InferenceRequest(modelName, request, cancellationToken);
    #endregion

    #region Validation and Modification Handling
    /// <summary>
    /// Validates a peptide sequence against model constraints and cleans it, in four steps: separate the
    /// modifications from the residues using the source format's own brackets; check the residues; resolve every
    /// modification and check it is allowed, and that any required one is present; write the cleaned sequence back
    /// in the source format. Incompatible modifications are handled according to <see cref="ModHandlingMode"/>.
    /// </summary>
    /// <param name="sequence">The raw input sequence string, in the format <paramref name="sourceParser"/> (or the
    /// model's default converter parser, when null) understands.</param>
    /// <param name="sourceParser">Parser for this input; null uses the model's own converter parser.</param>
    /// <param name="koinaSequence">The cleaned sequence with every modification resolved, which
    /// <see cref="SerializeKoinaSequence"/> turns into the sequence sent to Koina.</param>
    /// <returns>The cleaned sequence written back in the source format, with the same modifications as
    /// <paramref name="koinaSequence"/> (see <see cref="RetentionTimePredictionInput.ValidatedFullSequence"/>). Null when
    /// the sequence is invalid for this model.</returns>
    protected virtual string? TryCleanSequence(
        string sequence,
        ISequenceParser? sourceParser,
        out CanonicalSequence? koinaSequence,
        out WarningException? warning)
    {
        koinaSequence = null;
        warning = null;
        var parser = sourceParser ?? SequenceConverter.Parser;

        var sourceSerializer = GetSourceSerializer(parser);
        if (sourceSerializer == null)
        {
            var message = $"No sequence serializer is registered for the '{parser.FormatName}' format, " +
                          "so the validated sequence can't be written in it.";
            HandleFailure(ModHandlingMode, message);
            warning = new WarningException(message);
            return null;
        }

        var residues = SeparateResidues(sequence, parser.Schema);
        if (!IsValidBaseSequence(residues, AllowedAminoAcidPattern, MinPeptideLength, MaxPeptideLength))
        {
            var message = $"Invalid base sequence '{residues}': residues must match {AllowedAminoAcidPattern} " +
                          $"and be {MinPeptideLength}-{MaxPeptideLength} long.";
            HandleFailure(ModHandlingMode, message);
            warning = new WarningException(message);
            return null;
        }

        var conversionWarnings = new ConversionWarnings();
        CanonicalSequence? canonical;
        try
        {
            canonical = parser.Parse(sequence, conversionWarnings, ModHandlingMode);
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

        var cleaned = canonical.Value;
        if (ModHandlingMode == SequenceConversionHandlingMode.UsePrimarySequence && cleaned.HasModifications)
        {
            cleaned = cleaned.WithModifications(Array.Empty<CanonicalModification>());
            conversionWarnings.AddWarning("Sequence modifications were removed for prediction.");
        }

        // Resolve every modification here rather than during serialization, so modifications from any source
        // format (pre-identified ProForma UNIMOD:N tokens or mzLib names) pass through the same check.
        var accepted = new List<CanonicalModification>(cleaned.Modifications.Length);
        var acceptedAsParsed = new List<CanonicalModification>(cleaned.Modifications.Length);
        var incompatible = new List<CanonicalModification>();
        foreach (var mod in cleaned.Modifications)
        {
            var resolved = ResolveModification(mod);
            if (resolved.UnimodId is int id && (AcceptsAllUnimodModifications || AllowedUnimodIds.Contains(id))
                && AllowedModificationTokens?.Contains(ModificationToken(resolved, id, cleaned.BaseSequence)) != false)
            {
                accepted.Add(resolved);
                acceptedAsParsed.Add(mod);
            }
            else
                incompatible.Add(resolved);
        }

        if (incompatible.Count > 0)
        {
            foreach (var mod in incompatible)
                conversionWarnings.AddIncompatibleItem(mod.ToString());

            if (ModHandlingMode != SequenceConversionHandlingMode.RemoveIncompatibleElements)
            {
                HandleFailure(ModHandlingMode, $"Sequence contains unsupported modifications: {string.Join(", ", incompatible)}");
                warning = BuildWarning(conversionWarnings, null);
                return null;
            }

            foreach (var mod in incompatible)
                conversionWarnings.AddWarning($"Removing unsupported modification: {mod}");
        }
        var resolvedSequence = cleaned.WithModifications(accepted);

        if (RequiredNTerminalUnimodIds is { } required
            && (resolvedSequence.NTerminalModification?.UnimodId is not int nTermId || (required.Count > 0 && !required.Contains(nTermId))))
        {
            var message = required.Count == 0
                ? "Sequence must carry an N-terminal modification."
                : $"Sequence must carry one of these N-terminal modifications: {string.Join(", ", required.Order().Select(i => $"UNIMOD:{i}"))}.";
            HandleFailure(ModHandlingMode, message);
            warning = BuildWarning(conversionWarnings, message);
            return null;
        }

        // Written from the modifications as parsed rather than as resolved, which carry the lookup's names, and strictly,
        // so it can't silently lose a modification that Koina is sent.
        string? validated;
        try
        {
            validated = sourceSerializer.Serialize(cleaned.WithModifications(acceptedAsParsed), conversionWarnings,
                SequenceConversionHandlingMode.ThrowException);
        }
        catch (SequenceConversionException ex)
        {
            var message = $"Failed to write the cleaned sequence in the {sourceSerializer.FormatName} format: {ex.Message}";
            HandleFailure(ModHandlingMode, message);
            warning = BuildWarning(conversionWarnings, message);
            return null;
        }

        koinaSequence = resolvedSequence;
        warning = BuildWarning(conversionWarnings, null);
        return validated;
    }

    /// <summary>
    /// Serializes a sequence cleaned by <see cref="TryCleanSequence"/> into the sequence sent to Koina, with the
    /// model's own serializer. Returns null with a warning when it can't, and throws in ThrowException mode.
    /// </summary>
    private protected string? SerializeKoinaSequence(CanonicalSequence koinaSequence, out WarningException? warning)
    {
        var conversionWarnings = new ConversionWarnings();
        string? serialized;
        try
        {
            serialized = SequenceConverter.Serialize(koinaSequence, conversionWarnings, ModHandlingMode);
        }
        catch (SequenceConversionException ex)
        {
            HandleFailure(ModHandlingMode, ex.Message);
            warning = BuildWarning(conversionWarnings, ex.Message);
            return null;
        }

        warning = BuildWarning(conversionWarnings, serialized == null ? "Failed to serialize the sequence for Koina." : null);
        return serialized;
    }

    // mzLib and ProForma input is written with Koina's own serializers for those formats, whose lookups hold protein
    // modifications only (the registered ones also hold mzLib's RNA modifications); any other format with the serializer
    // registered for it.
    private static ISequenceSerializer? GetSourceSerializer(ISequenceParser parser) =>
        string.Equals(parser.FormatName, ProteinMzLibSerializer.FormatName, StringComparison.OrdinalIgnoreCase) ? ProteinMzLibSerializer
        : string.Equals(parser.FormatName, ProteinProFormaSerializer.FormatName, StringComparison.OrdinalIgnoreCase) ? ProteinProFormaSerializer
        : SequenceConversionService.Default.GetSerializer(parser.FormatName);

    private static readonly MzLibSequenceSerializer ProteinMzLibSerializer = new(GlobalModificationLookup.ProteinOnly);
    private static readonly ProFormaSequenceSerializer ProteinProFormaSerializer = ProFormaSequenceSerializer.WithLookup(MzLibModificationLookup.ProteinOnly);

    /// <summary>
    /// Returns what is left of a sequence once its modifications are lifted out using the source format's schema:
    /// each complete bracketed span, plus the terminal separator joining a terminal modification. Whether those
    /// modifications are well formed is the parser's call; anything that isn't a complete span stays in the
    /// result for the residue check.
    /// </summary>
    private static string SeparateResidues(string sequence, SequenceFormatSchema schema)
    {
        var residues = new StringBuilder(sequence.Length);
        var nTermSeparator = schema.NTermSeparator;
        var cTermSeparator = schema.CTermSeparator;

        int i = SkipModifications(sequence, 0, schema);
        if (i > 0 && !string.IsNullOrEmpty(nTermSeparator) && sequence.AsSpan(i).StartsWith(nTermSeparator))
            i += nTermSeparator.Length;

        while (i < sequence.Length)
        {
            int afterModifications = SkipModifications(sequence, i, schema);
            if (afterModifications > i)
            {
                i = afterModifications;
                continue;
            }

            if (!string.IsNullOrEmpty(cTermSeparator) && sequence.AsSpan(i).StartsWith(cTermSeparator))
            {
                int modificationsStart = i + cTermSeparator.Length;
                int afterCTerm = SkipModifications(sequence, modificationsStart, schema);
                if (afterCTerm > modificationsStart && afterCTerm == sequence.Length)
                    break;
            }

            residues.Append(sequence[i]);
            i++;
        }

        return residues.ToString();
    }

    /// <summary>
    /// Returns the index just past the complete bracketed modifications starting at <paramref name="start"/>, or
    /// <paramref name="start"/> itself when none starts there.
    /// </summary>
    private static int SkipModifications(string sequence, int start, SequenceFormatSchema schema)
    {
        int i = start;
        while (i < sequence.Length && sequence[i] == schema.ModOpenBracket)
        {
            int depth = 0;
            int end = -1;
            for (int k = i; k < sequence.Length && end < 0; k++)
            {
                if (sequence[k] == schema.ModOpenBracket)
                    depth++;
                else if (sequence[k] == schema.ModCloseBracket && --depth == 0)
                    end = k + 1;
            }

            if (end < 0)
                break;
            i = end;
        }
        return i;
    }

    /// <summary>
    /// Resolves a modification through the model's serializer lookup, merged the same way the serializer enriches
    /// the modifications it resolves. Returns the modification unchanged when it needs no resolution or none matches.
    /// </summary>
    private CanonicalModification ResolveModification(CanonicalModification mod)
    {
        var serializer = SequenceConverter.Serializer;
        if (!serializer.ShouldResolveMod(mod) || serializer.ModificationLookup?.TryResolve(mod) is not { } match)
            return mod;

        return match with
        {
            PositionType = mod.PositionType,
            ResidueIndex = mod.ResidueIndex,
            TargetResidue = mod.TargetResidue ?? match.TargetResidue,
            OriginalRepresentation = mod.OriginalRepresentation
        };
    }

    private static string ModificationToken(CanonicalModification mod, int unimodId, string residues) => mod.PositionType switch
    {
        ModificationPositionType.NTerminus => $"[UNIMOD:{unimodId}]-",
        ModificationPositionType.CTerminus => $"-[UNIMOD:{unimodId}]",
        _ => $"{residues[mod.ResidueIndex!.Value]}[UNIMOD:{unimodId}]"
    };

    protected static IReadOnlySet<int> UnimodIdsOf(IEnumerable<string> modificationTokens) =>
        modificationTokens.Select(token => int.Parse(Regex.Match(token, @"UNIMOD:(\d+)").Groups[1].Value)).ToHashSet();

    #endregion

    protected static ISequenceConverter CreateUnimodConverter(
        UnimodSequenceFormatSchema schema,
        IReadOnlySet<int> allowedUnimodIds)
    {
        var lookup = CreateLookup(allowedUnimodIds);
        var serializer = new UnimodSequenceSerializer(schema, lookup);
        return new SequenceConverter(MzLibSequenceParser.Instance, serializer);
    }

    /// <summary>
    /// Creates a sequence converter that accepts all UNIMOD modifications.
    /// Used by models like ms2pip and AlphaPeptDeep that accept any modification.
    /// </summary>
    protected static ISequenceConverter CreateUnimodConverterAcceptAll(UnimodSequenceFormatSchema schema)
    {
        var allMods = Mods.UnimodModifications.ToList();
        var lookup = new UnimodModificationLookup(allMods);
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

    protected static bool IsValidBaseSequence(string baseSequence, string allowedPattern, int minLength, int maxLength)
    {
        return Regex.IsMatch(baseSequence, allowedPattern)
               && baseSequence.Length <= maxLength
               && baseSequence.Length >= minLength;
    }

    protected static void HandleFailure(SequenceConversionHandlingMode mode, string message)
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

/// <summary>
/// Describes a single input field for a Koina API batch request.
/// </summary>
/// <param name="Name">Input field name as expected by the Koina model (e.g., "peptide_sequences")</param>
/// <param name="Datatype">Koina data type (e.g., "BYTES", "INT32", "FP32")</param>
/// <param name="Data">The data array for this batch. Must be a single batch chunk (already sliced to MaxBatchSize).</param>
public sealed record InputField(string Name, string Datatype, Array Data);
