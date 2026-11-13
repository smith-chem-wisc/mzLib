using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.Client;

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
    /// True when the model accepts every UNIMOD modification, whatever <see cref="AllowedUnimodIds"/> holds.
    /// </summary>
    public virtual bool AcceptsAllUnimodModifications => false;

    /// <summary>
    /// Modified residues and termini the model has a token for, as Koina writes them ("M[UNIMOD:35]", "[UNIMOD:737]-"); null = anywhere.
    /// </summary>
    public virtual IReadOnlySet<string>? AllowedModificationTokens => null;

    /// <summary>
    /// UNIMOD ids the N-terminal modification must be one of; null = none required, empty = any allowed one.
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
    /// - Peptide sequences (formatted according to model requirements)
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

    #region Sending a Session's Requests
    /// <summary>
    /// Sends <paramref name="validInputs"/> to Koina and parses the responses, keeping the prediction for
    /// every input Koina predicts even when the model refuses one of them.
    /// </summary>
    /// <param name="validInputs">The inputs to predict, in the order <paramref name="validInputsMask"/> marks them true.</param>
    /// <param name="validInputsMask">The session's mask over all of its inputs. Each input the model refused is set false here.</param>
    /// <param name="responseToPredictions">The model family's parser: responses, and the inputs they answer, to predictions.</param>
    /// <returns>
    /// The predictions, one per input that was predicted, in input order; and for each input the model refused, its
    /// index in <paramref name="validInputsMask"/> and a warning that carries Koina's message, with the
    /// <see cref="KoinaInputRejectedException"/> as its inner exception.
    /// </returns>
    /// <remarks>
    /// Requests go out as they always have: built by ToBatchedRequests, at most MaxNumberOfBatchesPerRequest at a time with
    /// ThrottlingDelayInMilliseconds between groups, under one SessionDeadline. When nothing is refused, that is all that
    /// happens, and the parser sees every response at once, as before.
    ///
    /// Koina answers the whole request with 400 when the model raises on one input, which used to throw away every
    /// prediction in the call. A request refused that way (a <see cref="KoinaInputRejectedException"/>) is now split in
    /// half, each half rebuilt by ToBatchedRequests and sent again, down to single inputs. The halves go out in rounds, each
    /// round throttled the same way and under its own SessionDeadline. The single input still refused is the one the model
    /// will not take. It is reported, and every other input keeps its prediction. Finding one such input among n costs at
    /// most 2*ceil(log2 n) extra requests.
    ///
    /// Nothing else is split. An outage (a 5xx, a timeout, a network error, a <see cref="KoinaServiceException"/>) or a
    /// malformed request throws exactly as before: sending those again in pieces would only multiply the requests.
    ///
    /// Splitting stops, and the first refusal is thrown just as it was before splitting existed, when the model has
    /// refused more inputs than <see cref="MaxRefusedInputs"/> allows, or when ToBatchedRequests does not batch its inputs
    /// in order, MaxBatchSize at a time, so that which inputs a request carried cannot be known.
    ///
    /// In <see cref="SequenceConversionHandlingMode.ThrowException"/> mode a refused input throws an
    /// <see cref="ArgumentException"/> naming it, as <see cref="HandleFailure"/> does for an input refused here.
    /// The other modes leave it unpredicted with the warning; RemoveIncompatibleElements cannot do better, because
    /// Koina does not say which part of the input it refused.
    /// </remarks>
    protected async Task<(List<TModelOutput> Predictions, IReadOnlyDictionary<int, WarningException> RefusedInputs)> RequestPredictionsAsync(
        List<TModelInput> validInputs,
        bool[] validInputsMask,
        Func<IReadOnlyList<string>, List<TModelInput>, List<TModelOutput>> responseToPredictions)
    {
        var requests = ToBatchedRequests(validInputs);
        var outcomes = await SendThrottledAsync(requests);

        var firstRefusal = outcomes.Select(o => o.Refusal).FirstOrDefault(r => r != null);
        if (firstRefusal == null)
            return (responseToPredictions(outcomes.Select(o => o.Response!).ToList(), validInputs), new Dictionary<int, WarningException>());

        // Splitting has to know which inputs each request carried. Every model batches its inputs in order, MaxBatchSize
        // at a time, so request k carries inputs k*MaxBatchSize onward. A model that batches otherwise fails as before.
        if (requests.Count != (validInputs.Count + MaxBatchSize - 1) / MaxBatchSize)
            Rethrow(firstRefusal);

        var answered = new List<(int Start, List<TModelInput> Inputs, string Response)>();
        var stillRefused = new List<(int Start, List<TModelInput> Inputs, KoinaInputRejectedException Refusal)>();
        for (int k = 0; k < outcomes.Count; k++)
        {
            int start = k * MaxBatchSize;
            var inputs = validInputs.GetRange(start, Math.Min(MaxBatchSize, validInputs.Count - start));
            if (outcomes[k].Refusal is { } refusal)
                stillRefused.Add((start, inputs, refusal));
            else
                answered.Add((start, inputs, outcomes[k].Response!));
        }

        var refused = new Dictionary<int, KoinaInputRejectedException>(); // keyed by index in validInputs
        int maxRefused = MaxRefusedInputs(validInputs.Count);
        while (stillRefused.Count > 0)
        {
            var halves = new List<(int Start, List<TModelInput> Inputs)>();
            foreach (var (start, inputs, refusal) in stillRefused)
            {
                if (inputs.Count == 1)
                {
                    if (ModHandlingMode == SequenceConversionHandlingMode.ThrowException)
                        throw new ArgumentException($"Koina refused input {inputs[0]}: {FirstLine(refusal.ServerError)}", refusal);
                    refused[start] = refusal;
                    continue;
                }

                int half = inputs.Count / 2;
                halves.Add((start, inputs.GetRange(0, half)));
                halves.Add((start + half, inputs.GetRange(half, inputs.Count - half)));
            }

            // Each request still refused holds at least one input the model will not take.
            if (refused.Count + halves.Count / 2 > maxRefused)
                Rethrow(firstRefusal);
            if (halves.Count == 0)
                break;

            var halfRequests = halves.Select(h => ToBatchedRequests(h.Inputs)).ToList();
            if (halfRequests.Any(r => r.Count != 1))
                Rethrow(firstRefusal);

            await Task.Delay(ThrottlingDelayInMilliseconds); // the previous group of requests has just gone
            var halfOutcomes = await SendThrottledAsync(halfRequests.Select(r => r[0]).ToList());

            stillRefused = new();
            for (int k = 0; k < halves.Count; k++)
            {
                if (halfOutcomes[k].Refusal is { } refusal)
                    stillRefused.Add((halves[k].Start, halves[k].Inputs, refusal));
                else
                    answered.Add((halves[k].Start, halves[k].Inputs, halfOutcomes[k].Response!));
            }
        }

        // Each answered request is parsed on its own: a split one holds fewer than MaxBatchSize inputs, and the
        // parsers read responses MaxBatchSize inputs at a time.
        var predictions = new List<TModelOutput>(validInputs.Count - refused.Count);
        foreach (var (_, inputs, response) in answered.OrderBy(a => a.Start))
            predictions.AddRange(responseToPredictions(new[] { response }, inputs));

        var refusedInputs = new Dictionary<int, WarningException>();
        for (int i = 0, validIndex = 0; i < validInputsMask.Length; i++)
        {
            if (!validInputsMask[i])
                continue;
            if (refused.TryGetValue(validIndex++, out var refusal))
            {
                validInputsMask[i] = false;
                refusedInputs[i] = new WarningException($"Koina refused this input: {FirstLine(refusal.ServerError)}", refusal);
            }
        }
        return (predictions, refusedInputs);
    }

    /// <summary>
    /// The most inputs one session lets the model refuse before it stops isolating them and fails as before: 1% of
    /// the session's inputs, and never fewer than 10.
    /// </summary>
    /// <remarks>
    /// A few refused inputs are bad inputs. Many look like something wrong for every input, such as a broken model or a
    /// notation the model no longer reads, and isolating each of those costs up to 2*ceil(log2 MaxBatchSize) requests
    /// for nothing. The limit bounds the extra requests in that case to at most 2*limit per round, for at most
    /// ceil(log2 MaxBatchSize) rounds, and keeps the failure loud.
    /// </remarks>
    private static int MaxRefusedInputs(int validInputCount) => Math.Max(10, validInputCount / 100);

    /// <summary>
    /// Sends <paramref name="requests"/> at most MaxNumberOfBatchesPerRequest at a time, with ThrottlingDelayInMilliseconds
    /// between groups, under the SessionDeadline for that many requests. Returns each request's response, or its refusal
    /// when the model refused one of its inputs. Anything else a request throws propagates.
    /// </summary>
    private async Task<List<(string? Response, KoinaInputRejectedException? Refusal)>> SendThrottledAsync(List<Dictionary<string, object>> requests)
    {
        var chunks = requests.Chunk(MaxNumberOfBatchesPerRequest).ToList();
        var outcomes = new List<(string? Response, KoinaInputRejectedException? Refusal)>(requests.Count);
        using var cts = new CancellationTokenSource(SessionDeadline(requests.Count, chunks.Count)); // Bound all of these requests together

        for (int i = 0; i < chunks.Count; i++)
        {
            outcomes.AddRange(await Task.WhenAll(chunks[i].Select(request => SendCatchingRefusalAsync(request, cts.Token))));

            if (i < chunks.Count - 1) // No need to throttle after the last batch
            {
                await Task.Delay(ThrottlingDelayInMilliseconds);
            }
        }
        return outcomes;
    }

    private async Task<(string? Response, KoinaInputRejectedException? Refusal)> SendCatchingRefusalAsync(
        Dictionary<string, object> request, CancellationToken cancellationToken)
    {
        try
        {
            return (await SendInferenceRequestAsync(ModelName, request, cancellationToken), null);
        }
        catch (KoinaInputRejectedException refusal)
        {
            return (null, refusal);
        }
    }

    /// <summary>
    /// Koina's message without the Python traceback that follows it.
    /// </summary>
    private static string FirstLine(string serverError) =>
        serverError.Split('\n').Select(line => line.Trim()).FirstOrDefault(line => line.Length > 0) ?? serverError;

    /// <summary>
    /// Throws <paramref name="refusal"/> again with the stack trace it was first thrown with, as the request would have
    /// failed before splitting existed.
    /// </summary>
    private static void Rethrow(KoinaInputRejectedException refusal) =>
        System.Runtime.ExceptionServices.ExceptionDispatchInfo.Throw(refusal);
    #endregion

    #region Validation and Modification Handling
    /// <summary>
    /// Validates a peptide sequence against model constraints for modifications and basic sequence requirements.
    /// Handles incompatible modifications according to the specified ModHandlingMode.
    /// </summary>
    /// <param name="apiSequence">The sequence as sent to Koina, in the model's API notation (for example UNIMOD).</param>
    /// <returns>
    /// The cleaned sequence in mzLib notation, exactly as predicted: the input's own modification labels, minus any
    /// modification the handling mode removed. Null when the sequence is invalid for this model.
    /// </returns>
    protected virtual string? TryCleanSequence(
        string sequence,
        ISequenceParser? sourceParser,
        out string? apiSequence,
        out WarningException? warning)
    {
        apiSequence = null;
        warning = null;
        var parser = sourceParser ?? SequenceConverter.Parser;

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

        // Resolved here, not in the serializer, so pre-identified UNIMOD ids face the same allow-list check.
        var accepted = new List<CanonicalModification>(cleaned.Modifications.Length);
        // Each accepted modification as the mzLib-notation sequence names it: the input's own label when it gave one, so a
        // peptide predicted from an mzLib sequence keeps that sequence's names; the resolved one otherwise (ProForma).
        var acceptedAsLabelled = new List<CanonicalModification>(cleaned.Modifications.Length);
        var incompatible = new List<CanonicalModification>();
        var refused = new List<string>();
        foreach (var mod in cleaned.Modifications)
        {
            var resolved = SequenceConverter.Serializer.ResolveModification(mod);
            if (resolved.UnimodId is int id && (AcceptsAllUnimodModifications || AllowedUnimodIds.Contains(id)))
            {
                var sent = SequenceConverter.Serializer.Schema.PlaceModification(resolved, cleaned);
                if (DecideModification(cleaned, sent) is { } reason)
                {
                    refused.Add($"{ModificationToken(sent, id, cleaned.BaseSequence)} ({reason})");
                    continue;
                }
                if (AllowedModificationTokens?.Contains(ModificationToken(sent, id, cleaned.BaseSequence)) != false)
                {
                    accepted.Add(sent);
                    acceptedAsLabelled.Add(mod.MzLibId != null || mod.MzLibModification != null ? mod : resolved);
                    continue;
                }
            }
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

        string? removedMessage = null;
        if (refused.Count > 0)
        {
            if (ModHandlingMode != SequenceConversionHandlingMode.RemoveIncompatibleElements)
            {
                var message = $"{ModelName} cannot place {string.Join("; ", refused)}.";
                HandleFailure(ModHandlingMode, message);
                warning = BuildWarning(conversionWarnings, message);
                return null;
            }
            removedMessage = $"Removing modifications {ModelName} cannot place: {string.Join("; ", refused)}.";
        }
        cleaned = cleaned.WithModifications(accepted);

        if (RequiredNTerminalUnimodIds is { } required
            && (cleaned.NTerminalModification?.UnimodId is not int nTermId || (required.Count > 0 && !required.Contains(nTermId))))
        {
            var message = required.Count == 0
                ? "Sequence must carry an N-terminal modification."
                : $"Sequence must carry one of these N-terminal modifications: {string.Join(", ", required.Order().Select(i => $"UNIMOD:{i}"))}.";
            HandleFailure(ModHandlingMode, message);
            warning = BuildWarning(conversionWarnings, message);
            return null;
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

        // The same sequence in mzLib notation, so a caller can rebuild the peptide that was predicted. When incompatible
        // modifications were removed, keep only those the model's serializer accepts on their own: its rule, not a second one.
        var predicted = cleaned.WithModifications(Enumerable.Range(0, accepted.Count)
            .Where(k => ModHandlingMode != SequenceConversionHandlingMode.RemoveIncompatibleElements
                || SequenceConverter.Serialize(cleaned.WithModifications(new[] { accepted[k] }), new ConversionWarnings(),
                    SequenceConversionHandlingMode.ReturnNull) != null)
            .Select(k => acceptedAsLabelled[k]));
        string? cleanedFullSequence = MzLibSequenceSerializer.Instance.Serialize(predicted, conversionWarnings, SequenceConversionHandlingMode.ThrowException);

        warning = BuildWarning(conversionWarnings, removedMessage);
        return cleanedFullSequence;
    }

    /// <summary>
    /// Decides one modification with an allowed UNIMOD id before its token is checked: null to keep it, or why Koina refuses it.
    /// </summary>
    /// <param name="modification">The modification where Koina will be sent it.</param>
    protected virtual string? DecideModification(CanonicalSequence sequence, CanonicalModification modification) => null;

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
        var serializer = new UnimodSequenceSerializer(WithPyroGluOnFirstResidue(schema), lookup);
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
        var serializer = new UnimodSequenceSerializer(WithPyroGluOnFirstResidue(schema), lookup);
        return new SequenceConverter(MzLibSequenceParser.Instance, serializer);
    }

    // Checked live 2026-11-12: every Koina model that takes pyro-Glu predicts it on its residue (E[UNIMOD:27], Q[UNIMOD:28]),
    // and the N-terminal forms give 400, are dropped or are misread.
    private static UnimodSequenceFormatSchema WithPyroGluOnFirstResidue(UnimodSequenceFormatSchema schema) => new(schema.LabelStyle)
    {
        ModOpenBracket = schema.ModOpenBracket,
        ModCloseBracket = schema.ModCloseBracket,
        NTermSeparator = schema.NTermSeparator,
        CTermSeparator = schema.CTermSeparator,
        PyroGluOnFirstResidue = true
    };

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
