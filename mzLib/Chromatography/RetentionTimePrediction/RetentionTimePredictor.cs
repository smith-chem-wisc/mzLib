using Omics.SequenceConversion;
using System.Collections.Concurrent;

namespace Chromatography.RetentionTimePrediction;

public abstract class RetentionTimePredictor : IRetentionTimePredictor
{
    /// <summary>
    /// The 20 standard proteinogenic amino-acid one-letter codes. Stored as a
    /// <see cref="HashSet{T}"/> so per-character membership tests in hot validation
    /// loops are amortized O(1) rather than an O(20) linear scan.
    /// </summary>
    protected static readonly HashSet<char> CanonicalAminoAcids = new("ACDEFGHIKLMNPQRSTVWY");

    /// <summary>
    /// Name/identifier for this predictor (e.g., "SSRCalc3", "Chronologer")
    /// </summary>
    public abstract string PredictorName { get; }

    /// <summary>
    /// Gets the separation type this predictor is designed for
    /// </summary>
    public abstract SeparationType SeparationType { get; }

    /// <summary>
    /// Determine what to do if we encounter incompatible modifications in an RT predictor where mods are important. 
    /// </summary>
    public virtual SequenceConversionHandlingMode SequenceHandlingMode { get; }
    protected virtual int MinSequenceLength => 7;
    protected virtual int MaxSequenceLength => int.MaxValue;

    protected RetentionTimePredictor(SequenceConversionHandlingMode sequenceHandlingMode = SequenceConversionHandlingMode.UsePrimarySequence)
    {
        SequenceHandlingMode = sequenceHandlingMode;
    }

    public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        if (peptide == null)
            throw new ArgumentNullException(nameof(peptide));

        if (!ValidateBasicConstraints(peptide, out failureReason))
            return null;

        string? sequenceToPredict = GetFormattedSequence(peptide, out failureReason);
        if (sequenceToPredict == null)
            return null;
        try
        {

            return PredictCore(peptide, sequenceToPredict);
        }
        catch (Exception)
        {
            failureReason = RetentionTimeFailureReason.PredictionError;
            return null;
        }
    }

    /// <inheritdoc cref="PredictRetentionTimeEquivalent"/>
    [Obsolete("Use PredictRetentionTimeEquivalent instead.")]
    public double? PredictRetentionTime(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        => PredictRetentionTimeEquivalent(peptide, out failureReason);

    /// <summary>
    /// Core prediction logic - called when peptide passes all validation
    /// </summary>
    /// <remarks>
    /// Implementations are responsible for their own thread-safety. The batch pipeline
    /// in <see cref="PredictRetentionTimeEquivalents"/> may invoke this method
    /// concurrently from multiple worker threads when <c>maxThreads &gt; 1</c>.
    /// </remarks>
    protected abstract double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null);

    /// <summary>
    /// Parallel producer — emits results as they complete, unordered.
    /// Peptides are partitioned into at most <paramref name="maxThreads"/> range-based chunks
    /// (using ceiling division) so each worker processes a contiguous slice of roughly equal size
    /// (e.g. 100 000 peptides / 10 threads = 10 000 each). When the peptide count is smaller than
    /// <paramref name="maxThreads"/>, the number of partitions equals the peptide count.
    /// <see cref="PredictRetentionTimeEquivalents"/> delegates here. Predictors that support
    /// true batched inference (e.g. Chronologer) should override
    /// <see cref="PredictRetentionTimeEquivalents"/> instead of this method.
    /// </summary>
    /// <remarks>
    /// Exception handling contract:
    /// <list type="bullet">
    ///   <item>Per-peptide prediction errors are already absorbed by
    ///     <see cref="PredictRetentionTimeEquivalent"/>, which returns <c>null</c> with a
    ///     <see cref="RetentionTimeFailureReason"/> — those never reach this method.</item>
    ///   <item>Structural faults (enumerating <paramref name="peptides"/>, partitioning, adding to
    ///     the blocking collection, etc.) are captured and surfaced to the caller as an
    ///     <see cref="AggregateException"/> after the enumerator drains. This guarantees the
    ///     consumer cannot deadlock on <see cref="BlockingCollection{T}.GetConsumingEnumerable"/>
    ///     when the background producer faults, and ensures failures are never silently swallowed.</item>
    ///   <item><see cref="BlockingCollection{T}.CompleteAdding"/> is always invoked exactly once,
    ///     in a <c>finally</c> block, regardless of success or failure.</item>
    /// </list>
    /// </remarks>
    /// <exception cref="AggregateException">
    /// Thrown from the enumerator after all successfully-produced items have been yielded,
    /// when the background producer encountered one or more unhandled exceptions.
    /// </exception>
    private IEnumerable<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> ProduceResults(
        IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
    {
        int effectiveThreads = maxThreads < 1 ? 1 : maxThreads;

        // Fast path: when there is no parallelism to exploit, bypass BlockingCollection +
        // Task.Run + Partitioner + Parallel.ForEach entirely. For small/medium batches this
        // avoids a measurable per-item synchronization cost (Monitor.Enter/Exit inside
        // BlockingCollection.Add plus the producer/consumer handshake) and also streams the
        // source enumerable lazily rather than materializing it into a List<T> up front.
        if (effectiveThreads <= 1)
        {
            foreach (var peptide in peptides)
            {
                var value = PredictRetentionTimeEquivalent(peptide, out var reason);
                yield return (value, peptide, reason);
            }
            yield break;
        }

        // Materialize on the caller's thread so exceptions from enumerating non-rewindable
        // or thread-affine sources surface synchronously, and the source is enumerated
        // exactly once.
        var peptideList = peptides.ToList();
        if (peptideList.Count == 0)
            yield break;

        // BlockingCollection wraps a SemaphoreSlim (finalizable). Use a `using` declaration
        // so it's disposed when the iterator state machine is disposed — i.e. when the caller's
        // foreach drains or the enumerator is otherwise disposed — rather than leaking one
        // finalizable instance per batch call.
        using var collection = new BlockingCollection<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>();
        var producerExceptions = new ConcurrentQueue<Exception>();

        Task.Run(() =>
        {
            try
            {
                // Ceiling division so we produce at most `effectiveThreads` partitions.
                // Floor division would yield rangeSize=1 whenever peptideList.Count < 2*effectiveThreads,
                // creating one partition per peptide and adding needless scheduling overhead.
                int rangeSize = Math.Max(1, (peptideList.Count + effectiveThreads - 1) / effectiveThreads);

                Parallel.ForEach(
                    Partitioner.Create(0, peptideList.Count, rangeSize),
                    new ParallelOptions { MaxDegreeOfParallelism = effectiveThreads },
                    range =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            var value = PredictRetentionTimeEquivalent(peptideList[i], out var reason);
                            collection.Add((value, peptideList[i], reason));
                        }
                    });
            }
            catch (Exception ex)
            {
                producerExceptions.Enqueue(ex);
            }
            finally
            {
                // Guarantees the consumer's GetConsumingEnumerable() can terminate even if
                // the producer task faulted before finishing its work.
                collection.CompleteAdding();
            }
        });

        foreach (var item in collection.GetConsumingEnumerable())
            yield return item;

        if (!producerExceptions.IsEmpty)
        {
            throw new AggregateException(
                "One or more retention time prediction producer threads failed.",
                producerExceptions);
        }
    }

    /// <summary>
    /// Default batch implementation — parallel production, materialized result.
    /// Predictors that support true batched inference (e.g. Chronologer) should override this method.
    /// </summary>
    /// <remarks>
    /// <b>Result order is not guaranteed.</b> Items in the returned list reflect worker
    /// completion order, not the order of <paramref name="peptides"/>. Each tuple carries
    /// its source <c>Peptide</c> reference, so callers should pair predictions to inputs
    /// via the tuple element rather than by index. Do not assume
    /// <c>results[i].Peptide == peptides.ElementAt(i)</c>.
    /// </remarks>
    /// <param name="maxThreads">Degree of parallelism; values less than 1 are clamped to 1.</param>
    /// <exception cref="AggregateException">
    /// Thrown if the background producer encountered one or more unhandled exceptions while
    /// enumerating or partitioning the input. Per-peptide prediction errors are reported via the
    /// <see cref="RetentionTimeFailureReason"/> element of each tuple, not thrown.
    /// </exception>
    public virtual IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> PredictRetentionTimeEquivalents(
        IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
    {
        if (maxThreads < 1)
            maxThreads = 1;
        return ProduceResults(peptides, maxThreads).ToList();
    }

    /// <summary>
    /// Releases resources held by the predictor. The base implementation is a no-op
    /// for predictors that hold no unmanaged state. Derived classes that own
    /// unmanaged resources (e.g. a TorchSharp model) MUST use <c>override</c> —
    /// not <c>new</c> / method hiding — so that <see cref="Dispose"/> invoked
    /// through a base-class or <see cref="IDisposable"/> reference still reaches
    /// the derived cleanup via virtual dispatch.
    /// </summary>
    public virtual void Dispose() { }

    /// <summary>
    /// Format the peptide sequence appropriately for this predictor
    /// </summary>
    public abstract string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason);

    /// <summary>
    /// Validate basic constraints (length, amino acid content, etc.) and NOT modifications
    /// </summary>
    protected virtual bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;

        if (string.IsNullOrEmpty(peptide.BaseSequence))
        {
            failureReason = RetentionTimeFailureReason.EmptySequence;
            return false;
        }

        if (peptide.BaseSequence.Length < MinSequenceLength)
        {
            failureReason = RetentionTimeFailureReason.SequenceTooShort;
            return false;
        }

        if (peptide.BaseSequence.Length > MaxSequenceLength)
        {
            failureReason = RetentionTimeFailureReason.SequenceTooLong;
            return false;
        }

        return true;
    }
}