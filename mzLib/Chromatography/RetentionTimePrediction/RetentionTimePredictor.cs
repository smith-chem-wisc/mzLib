using Omics.SequenceConversion;
using System.Collections.Concurrent;

namespace Chromatography.RetentionTimePrediction;

public abstract class RetentionTimePredictor : IRetentionTimePredictor, IDisposable
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

    /// <summary>
    /// Indicates whether <see cref="PredictCore"/> and <see cref="GetFormattedSequence"/>
    /// are safe to invoke concurrently on the same instance from multiple threads.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Default is <c>false</c></b> — a safe default for any subclass that has not been
    /// explicitly audited for concurrency. When this property is <c>false</c>, the batch
    /// pipeline (<see cref="PredictRetentionTimeEquivalents"/> and
    /// <see cref="StreamRetentionTimeEquivalents"/>) transparently collapses the degree of
    /// parallelism to 1, regardless of the <c>maxThreads</c> argument. This prevents
    /// silently corrupt predictions from races on mutable predictor state (models, caches,
    /// accumulators) without forcing the caller to know the internals of each predictor.
    /// </para>
    /// <para>
    /// <b>Subclasses should override this property to return <c>true</c></b> only when
    /// every invocation of <see cref="PredictCore"/> and <see cref="GetFormattedSequence"/>
    /// either reads exclusively from immutable / readonly state or is otherwise fully
    /// reentrant (no shared mutable fields, no non-reentrant native handles, no locks that
    /// would merely re-serialize threads defeating the parallel path).
    /// </para>
    /// <para>
    /// This gate affects batch parallelism only. Single-call predictions via
    /// <see cref="PredictRetentionTimeEquivalent"/> are unaffected: callers always control
    /// thread use at that entry point.
    /// </para>
    /// </remarks>
    protected virtual bool IsConcurrentPredictionSafe => false;

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
    /// Core prediction logic — called when a peptide has passed all validation.
    /// </summary>
    /// <remarks>
    /// Implementations that are safe to invoke concurrently on the same instance
    /// from multiple threads MUST also override <see cref="IsConcurrentPredictionSafe"/>
    /// to return <c>true</c>. Otherwise the batch pipeline will run this method
    /// single-threaded on this instance, regardless of the caller's <c>maxThreads</c>.
    /// </remarks>
    protected abstract double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null);

    /// <summary>
    /// Parallel producer — emits results as they complete, unordered.
    /// Peptides are partitioned into at most <paramref name="maxThreads"/> range-based chunks
    /// (using ceiling division) so each worker processes a contiguous slice of roughly equal size
    /// (e.g. 100 000 peptides / 10 threads = 10 000 each). When the peptide count is smaller than
    /// <paramref name="maxThreads"/>, the number of partitions equals the peptide count.
    /// Both <see cref="PredictRetentionTimeEquivalents"/> and <see cref="StreamRetentionTimeEquivalents"/>
    /// delegate here. Predictors that support true batched inference (e.g. Chronologer) should
    /// override <see cref="PredictRetentionTimeEquivalents"/> instead of this method.
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
        // Gate: predictors that have not opted in to concurrent prediction are run
        // single-threaded regardless of caller-supplied maxThreads. See
        // IsConcurrentPredictionSafe for the rationale.
        int effectiveThreads = IsConcurrentPredictionSafe ? maxThreads : 1;

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
                var peptideList = peptides.ToList();
                if (peptideList.Count == 0)
                    return;

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
            catch (AggregateException aex)
            {
                foreach (var inner in aex.Flatten().InnerExceptions)
                    producerExceptions.Enqueue(inner);
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
    /// <param name="maxThreads">Degree of parallelism; must be at least 1.</param>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="maxThreads"/> is less than 1.
    /// </exception>
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