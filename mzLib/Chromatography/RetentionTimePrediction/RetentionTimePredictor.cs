using Omics.SequenceConversion;
using System.Collections.Concurrent;

namespace Chromatography.RetentionTimePrediction;

public abstract class RetentionTimePredictor : IRetentionTimePredictor
{
    protected static readonly char[] CanonicalAminoAcids = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();

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

    /// <summary>
    /// Core prediction logic - called when peptide passes all validation
    /// </summary>
    protected abstract double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null);

    /// <summary>
    /// Parallel producer — emits results as they complete, unordered.
    /// Peptides are partitioned into <paramref name="maxThreads"/> range-based chunks so each
    /// thread processes a contiguous slice (e.g. 100 000 peptides / 10 threads = 10 000 each).
    /// Both PredictRetentionTimeEquivalents and StreamRetentionTimeEquivalents delegate here.
    /// Predictors that support true batched inference (e.g. Chronologer) should override
    /// PredictRetentionTimeEquivalents instead of this method.
    /// </summary>
    private IEnumerable<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> ProduceResults(
        IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
    {
        var collection = new System.Collections.Concurrent.BlockingCollection<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>();

        Task.Run(() =>
        {
            var peptideList = peptides.ToList();
            if (peptideList.Count == 0)
            {
                collection.CompleteAdding();
                return;
            }
            int rangeSize = Math.Max(1, peptideList.Count / maxThreads);

            Parallel.ForEach(
                Partitioner.Create(0, peptideList.Count, rangeSize),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                range =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var value = PredictRetentionTimeEquivalent(peptideList[i], out var reason);
                        collection.Add((value, peptideList[i], reason));
                    }
                });
            collection.CompleteAdding();
        });

        foreach (var item in collection.GetConsumingEnumerable())
            yield return item;
    }

    /// <summary>
    /// Streams results as they are produced in parallel. Results may be unordered.
    /// </summary>
    public virtual IEnumerable<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> StreamRetentionTimeEquivalents(
        IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
        => ProduceResults(peptides, maxThreads);

    /// <summary>
    /// Default batch implementation — parallel production, materialized result.
    /// Predictors that support true batched inference (e.g. Chronologer) should override this method.
    /// </summary>
    // TODO: Chronologer should override this with a true batched tensor call
    public virtual IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> PredictRetentionTimeEquivalents(
        IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
        => ProduceResults(peptides, maxThreads).ToList();

    /// <summary>
    /// No-op dispose for predictors that hold no unmanaged resources.
    /// Chronologer overrides this to release its TorchSharp model.
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