using System;
using System.Buffers;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Takes the median, averaging the two middle values for an even count.
    /// An empty span aggregates to 0.
    /// </summary>
    public class MedianAggregation : IAggregationStrategy
    {
        public string Name => "Median";

        /// <summary>Sorting needs a mutable copy, so small inputs go on the stack.</summary>
        private const int StackAllocThreshold = 64;

        public double Aggregate(ReadOnlySpan<double> values)
        {
            if (values.Length == 0) return 0.0;
            if (values.Length == 1) return values[0];

            if (values.Length <= StackAllocThreshold)
            {
                Span<double> sorted = stackalloc double[values.Length];
                values.CopyTo(sorted);
                sorted.Sort();
                return MiddleOf(sorted);
            }

            double[] rented = ArrayPool<double>.Shared.Rent(values.Length);
            try
            {
                var sorted = rented.AsSpan(0, values.Length);
                values.CopyTo(sorted);
                sorted.Sort();
                return MiddleOf(sorted);
            }
            finally
            {
                ArrayPool<double>.Shared.Return(rented);
            }
        }

        private static double MiddleOf(Span<double> sorted)
        {
            int mid = sorted.Length / 2;
            return sorted.Length % 2 == 0
                ? (sorted[mid - 1] + sorted[mid]) / 2.0
                : sorted[mid];
        }
    }
}
