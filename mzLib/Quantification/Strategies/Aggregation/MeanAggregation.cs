using System;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>Takes the arithmetic mean. An empty span aggregates to 0.</summary>
    public class MeanAggregation : IAggregationStrategy
    {
        public string Name => "Mean";

        public double Aggregate(ReadOnlySpan<double> values)
        {
            if (values.Length == 0) return 0.0;

            double sum = 0;
            for (int i = 0; i < values.Length; i++) sum += values[i];
            return sum / values.Length;
        }
    }
}
