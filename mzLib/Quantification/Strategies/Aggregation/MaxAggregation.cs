using System;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>Takes the largest value. An empty span aggregates to 0.</summary>
    public class MaxAggregation : IAggregationStrategy
    {
        public string Name => "Max";

        public double Aggregate(ReadOnlySpan<double> values)
        {
            if (values.Length == 0) return 0.0;

            double max = values[0];
            for (int i = 1; i < values.Length; i++)
                if (values[i] > max) max = values[i];
            return max;
        }
    }
}
