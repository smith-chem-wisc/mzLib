using System;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>Adds the values together.</summary>
    public class SumAggregation : IAggregationStrategy
    {
        public string Name => "Sum";

        public double Aggregate(ReadOnlySpan<double> values)
        {
            double sum = 0;
            for (int i = 0; i < values.Length; i++) sum += values[i];
            return sum;
        }
    }
}
