using System;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Averages the three largest values, the "Top3" estimator of Silva et al. (2006),
    /// doi:10.1074/mcp.M500230-MCP200.
    ///
    /// Rolling peptides up to a protein by summing makes a protein's intensity depend on how many of
    /// its peptides were identified, which is a property of the search as much as of the sample.
    /// Averaging every peptide removes that bias but lets one poorly-ionizing peptide drag the
    /// estimate down. Top3 takes the three best-responding peptides instead, which is the basis for
    /// comparing abundance across different proteins rather than one protein across samples.
    ///
    /// Fewer than three values average what is there, rather than dividing by three -- dividing by
    /// three would report a protein seen twice as smaller than the same protein seen three times,
    /// which is the bias Top3 exists to avoid. An empty span aggregates to 0.
    ///
    /// NOT recommended as a <see cref="QuantificationParameters.CollapseAggregationStrategy"/>. This
    /// is an estimator of a protein's abundance from its peptides, so it only means anything when the
    /// values being combined are peptides of one protein. Collapsing fractions or technical replicates
    /// combines repeated measurements of the SAME thing, where taking the three largest is a bias, not
    /// an estimate -- and the collapse path does not exclude zeros first, so unobserved samples would
    /// count as values too. Use it through <see cref="Top3RollUp"/>.
    /// </summary>
    public class Top3Aggregation : IAggregationStrategy
    {
        private const int N = 3;

        public string Name => "Top 3";

        public double Aggregate(ReadOnlySpan<double> values)
        {
            if (values.Length == 0) return 0.0;

            if (values.Length <= N)
            {
                double total = 0.0;
                foreach (double value in values)
                {
                    total += value;
                }
                return total / values.Length;
            }

            // Only the top three matter, so track them directly rather than sorting the whole span.
            // Buffers from ArrayPool are not zeroed, and negative intensities are not meaningful, so
            // the running maxima start at the first three values rather than at 0.
            double first = double.NegativeInfinity;
            double second = double.NegativeInfinity;
            double third = double.NegativeInfinity;

            foreach (double value in values)
            {
                if (value > first)
                {
                    third = second;
                    second = first;
                    first = value;
                }
                else if (value > second)
                {
                    third = second;
                    second = value;
                }
                else if (value > third)
                {
                    third = value;
                }
            }

            return (first + second + third) / N;
        }
    }
}
