using System;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Leaves every sample as its own column. The aggregation strategy is unused, because no group
    /// ever has more than one member.
    /// </summary>
    public class NoCollapse : ICollapseStrategy
    {
        public string Name => "No Collapse";

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix, IAggregationStrategy aggregation)
            where T : IEquatable<T>
        {
            return quantMatrix;
        }
    }
}
