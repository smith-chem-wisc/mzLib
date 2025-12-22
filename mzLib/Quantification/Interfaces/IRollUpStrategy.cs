using MassSpectrometry;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics.BioPolymerGroup;
using System.Numerics;

namespace Quantification.Interfaces
{
    /// <summary>
    /// Roll-up strategies take in lower-level quantification data of type TLow (e.g., PSMs or peptides) and aggregate to a 
    /// higher level of type THigh (e.g., peptides or proteins).
    /// This is a row-wise operation that reduces the number of rows in the matrix.
    /// Incoming matrix has dimensions (p x s) where p = number of lower-level entities (TLow) and s = number of samples.
    /// Outgoing matrix has dimensions (P x s) where P = number of higher-level entities (THigh), P <= p, and s = number of samples.
    /// Roll-up requires two things
    ///  1) A Quant Matrix 
    ///  2) A mapping from higher-level entities to the indices of their corresponding lower-level entities. 
    ///     e.g., Dictionary<IBioPolymerWithSetMods, List<int>> mapping peptides to their PSM indices in the QuantMatrix

    /// </summary>
    public interface IRollUpStrategy
    {
        string Name { get; }

        /// <summary>
        /// Aggregates the values of the specified matrix into a new matrix with higher-level keys, according to the
        /// provided mapping.
        /// </summary>
        /// <remarks>The aggregation is performed based on the row indices specified in the mapping. The output matrix will have
        /// one row for every key in the mapping dictionary, with values computed by aggregating the corresponding rows </remarks>
        /// <typeparam name="TLow">The type of the elements in the input matrix. Must implement <see cref="IEquatable{TLow}"/>.</typeparam>
        /// <typeparam name="THigh">The type of the higher-level keys in the resulting matrix. Must implement <see cref="IEquatable{THigh}"/>.</typeparam>
        /// <param name="matrix">The input matrix whose values are to be aggregated.</param>
        /// <param name="map">A dictionary that maps each higher-level key to a list of column indices in the input matrix to aggregate
        /// under that key. Each key in the dictionary represents a group in the resulting matrix.</param>
        /// <returns>A new matrix of type <see cref="QuantMatrix{THigh}"/> containing the aggregated values grouped by the
        /// specified higher-level keys.</returns>
        public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map) 
            where TLow : IEquatable<TLow>
            where THigh: IEquatable<THigh>;
    }
}
