using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Quantification.Interfaces
{
    /// <summary>
    /// Collapse strategies take in a QuantMatrix and reduce the number of columns (samples) by combining them according to some rule.
    /// An example would be collapsing technical replicates into a single sample by averaging their values.
    /// This is a column-wise operation that reduces the number of columns in the matrix.
    /// Incoming matrix has dimensions (p x s) where s = number of samples and p = number of peptide/proteins 
    /// Outgoing matrix has dimensions (p x S) where S = number of collapsed sample, S <= s, and p = number of peptide/proteins 
    /// </summary>
    public interface ICollapseStrategy
    {
        string Name { get; }

        /// <summary>
        /// Example would be to add all technical replicates together by summing their values.
        /// Then add all fractions together by summing their values.
        /// </summary>
        QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>;
    }
}
