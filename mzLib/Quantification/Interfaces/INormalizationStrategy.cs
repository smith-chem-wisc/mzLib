using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification.Interfaces
{


    /// <summary>
    /// Takes in a QuantMatrix and outputs a new QuantMatrix with normalized values. 
    /// Normalization does not change the dimensions of the matrix.
    /// Normalization generally operates on columns (samples) of the matrix, comparing one sample/channel to another.
    /// Strategies to implement include:
    /// Median, Mean, Quantile, None, Log Fold Change (see FlashLFQ), etc.
    /// </summary>
    public interface INormalizationStrategy
    {
        string Name { get; }

        /// <summary>
        /// Normalize intensities in any type of QuantMatrix.
        /// Creates a new matrix with the same structure but normalized values.
        /// </summary>
        /// <typeparam name="T">The type of entity in the matrix (e.g., IBioPolymerWithSetMods, IBioPolymerGroup, ISpectralMatch)</typeparam>
        /// <param name="quantMatrix">The matrix to normalize</param>
        /// <returns>A new normalized QuantMatrix of the same type</returns>
        QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>;

    }

}
