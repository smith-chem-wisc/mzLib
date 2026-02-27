using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Provides a normalization strategy that leaves intensity values unchanged.
    /// </summary>
    /// <remarks>Use this strategy when no normalization is required for quantitative data. The input matrix
    /// is returned without modification, preserving original values. This can be useful for baseline comparisons or
    /// when normalization is handled externally.</remarks>
    public class NoNormalization : INormalizationStrategy
    {
        public string Name => "No Normalization";

        public QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            // Return the input matrix as is, without any normalization
            return quantMatrix;
        }
    }
}
