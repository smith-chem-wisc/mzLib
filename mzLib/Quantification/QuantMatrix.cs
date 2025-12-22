using System.Collections.Immutable;
using System.Runtime.CompilerServices;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Double;
using Omics;
using Omics.BioPolymerGroup;


// Make the internals visible to the Test project
[assembly: InternalsVisibleTo("Test")]

namespace Quantification
{

    /// <summary>
    /// The QuantMatrix class represents a matrix structure for quantification data.
    /// It uses the MathNet.Numerics library to handle matrix operations.
    /// The DenseMatrix class from MathNet.Numerics is used to store the quantification values.
    /// This enables really efficient numerical computations and manipulations of the quantification data.
    /// (NOTE: This may or may not be true: The MathNet.Numerics claims that it uses hardware acceleration, but it's mostly
    /// designed for 3D matrices, like those you would encounter when working with computer graphics. It's likely that the DenseMatrix
    /// is slower than some other implementations for small datasets, but becomes more efficient with larger datasets. We'll need to benchmark this later to be sure.)
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class QuantMatrix<T> where T : IEquatable<T>
    {
        internal ImmutableList<T> RowKeys { get; }
        internal ImmutableList<ISampleInfo> ColumnKeys { get; }
        internal IExperimentalDesign ExperimentalDesign { get; }
        internal DenseMatrix Matrix { get; }

        /// <summary>
        /// A matrix that holds quantification values.
        /// Rows represent different types of biological entities (e.g., peptides, proteins) identified by their keys of type T.
        /// Columns represent the different samples in which they were observed
        /// </summary>
        /// <param name="rowKeys"></param>
        /// <param name="columnKeys"></param>
        public QuantMatrix(
            ICollection<T> rowKeys,
            ICollection<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign)
        {
            RowKeys = rowKeys.ToImmutableList();
            ColumnKeys = columnKeys.ToImmutableList();
            Matrix = new DenseMatrix(RowKeys.Count, ColumnKeys.Count);
            ExperimentalDesign = experimentalDesign;
        }

        public void SetRow(T rowKey, double[] values)
        {
            int rowIndex = RowKeys.IndexOf(rowKey);
            if (rowIndex == -1)
            {
                throw new ArgumentException("Row key not found in matrix.");
            }
            if (values.Length != ColumnKeys.Count)
            {
                throw new ArgumentException("Values array length does not match number of columns.");
            }
            for (int colIndex = 0; colIndex < ColumnKeys.Count; colIndex++)
            {
                Matrix[rowIndex, colIndex] = values[colIndex];
            }
        }

        public double[] GetRow(T rowKey)
        {
            int rowIndex = RowKeys.IndexOf(rowKey);
            if (rowIndex == -1)
            {
                throw new ArgumentException("Row key not found in matrix.");
            }
            double[] values = new double[ColumnKeys.Count];
            for (int colIndex = 0; colIndex < ColumnKeys.Count; colIndex++)
            {
                values[colIndex] = Matrix[rowIndex, colIndex];
            }
            return values;
        }

        public double[] GetRow(int rowIndex)
        {
            if (rowIndex < 0 || rowIndex >= RowKeys.Count)
            {
                throw new ArgumentOutOfRangeException("Row index is out of range.");
            }
            double[] values = new double[ColumnKeys.Count];
            for (int colIndex = 0; colIndex < ColumnKeys.Count; colIndex++)
            {
                values[colIndex] = Matrix[rowIndex, colIndex];
            }
            return values;
        }

        public void SetColumn(ISampleInfo columnKey, double[] values)
        {
            int colIndex = ColumnKeys.IndexOf(columnKey);
            if (colIndex == -1)
            {
                throw new ArgumentException("Column key not found in matrix.");
            }
            if (values.Length != RowKeys.Count)
            {
                throw new ArgumentException("Values array length does not match number of rows.");
            }
            for (int rowIndex = 0; rowIndex < RowKeys.Count; rowIndex++)
            {
                Matrix[rowIndex, colIndex] = values[rowIndex];
            }
        }   
    }

    public class PeptideMatrix : QuantMatrix<IBioPolymerWithSetMods>
    {
        public PeptideMatrix(
            ICollection<IBioPolymerWithSetMods> rowKeys,
            ICollection<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign = null) 
            : base(rowKeys, columnKeys, experimentalDesign)
        {
        }
    }

    public class ProteinMatrix : QuantMatrix<IBioPolymerGroup>
    {
        public ProteinMatrix(
            ICollection<IBioPolymerGroup> rowKeys,
            ICollection<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign = null) 
            : base(rowKeys, columnKeys, experimentalDesign)
        {
        }
    }   

    public class SpectralMatchMatrix : QuantMatrix<ISpectralMatch>
    {
        public SpectralMatchMatrix(
            ICollection<ISpectralMatch> rowKeys,
            ICollection<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign = null) 
            : base(rowKeys, columnKeys, experimentalDesign)
        {
        }
    }
}
