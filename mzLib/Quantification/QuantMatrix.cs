using System.Collections.Immutable;
using System.Runtime.CompilerServices;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;


// Make the internals visible to the Test project
[assembly: InternalsVisibleTo("Test")]

namespace Quantification
{

    /// <summary>
    /// The QuantMatrix class represents a matrix structure for quantification data.
    /// Rows represent different types of biological entities (e.g., peptides, proteins) identified by their keys of type T.
    /// Columns represent the different samples in which they were observed.
    /// Uses a raw double[,] array in row-major layout for cache-friendly access and minimal overhead.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class QuantMatrix<T> where T : IEquatable<T>
    {
        internal ImmutableList<T> RowKeys { get; }
        internal ImmutableList<ISampleInfo> ColumnKeys { get; }
        internal IExperimentalDesign ExperimentalDesign { get; }
        internal double[,] Matrix { get; }
        internal int RowCount { get; }
        internal int ColumnCount { get; }

        private readonly Dictionary<T, int> _rowIndex;
        private readonly Dictionary<ISampleInfo, int> _colIndex;

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
            RowCount = RowKeys.Count;
            ColumnCount = ColumnKeys.Count;
            Matrix = new double[RowCount, ColumnCount];
            ExperimentalDesign = experimentalDesign;

            _rowIndex = new Dictionary<T, int>(RowCount);
            for (int i = 0; i < RowCount; i++)
                _rowIndex[RowKeys[i]] = i;

            _colIndex = new Dictionary<ISampleInfo, int>(ColumnCount);
            for (int i = 0; i < ColumnCount; i++)
                _colIndex[ColumnKeys[i]] = i;
        }

        public void SetRow(T rowKey, double[] values)
        {
            if (!_rowIndex.TryGetValue(rowKey, out int rowIdx))
            {
                throw new ArgumentException("Row key not found in matrix.");
            }
            if (values.Length < ColumnCount)
            {
                throw new ArgumentException("Values array length is smaller than number of columns."); // Due to use of array pools, array may be larger than column count, but not smaller
            }
            for (int colIndex = 0; colIndex < ColumnCount; colIndex++)
            {
                Matrix[rowIdx, colIndex] = values[colIndex];
            }
        }

        public double[] GetRow(T rowKey)
        {
            if (!_rowIndex.TryGetValue(rowKey, out int rowIdx))
            {
                throw new ArgumentException("Row key not found in matrix.");
            }
            double[] values = new double[ColumnCount];
            for (int colIndex = 0; colIndex < ColumnCount; colIndex++)
            {
                values[colIndex] = Matrix[rowIdx, colIndex];
            }
            return values;
        }

        public double[] GetRow(int rowIndex)
        {
            if (rowIndex < 0 || rowIndex >= RowCount)
            {
                throw new ArgumentOutOfRangeException("Row index is out of range.");
            }
            double[] values = new double[ColumnCount];
            for (int colIndex = 0; colIndex < ColumnCount; colIndex++)
            {
                values[colIndex] = Matrix[rowIndex, colIndex];
            }
            return values;
        }

        public void SetColumn(ISampleInfo columnKey, double[] values)
        {
            if (!_colIndex.TryGetValue(columnKey, out int colIdx))
            {
                throw new ArgumentException("Column key not found in matrix.");
            }
            if (values.Length < RowCount)
            {
                throw new ArgumentException("Values array length is smaller than number of rows.");
            }
            for (int rowIndex = 0; rowIndex < RowCount; rowIndex++)
            {
                Matrix[rowIndex, colIdx] = values[rowIndex];
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
