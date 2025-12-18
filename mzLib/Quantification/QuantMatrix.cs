using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MathNet.Numerics.LinearAlgebra.Double;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification
{

    /// <summary>
    /// The QuantMatrix class represents a matrix structure for quantification data.
    /// It uses the MathNet.Numerics library to handle matrix operations.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class QuantMatrix<T> where T : IEquatable<T>
    {
        internal List<T> RowKeys { get; }
        internal List<ISampleInfo> ColumnKeys { get; }
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
            List<T> rowKeys,
            List<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign)
        {
            RowKeys = rowKeys;
            ColumnKeys = columnKeys;
            Matrix = new DenseMatrix(RowKeys.Count, ColumnKeys.Count);
            ExperimentalDesign = experimentalDesign;
        }
    }

    public class PeptideMatrix : QuantMatrix<IBioPolymerWithSetMods>
    {
        public PeptideMatrix(
            List<IBioPolymerWithSetMods> rowKeys,
            List<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign = null) 
            : base(rowKeys, columnKeys, experimentalDesign)
        {
        }
    }

    public class ProteinMatrix : QuantMatrix<IBioPolymerGroup>
    {
        public ProteinMatrix(
            List<IBioPolymerGroup> rowKeys,
            List<ISampleInfo> columnKeys,
            IExperimentalDesign experimentalDesign = null) 
            : base(rowKeys, columnKeys, experimentalDesign)
        {
        }
    }   
}
