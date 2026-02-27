using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Quantification.Strategies
{
    public class SumRollUp : IRollUpStrategy
    {
        public string Name => "Sum Roll-Up";
        public ArrayPool<double> ArrayPool => ArrayPool<double>.Shared;

        public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map)
            where TLow : IEquatable<TLow>
            where THigh : IEquatable<THigh>
        {
            // Create a new QuantMatrix for the rolled-up data
            QuantMatrix<THigh> rolledUpMatrix = new QuantMatrix<THigh>(map.Keys, matrix.ColumnKeys, matrix.ExperimentalDesign);
            // Iterate over each higher-level key in the mapping
            foreach (var kvp in map)
            {
                THigh highKey = kvp.Key;
                List<int> lowIndices = kvp.Value;
                // Initialize an array to hold the summed values for this higher-level key
                double[] summedValues = ArrayPool.Rent(matrix.ColumnCount);
                // Sum the values from the lower-level rows corresponding to this higher-level key
                foreach (int lowIndex in lowIndices)
                {
                    for (int col = 0; col < matrix.ColumnCount; col++)
                    {
                        summedValues[col] += matrix.Matrix[lowIndex, col];
                    }
                }
                // Add the summed values as a new row in the rolled-up matrix
                rolledUpMatrix.SetRow(highKey, summedValues);
                ArrayPool.Return(summedValues, clearArray: true);
            }
            return rolledUpMatrix;
        }
    }
}
