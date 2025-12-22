using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    public class SumRollUp : IRollUpStrategy
    {
        public string Name => "Sum Roll-Up";
        // Implement roll-up methods here

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
                double[] summedValues = new double[matrix.Matrix.ColumnCount];
                // Sum the values from the lower-level rows corresponding to this higher-level key
                foreach (int lowIndex in lowIndices)
                {
                    var lowRow = matrix.GetRow(lowIndex);
                    for (int sampleIndex = 0; sampleIndex < summedValues.Length; sampleIndex++)
                    {
                        summedValues[sampleIndex] += lowRow[sampleIndex];
                    }
                }
                // Add the summed values as a new row in the rolled-up matrix
                rolledUpMatrix.SetRow(highKey, summedValues);
            }
            return rolledUpMatrix;
        }
    }
}
