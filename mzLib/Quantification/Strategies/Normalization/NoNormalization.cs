using Quantification.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification.Strategies
{
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
