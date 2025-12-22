using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Quantification.Interfaces;

namespace Quantification.Strategies.Collapse
{
    public class NoCollapse : ICollapseStrategy
    {
        public string Name => "No Collapse";
        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            return quantMatrix;
        }
    }
}
