using Quantification.Interfaces;

namespace Quantification.Strategies
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
