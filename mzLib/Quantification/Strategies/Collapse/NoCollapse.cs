using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// The NoCollapse strategy does not collapse any samples and returns the original quant matrix.
    /// This is useful for when the user does not want to perform any collapsing and wants to work with the original quant matrix.
    /// In the Quant workflow, a collapse strategy is required to be selected, so this strategy allows the user to bypass collapsing if they choose to do so.
    /// </summary>
    public class NoCollapse : ICollapseStrategy
    {
        public string Name => "No Collapse";
        
        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            return quantMatrix;
        }
    }
}
