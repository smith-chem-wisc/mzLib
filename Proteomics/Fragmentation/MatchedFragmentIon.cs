using Chemistry;

namespace Proteomics.Fragmentation
{
    public class MatchedFragmentIon
    {
        public readonly Product TheoreticalProduct;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly int Charge;

        /// <summary>
        /// Constructs a new MatchedFragmentIon given information about a theoretical and an experimental fragment mass spectral peak
        /// </summary>
        public MatchedFragmentIon(Product theoreticalProduct, double experMz, double experIntensity, int charge)
        {
            TheoreticalProduct = theoreticalProduct;
            Mz = experMz;
            Intensity = experIntensity;
            Charge = charge;
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// TODO: Convert to a usable format for output
        /// </summary>
        public override string ToString()
        {
            return TheoreticalProduct.ProductType + TheoreticalProduct.TerminusFragment.FragmentNumber + "+" + Charge + "\t;" + TheoreticalProduct.NeutralMass;
        }
    }
}
