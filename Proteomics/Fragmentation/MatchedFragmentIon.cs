using Chemistry;

namespace Proteomics.Fragmentation
{
    public class MatchedFragmentIon
    {
        public readonly Product TheoreticalFragmentIon;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly int Charge;

        /// <summary>
        /// Constructs a new MatchedFragmentIon given information about a theoretical and an experimental fragment mass spectral peak
        /// </summary>
        public MatchedFragmentIon(Product theoreticalFragmentIon, double experMz, double experIntensity, int charge)
        {
            TheoreticalFragmentIon = theoreticalFragmentIon;
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
            return TheoreticalFragmentIon.ProductType + TheoreticalFragmentIon.TerminusFragment.FragmentNumber + "+" + Charge + "\t;" + TheoreticalFragmentIon.NeutralMass;
        }
    }
}
