using Chemistry;

namespace Proteomics.Fragmentation
{
    public class MatchedFragmentIon
    {
        public readonly TheoreticalFragmentIon TheoreticalFragmentIon;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly double PpmMassError;

        /// <summary>
        /// Constructs a new MatchedFragmentIon given information about a theoretical and an experimental fragment mass spectral peak
        /// </summary>
        public MatchedFragmentIon(TheoreticalFragmentIon theoreticalFragmentIon, double experMz, double experIntensity)
        {
            TheoreticalFragmentIon = theoreticalFragmentIon;
            Mz = experMz;
            Intensity = experIntensity;
            PpmMassError = ((experMz.ToMass(theoreticalFragmentIon.Charge) - theoreticalFragmentIon.Mass) / theoreticalFragmentIon.Mass) * 1e6;
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// TODO: Convert to a usable format for output
        /// </summary>
        public override string ToString()
        {
            return TheoreticalFragmentIon.ProductType.ToString().ToLowerInvariant() + TheoreticalFragmentIon.IonNumber + "+" + TheoreticalFragmentIon.Charge + "\t;" + TheoreticalFragmentIon.Mass;
        }
    }
}
