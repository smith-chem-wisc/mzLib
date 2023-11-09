using Chemistry;
using System.Text;

namespace Proteomics.Fragmentation
{
    public class MatchedFragmentIon
    {
        public readonly Product NeutralTheoreticalProduct;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly int Charge;

        /// <summary>
        /// Constructs a new MatchedFragmentIon given information about a theoretical and an experimental fragment mass spectral peak
        /// </summary>
        public MatchedFragmentIon(ref Product neutralTheoreticalProduct, double experMz, double experIntensity, int charge)
        {
            NeutralTheoreticalProduct = neutralTheoreticalProduct;
            Mz = experMz;
            Intensity = experIntensity;
            Charge = charge;
        }

        public double MassErrorDa
        {
            get
            {
                return Mz.ToMass(Charge) - NeutralTheoreticalProduct.NeutralMass;
            }
        }

        public double MassErrorPpm
        {
            get
            {
                return (MassErrorDa / NeutralTheoreticalProduct.NeutralMass) * 1e6;
            }
        }

        public string Annotation
        {
            get
            {
                StringBuilder sb = new StringBuilder();

                bool containsNeutralLoss = NeutralTheoreticalProduct.NeutralLoss != 0;

                if (containsNeutralLoss)
                {
                    sb.Append("(");
                }

                sb.Append(NeutralTheoreticalProduct.Annotation);

                if (containsNeutralLoss)
                {
                    sb.Append(")");
                }

                sb.Append("+");
                sb.Append(Charge);

                return sb.ToString();
            }
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            // we add the blank space in the tostring because the values are treated like integers and looked up as index in the enum instead of being converted to just string and concatenated
            return NeutralTheoreticalProduct.ProductType + "" + NeutralTheoreticalProduct.FragmentNumber + "+" + Charge + "\t;" + NeutralTheoreticalProduct.NeutralMass;
        }

        public override bool Equals(object obj)
        {
            MatchedFragmentIon other = (MatchedFragmentIon)obj;

            return this.NeutralTheoreticalProduct.Equals(other.NeutralTheoreticalProduct)
                && this.Charge == other.Charge
                && this.Mz == other.Mz
                && this.Intensity == other.Intensity;
        }

        public override int GetHashCode()
        {
            return Mz.GetHashCode();
        }
    }
}
