using Chemistry;
using System.Text;

namespace Omics.Fragmentation
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
        public MatchedFragmentIon(Product neutralTheoreticalProduct, double experMz, double experIntensity, int charge)
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

        // Doubles are accurate to within 15-17 significant digits. Rounding to 10 decimal places ensures accurate comparison up to 100,000 m/z (non-inclusive)
        internal const int MzDecimalDigits = 10;
        // Rounding to 6 decimal places ensures accurate comparison up to 1,000,000,000 AU (non-inclusive)
        internal const int IntensityDecimalDigits = 6;

        public override bool Equals(object obj)
        {
            return obj is MatchedFragmentIon otherIon && this.Equals(otherIon);
        }


        public bool Equals(MatchedFragmentIon other)
        {
            return this.NeutralTheoreticalProduct.Equals(other.NeutralTheoreticalProduct)
                && this.Charge == other.Charge
                && Math.Round(Mz, MzDecimalDigits) - Math.Round(other.Mz, MzDecimalDigits) == 0
                && Math.Round(Intensity, IntensityDecimalDigits) - Math.Round(other.Intensity, IntensityDecimalDigits) == 0;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(
                NeutralTheoreticalProduct.GetHashCode(),
                Charge.GetHashCode(),
                Math.Round(Mz, MzDecimalDigits).GetHashCode(),
                Math.Round(Intensity, IntensityDecimalDigits).GetHashCode());
        }
    }
}
