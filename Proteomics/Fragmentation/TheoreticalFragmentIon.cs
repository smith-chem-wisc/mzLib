using Chemistry;

namespace Proteomics.Fragmentation
{
    public class TheoreticalFragmentIon
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="charge"></param>
        /// <param name="productType"></param>
        /// <param name="terminusFragment"></param>
        public TheoreticalFragmentIon(int charge, ProductType productType, NeutralTerminusFragment terminusFragment)
        {
            TerminusFragment = terminusFragment;
            Charge = charge;
            ProductType = productType;
            Mz = TerminusFragment.NeutralMass.ToMz(charge);
        }

        public readonly int Charge;
        public readonly double Mz;
        public readonly ProductType ProductType;
        public readonly NeutralTerminusFragment TerminusFragment;

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + TerminusFragment.AminoAcidPositionNumber + ";" + Mz + "-" + TerminusFragment.NeutralLossLabel;
        }
    }
}