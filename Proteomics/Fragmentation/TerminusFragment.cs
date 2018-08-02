namespace Proteomics.Fragmentation
{
    public class NeutralTerminusFragment
    {
        public readonly FragmentationTerminus N_or_C_terminus;
        public readonly double NeutralMass;
        public readonly string NeutralLossLabel;
        public readonly int AminoAcidPositionNumber;

        /// <summary>
        ///This object has the neutral mass, the fragment terminus, the description of the netural loss (if any) and the amino acid position of the last amino acid in the fragment next to the break in the peptide backbone. Terminal amino acid is numbered 1.
        /// </summary>
        /// <param name="n_or_c_terminus"></param>
        /// <param name="neutralMass"></param>
        /// <param name="neutralLossLable"></param>
        /// <param name="aminoAcidPositionNumber"></param>
        public NeutralTerminusFragment(FragmentationTerminus n_or_c_terminus, double neutralMass, string neutralLossLable, int aminoAcidPositionNumber)
        {
            this.N_or_C_terminus = n_or_c_terminus;
            this.NeutralMass = neutralMass;
            this.NeutralLossLabel = neutralLossLable;
            this.AminoAcidPositionNumber = aminoAcidPositionNumber;
        }
    }
}