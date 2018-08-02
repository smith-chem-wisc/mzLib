namespace Proteomics.Fragmentation
{
    public class NeutralTerminusFragment
    {
        public readonly FragmentationTerminus Terminus;
        public readonly double NeutralMass;
        public readonly int AminoAcidPositionNumber;

        /// <summary>
        ///This object has the neutral mass, the fragment terminus, the description of the netural loss (if any) and the amino acid position of the last amino acid in the fragment next to the break in the peptide backbone. Terminal amino acid is numbered 1.
        /// </summary>
        public NeutralTerminusFragment(FragmentationTerminus terminus, double neutralMass,  int aminoAcidPositionNumber)
        {
            this.Terminus = terminus;
            this.NeutralMass = neutralMass;
            this.AminoAcidPositionNumber = aminoAcidPositionNumber;
        }

        public override bool Equals(object obj)
        {
            return NeutralMass.Equals(((NeutralTerminusFragment) obj).NeutralMass);
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }
    }
}