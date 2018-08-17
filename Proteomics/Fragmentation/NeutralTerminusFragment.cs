using Chemistry;
using System;

namespace Proteomics.Fragmentation
{
    /// <summary>
    /// This object has the neutral mass, the fragment terminus, and the amino acid position of the last amino acid in the fragment next to the 
    /// break in the peptide backbone. N-terminal amino acid is numbered 1.
    /// See for Reference the following two web-pages.
    /// http://www.matrixscience.com/help/fragmentation_help.html
    /// http://www.matrixscience.com/help/aa_help.html
    /// </summary>
    [Serializable]
    public class NeutralTerminusFragment
    {
        public readonly FragmentationTerminus Terminus;
        public readonly double NeutralMass;
        public readonly int FragmentNumber;
        public readonly int AminoAcidPosition;
        
        public NeutralTerminusFragment(FragmentationTerminus terminus, double neutralMass, int fragmentNumber, int aminoAcidPosition)
        {
            this.Terminus = terminus;
            this.NeutralMass = (double)ClassExtensions.RoundedDouble(neutralMass);
            this.FragmentNumber = fragmentNumber;
            this.AminoAcidPosition = aminoAcidPosition;
        }

        public override bool Equals(object obj)
        {
            return NeutralMass.Equals(((NeutralTerminusFragment) obj).NeutralMass);
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }

        public override string ToString()
        {
            return "Term: " + Terminus + "; Mass: " + NeutralMass + "; FragNum:" + FragmentNumber + "; AA: " + AminoAcidPosition;
        }
    }
}