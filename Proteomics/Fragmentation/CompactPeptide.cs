using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.Fragmentation
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, FragmentationTerminus fragmentationTerminus, DissociationType dissociationType)
        {
            NeutralTerminusFragment[] NTerminalMasses;
            NeutralTerminusFragment[] CTerminalMasses;
            if (fragmentationTerminus == FragmentationTerminus.Both || fragmentationTerminus == FragmentationTerminus.N)
            {
                //sum of amino acid masses beginning at the N-terminus side. Modification masses are added/subtracted accordingly. The mass of water is subtracted from the total.
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, 0, 1, dissociationType).ToArray();
            }
            if (fragmentationTerminus == FragmentationTerminus.Both || fragmentationTerminus == FragmentationTerminus.C)
            {
                //sum of amino acid masses beginning at the C-terminus side.  Modification masses are added/subtracted accordingly. The mass of water is subtracted from the total.
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, peptideWithSetModifications.Length, -1, dissociationType).ToArray(); //(PeptideWithSetModifications peptide, double prevMass, int residue, int direction, DissociationType dissociationType)
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
            TerminalMasses = NTerminalMasses.Concat(CTerminalMasses);
        }
    }
}