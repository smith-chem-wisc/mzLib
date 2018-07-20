using System;
using System.Linq;
using MassSpectrometry;
using Proteomics.Fragmentation;

namespace Proteomics.ProteolyticDigestion
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, FragmentationTerminus fragmentationTerminus, DissociationType dissociationType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (fragmentationTerminus == FragmentationTerminus.Both || fragmentationTerminus == FragmentationTerminus.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, 0, 1, dissociationType).ToArray();
            }
            if (fragmentationTerminus == FragmentationTerminus.Both || fragmentationTerminus == FragmentationTerminus.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1, dissociationType).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }
    }
}