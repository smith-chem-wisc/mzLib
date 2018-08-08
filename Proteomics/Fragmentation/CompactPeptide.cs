using System;
using System.Linq;
using Proteomics.ProteolyticDigestion;

namespace Proteomics.Fragmentation
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, FragmentationTerminus fragmentationTerminus)
        {
            TerminalMasses = ComputeNeutralTerminusFragments(peptideWithSetModifications, fragmentationTerminus).ToArray();
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }
    }
}