using System;
using System.Linq;
using MassSpectrometry;

namespace Proteomics.ProteolyticDigestion
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, TerminusType terminusType, DissociationType dissociationType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, 0, 1, dissociationType).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1, dissociationType).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }
    }
}