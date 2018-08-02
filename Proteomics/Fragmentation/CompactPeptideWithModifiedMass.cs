using System;
using System.Collections.Generic;

namespace Proteomics.Fragmentation
{
    [Serializable]
    public class CompactPeptideWithModifiedMass : CompactPeptideBase
    {
        public CompactPeptideWithModifiedMass(CompactPeptideBase cp, double monoisotopicMassIncludingFixedMods)
        {
            CTerminalMasses = cp.CTerminalMasses;
            NTerminalMasses = cp.NTerminalMasses;
            MonoisotopicMassIncludingFixedMods = cp.MonoisotopicMassIncludingFixedMods;
            ModifiedMass = monoisotopicMassIncludingFixedMods;
        }

        public double ModifiedMass { get; set; }

        public void SwapMonoisotopicMassWithModifiedMass()
        {
            double tempDouble = MonoisotopicMassIncludingFixedMods;
            MonoisotopicMassIncludingFixedMods = ModifiedMass;
            ModifiedMass = tempDouble;
        }

        public void CropTerminalMasses(FragmentationTerminus terminusType)
        {
            List<double> tempList = new List<double>();
            double[] masses = terminusType == FragmentationTerminus.N ? NTerminalMasses : CTerminalMasses;
            for (int i = 0; i < masses.Length; i++)
            {
                if (masses[i] < MonoisotopicMassIncludingFixedMods)
                {
                    tempList.Add(masses[i]);
                }
                else if (terminusType == FragmentationTerminus.N)
                {
                    NTerminalMasses = tempList.ToArray();
                    break;
                }
                else
                {
                    CTerminalMasses = tempList.ToArray();
                    break;
                }
            }
        }
    }
}