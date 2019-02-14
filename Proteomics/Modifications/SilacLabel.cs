using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class SilacLabel
    {
        public char OriginalAminoAcid { get; private set; }
        public char AminoAcidLabel { get; private set; }
        public string LabelChemicalFormula { get; private set; }
        public string MassDifference { get; private set; }

        public SilacLabel(char originalAminoAcid, char aminoAcidLabel, string labelChemicalFormula, double massDifference)
        {
            OriginalAminoAcid = originalAminoAcid;
            AminoAcidLabel = aminoAcidLabel;
            LabelChemicalFormula = labelChemicalFormula;
            MassDifference = Math.Round(massDifference, 3).ToString("F3");
            if (massDifference > 0)//if not negative, add a plus
            {
                MassDifference = "+" + MassDifference;
            }
        }
    }
}
