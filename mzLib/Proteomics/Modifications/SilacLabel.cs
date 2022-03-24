using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace Proteomics
{
    /// <summary>
    /// Silac labels used to modify unlabeled proteins
    /// </summary>
    public class SilacLabel
    {
        public char OriginalAminoAcid { get; private set; }
        public char AminoAcidLabel { get; private set; }
        public string LabelChemicalFormula { get; private set; }
        public string MassDifference { get; private set; }
        public List<SilacLabel> AdditionalLabels { get; private set; }

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

        public void AddAdditionalSilacLabel(SilacLabel label)
        {
            if (AdditionalLabels == null)
            {
                AdditionalLabels = new List<SilacLabel> { label };
            }
            else
            {
                AdditionalLabels.Add(label);
            }
        }

        /// <summary>
        /// This method exists for conversion of Silac labels, which take double inputs
        /// Although a double object could be saved, it clutters tomls
        /// </summary>
        /// <returns></returns>
        public double ConvertMassDifferenceToDouble()
        {
            string substring = MassDifference.Substring(1);
            double value = Convert.ToDouble(substring, CultureInfo.InvariantCulture);
            if (MassDifference[0] == '-')
            {
                value *= -1;
            }
            return value;
        }

        /// this parameterless constructor needs to exist to read the toml.
        /// if you can figure out a way to get rid of it, feel free...
        /// this is also encountered in MetaMorpheus's "CommonParameters.cs" if you find a solution.
        public SilacLabel()
        {
        }
    }
}
