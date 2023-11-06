using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics.PSM
{
    public enum LocalizationLevel
    {
        Level1,
        Level1b,
        Level2,
        Level3
    }
    public class PsmTsvHeader_Glyco
    {
        public const string GlycanMass = "GlycanMass";
        public const string GlycanComposition = "Plausible GlycanComposition";
        public const string GlycanStructure = "Plausible GlycanStructure";
        public const string GlycanLocalizationLevel = "GlycanLocalizationLevel";
        public const string LocalizedGlycan = "Localized Glycans with Peptide Site Specific Probability";
    }
}
