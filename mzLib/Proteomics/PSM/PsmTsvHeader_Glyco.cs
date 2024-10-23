namespace Proteomics.PSM
{
    public enum LocalizationLevel
    {
        Level1 = 0,
        Level1b = 1,
        Level2 = 2,
        Level3 = 3
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
