using Proteomics.Fragmentation;

namespace Proteomics.ProteolyticDigestion
{
    public class DigestionParams
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this("trypsin")
        {
        }

        public DigestionParams(string protease = "trypsin", int maxMissedCleavages = 2, int minPeptideLength = 7, int maxPeptideLength = int.MaxValue,
            int maxModificationIsoforms = 1024, InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
            int maxModsForPeptides = 2, bool semiProteaseDigestion = false, FragmentationTerminus terminusTypeSemiProtease = FragmentationTerminus.N)
        {
            Protease = ProteaseDictionary.Dictionary[protease];
            MaxMissedCleavages = maxMissedCleavages;
            MinPeptideLength = minPeptideLength;
            MaxPeptideLength = maxPeptideLength;
            MaxModificationIsoforms = maxModificationIsoforms;
            InitiatorMethionineBehavior = initiatorMethionineBehavior;
            MaxModsForPeptide = maxModsForPeptides;
            SemiProteaseDigestion = semiProteaseDigestion;
            TerminusTypeSemiProtease = terminusTypeSemiProtease;
        }

        public int MaxMissedCleavages { get; private set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MinPeptideLength { get; private set; }
        public int MaxPeptideLength { get; private set; }
        public int MaxModificationIsoforms { get; private set; }
        public int MaxModsForPeptide { get; private set; }
        public Protease Protease { get; private set; }
        public bool SemiProteaseDigestion { get; private set; } //for nonspecific searching of proteases
        public FragmentationTerminus TerminusTypeSemiProtease { get; private set; }

        public override bool Equals(object obj)
        {
            DigestionParams a = obj as DigestionParams;

            return a != null
                && this.MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && this.MinPeptideLength.Equals(a.MinPeptideLength)
                && this.MaxPeptideLength.Equals(a.MaxPeptideLength)
                && this.InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && this.MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && this.MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && this.Protease.Equals(a.Protease)
                && this.SemiProteaseDigestion.Equals(a.SemiProteaseDigestion)
                && this.TerminusTypeSemiProtease.Equals(a.TerminusTypeSemiProtease);
        }

        public override int GetHashCode()
        {
            return
                MaxMissedCleavages.GetHashCode()
                ^ InitiatorMethionineBehavior.GetHashCode()
                ^ MaxModificationIsoforms.GetHashCode()
                ^ MaxModsForPeptide.GetHashCode();
        }
    }
}