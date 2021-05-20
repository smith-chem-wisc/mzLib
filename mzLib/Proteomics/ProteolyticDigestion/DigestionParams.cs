using Proteomics.Fragmentation;
using System;

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
            int maxModsForPeptides = 2, CleavageSpecificity searchModeType = CleavageSpecificity.Full, FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            bool generateUnlabeledProteinsForSilac = true)
        {
            Protease = ProteaseDictionary.Dictionary[protease];
            MaxMissedCleavages = maxMissedCleavages;
            MinPeptideLength = minPeptideLength;
            MaxPeptideLength = maxPeptideLength;
            MaxModificationIsoforms = maxModificationIsoforms;
            InitiatorMethionineBehavior = initiatorMethionineBehavior;
            MaxModsForPeptide = maxModsForPeptides;
            SearchModeType = searchModeType;
            FragmentationTerminus = fragmentationTerminus;
            RecordSpecificProtease();
            GeneratehUnlabeledProteinsForSilac = generateUnlabeledProteinsForSilac;
        }

        public int MaxMissedCleavages { get; private set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MinPeptideLength { get; private set; }
        public int MaxPeptideLength { get; private set; }
        public int MaxModificationIsoforms { get; private set; }
        public int MaxModsForPeptide { get; private set; }
        public Protease Protease { get; private set; }
        public CleavageSpecificity SearchModeType { get; private set; } //for fast semi and nonspecific searching of proteases
        public FragmentationTerminus FragmentationTerminus { get; private set; } //for fast semi searching of proteases
        public Protease SpecificProtease { get; private set; } //for fast semi and nonspecific searching of proteases
        public bool GeneratehUnlabeledProteinsForSilac { get; private set; } //used to look for unlabeled proteins (in addition to labeled proteins) for SILAC experiments

        public override bool Equals(object obj)
        {
            return obj is DigestionParams a
                && MaxMissedCleavages.Equals(a.MaxMissedCleavages)
                && MinPeptideLength.Equals(a.MinPeptideLength)
                && MaxPeptideLength.Equals(a.MaxPeptideLength)
                && InitiatorMethionineBehavior.Equals(a.InitiatorMethionineBehavior)
                && MaxModificationIsoforms.Equals(a.MaxModificationIsoforms)
                && MaxModsForPeptide.Equals(a.MaxModsForPeptide)
                && Protease.Equals(a.Protease)
                && SearchModeType.Equals(a.SearchModeType)
                && FragmentationTerminus.Equals(a.FragmentationTerminus)
                && GeneratehUnlabeledProteinsForSilac.Equals(a.GeneratehUnlabeledProteinsForSilac);
        }

        public override int GetHashCode()
        {
            return
                MaxMissedCleavages.GetHashCode()
                ^ InitiatorMethionineBehavior.GetHashCode()
                ^ MaxModificationIsoforms.GetHashCode()
                ^ MaxModsForPeptide.GetHashCode();
        }

        public override string ToString()
        {
            return MaxMissedCleavages + "," + InitiatorMethionineBehavior + "," + MinPeptideLength + "," + MaxPeptideLength + ","
                + MaxModificationIsoforms + "," + MaxModsForPeptide + "," + SpecificProtease.Name + "," + SearchModeType + "," + FragmentationTerminus + ","
                + GeneratehUnlabeledProteinsForSilac;
        }

        private void RecordSpecificProtease()
        {
            SpecificProtease = Protease;
            if (SearchModeType == CleavageSpecificity.None) //nonspecific searches, which might have a specific protease
            {
                Protease = FragmentationTerminus == FragmentationTerminus.N ?
                   ProteaseDictionary.Dictionary["singleN"] :
                   ProteaseDictionary.Dictionary["singleC"];
            }
        }
    }
}