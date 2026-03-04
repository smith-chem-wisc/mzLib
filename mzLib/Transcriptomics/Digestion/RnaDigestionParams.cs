using Omics.Digestion;
using Omics.Fragmentation;

namespace Transcriptomics.Digestion
{
    public class RnaDigestionParams : IDigestionParams
    {

        // this parameterless constructor needs to exist to read the toml.
        public RnaDigestionParams() : this("top-down")
        {
        }

        public RnaDigestionParams(string rnase = "top-down", int maxMissedCleavages = 0, int minLength = 3,
            int maxLength = int.MaxValue, int maxModificationIsoforms = 1024, int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both)
        {
            Rnase = RnaseDictionary.Dictionary[rnase];
            MaxMissedCleavages = maxMissedCleavages;
            MinLength = minLength;
            MaxLength = maxLength;
            MaxMods = maxMods;
            MaxModificationIsoforms = maxModificationIsoforms;
            FragmentationTerminus = fragmentationTerminus;
        }

        public int MaxMissedCleavages { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MaxMods { get; set; }
        public DigestionAgent DigestionAgent => Rnase;
        public Rnase Rnase { get; private set; }
        public FragmentationTerminus FragmentationTerminus { get; set; }
        public CleavageSpecificity SearchModeType { get; set; } = CleavageSpecificity.Full;
        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            return newTerminus.HasValue
                ? new RnaDigestionParams(Rnase.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, MaxMods, newTerminus.Value)
                : new RnaDigestionParams(Rnase.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, MaxMods, FragmentationTerminus);
        }
    }
}
