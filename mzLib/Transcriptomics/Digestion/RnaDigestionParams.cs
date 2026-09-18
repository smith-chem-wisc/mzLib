using Omics.Digestion;
using Omics.Fragmentation;

namespace Transcriptomics.Digestion
{
    public class RnaDigestionParams : IDigestionParams, IEquatable<RnaDigestionParams>
    {

        // this parameterless constructor needs to exist to read the toml.
        public RnaDigestionParams() : this("top-down")
        {
        }

        public RnaDigestionParams(string rnase = "top-down", int maxMissedCleavages = 0, int minLength = 3,
            int maxLength = int.MaxValue, int maxModificationIsoforms = 1024, int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full)
        {
            Rnase = RnaseDictionary.Dictionary[rnase];
            MaxMissedCleavages = maxMissedCleavages;
            MinLength = minLength;
            MaxLength = maxLength;
            MaxMods = maxMods;
            MaxModificationIsoforms = maxModificationIsoforms;
            FragmentationTerminus = fragmentationTerminus;
            SearchModeType = searchModeType;

            SpecificRnase = Rnase;
            if (SearchModeType == CleavageSpecificity.None) //nonspecific searches, which might have a specific protease
            {
                Rnase = FragmentationTerminus == FragmentationTerminus.N ?
                   RnaseDictionary.Dictionary["singleN"] :
                   RnaseDictionary.Dictionary["singleC"];
            }
        }

        public int MaxMissedCleavages { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MaxMods { get; set; }
        public DigestionAgent DigestionAgent => Rnase;
        public DigestionAgent SpecificDigestionAgent  => SpecificRnase;
        public Rnase Rnase { get; private set; }
        public DigestionAgent SpecificRnase { get; private set; }
        public FragmentationTerminus FragmentationTerminus { get; set; }
        public CleavageSpecificity SearchModeType { get; set; }
        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            var terminus = newTerminus ?? FragmentationTerminus;
            if (SearchModeType == CleavageSpecificity.None)
                return new RnaDigestionParams(SpecificDigestionAgent.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, MaxMods, terminus, SearchModeType);
            return new RnaDigestionParams(Rnase.Name, MaxMissedCleavages, MinLength, MaxLength,
                MaxModificationIsoforms, MaxMods, terminus, SearchModeType);
        }

        #region Equality

        public override bool Equals(object? obj)
            => obj is RnaDigestionParams rdp && Equals(rdp);

        bool IEquatable<IDigestionParams>.Equals(IDigestionParams? other)
            => other is RnaDigestionParams rdp && Equals(rdp);

        public bool Equals(RnaDigestionParams? other)
        {
            if (other is null) return false;
            return MaxMissedCleavages == other.MaxMissedCleavages
                   && MinLength == other.MinLength
                   && MaxLength == other.MaxLength
                   && MaxModificationIsoforms == other.MaxModificationIsoforms
                   && MaxMods == other.MaxMods
                   && Rnase.Equals(other.Rnase)
                   && SpecificDigestionAgent.Equals(other.SpecificDigestionAgent)
                   && FragmentationTerminus == other.FragmentationTerminus
                   && SearchModeType == other.SearchModeType;
        }

        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(MaxMissedCleavages);
            hash.Add(MinLength);
            hash.Add(MaxLength);
            hash.Add(MaxModificationIsoforms);
            hash.Add(MaxMods);
            hash.Add(Rnase);
            hash.Add(SpecificDigestionAgent);
            hash.Add((int)FragmentationTerminus);
            hash.Add((int)SearchModeType);
            return hash.ToHashCode();
        }

        #endregion

    }
}
