using Omics.Digestion;
using Omics.Fragmentation;

namespace Transcriptomics.Digestion
{
    public class RnaDigestionParams : DigestionParamsBase, IDigestionParams, IEquatable<RnaDigestionParams>
    {

        // this parameterless constructor needs to exist to read the toml.
        public RnaDigestionParams() : this("top-down")
        {
        }

        public RnaDigestionParams(string rnase = "top-down", int maxMissedCleavages = 0, int minLength = 3,
            int maxLength = int.MaxValue, int maxModificationIsoforms = 1024, int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full) : base(RnaseDictionary.Dictionary[rnase], maxMissedCleavages, minLength, maxLength, maxModificationIsoforms, maxMods, fragmentationTerminus, searchModeType)
        {

        }


        public Rnase Rnase => (Rnase)DigestionAgent;
        public Rnase SpecificRnase => (Rnase)SpecificDigestionAgent;

        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            var terminus = newTerminus ?? FragmentationTerminus;
            if (SearchModeType == CleavageSpecificity.None)
                return new RnaDigestionParams(SpecificRnase.Name, MaxMissedCleavages, MinLength, MaxLength,
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
                   && SpecificRnase.Equals(other.SpecificRnase)
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
            hash.Add(SpecificRnase);
            hash.Add((int)FragmentationTerminus);
            hash.Add((int)SearchModeType);
            return hash.ToHashCode();
        }

        protected override DigestionAgent GetSingleTerminusAgent(FragmentationTerminus fragmentationTerminus)
        {
            return fragmentationTerminus == FragmentationTerminus.FivePrime
                ? RnaseDictionary.Dictionary["singleN"]
                : RnaseDictionary.Dictionary["singleC"];
        }

        #endregion
    }
}
