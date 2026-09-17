using System;
using Omics.Digestion;
using Omics.Fragmentation;

namespace Proteomics.ProteolyticDigestion 
{
    public class DigestionParams : DigestionParamsBase, IDigestionParams, IEquatable<DigestionParams>
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this("trypsin")
        {
        }

        public DigestionParams(string protease = "trypsin", int maxMissedCleavages = 2, int minPeptideLength = 7, int maxPeptideLength = int.MaxValue,
            int maxModificationIsoforms = 1024, InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
            int maxModsForPeptides = 2, CleavageSpecificity searchModeType = CleavageSpecificity.Full, FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            bool generateUnlabeledProteinsForSilac = true, bool keepNGlycopeptide = false, bool keepOGlycopeptide = false)
            : base(ProteaseDictionary.Dictionary[protease], maxMissedCleavages, minPeptideLength, maxPeptideLength,
                  maxModificationIsoforms, maxModsForPeptides, fragmentationTerminus, searchModeType)
        {
            InitiatorMethionineBehavior = initiatorMethionineBehavior;

            GeneratehUnlabeledProteinsForSilac = generateUnlabeledProteinsForSilac;
            KeepNGlycopeptide = keepNGlycopeptide;
            KeepOGlycopeptide = keepOGlycopeptide;
        }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }

        /// <summary>
        /// The protease the caller named. It is the same object as <see cref="Protease"/> except when
        /// <see cref="SearchModeType"/> is None: then <see cref="Protease"/> becomes singleN (FragmentationTerminus N) or
        /// singleC (anything else) to make the non-specific seeds, and this keeps the named protease, whose sites still
        /// limit missed cleavages and are written to settings and output.
        /// </summary>
        public Protease SpecificProtease => (Protease)SpecificDigestionAgent;
        public bool GeneratehUnlabeledProteinsForSilac { get; private set; } //used to look for unlabeled proteins (in addition to labeled proteins) for SILAC experiments
        public bool KeepNGlycopeptide { get; private set; }
        public bool KeepOGlycopeptide { get; private set; }

        #region Properties overridden by more generic interface

        /// <summary>
        /// The protease digestion runs: the one the caller named, or singleN/singleC when <see cref="SearchModeType"/> is
        /// None (see <see cref="SpecificProtease"/>).
        /// </summary>
        public Protease Protease => (Protease)DigestionAgent;

        public int MinPeptideLength
        {
            get => MinLength;
            set => MinLength = value;
        }

        public int MaxPeptideLength
        {
            get => MaxLength;
            set => MaxLength = value;
        }

        public int MaxModsForPeptide
        {
            get => MaxMods;
            set => MaxMods = value;
        }

        #endregion

        #region Equality

        public override bool Equals(object? obj) 
            => obj is DigestionParams dp && Equals(dp);

        bool IEquatable<IDigestionParams>.Equals(IDigestionParams? other)
            => other is DigestionParams dp && Equals(dp);

        public bool Equals(DigestionParams? other)
        {
            if (other is null) return false;
            return MaxMissedCleavages == other.MaxMissedCleavages
                   && MinLength == other.MinLength
                   && MaxLength == other.MaxLength
                   && InitiatorMethionineBehavior == other.InitiatorMethionineBehavior
                   && MaxModificationIsoforms == other.MaxModificationIsoforms
                   && MaxMods == other.MaxMods
                   && Protease.Equals(other.Protease)
                   && SearchModeType == other.SearchModeType
                   && FragmentationTerminus == other.FragmentationTerminus
                   && SpecificProtease.Equals(other.SpecificProtease)
                   && GeneratehUnlabeledProteinsForSilac == other.GeneratehUnlabeledProteinsForSilac
                   && KeepNGlycopeptide == other.KeepNGlycopeptide
                   && KeepOGlycopeptide == other.KeepOGlycopeptide;
        }

        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(MaxMissedCleavages);
            hash.Add(MinLength);
            hash.Add(MaxLength);
            hash.Add(MaxModificationIsoforms);
            hash.Add(MaxMods);
            hash.Add((int)InitiatorMethionineBehavior);
            hash.Add(Protease);
            hash.Add((int)SearchModeType);
            hash.Add((int)FragmentationTerminus);
            hash.Add(SpecificProtease);
            hash.Add(GeneratehUnlabeledProteinsForSilac);
            hash.Add(KeepNGlycopeptide);
            hash.Add(KeepOGlycopeptide);
            return hash.ToHashCode();
        }

        #endregion

        public override string ToString()
        {
            return MaxMissedCleavages + "," + InitiatorMethionineBehavior + "," + MinLength + "," + MaxLength + ","
                   + MaxModificationIsoforms + "," + MaxMods + "," + SpecificProtease.Name + "," + SearchModeType + "," + FragmentationTerminus + ","
                   + GeneratehUnlabeledProteinsForSilac + "," + KeepNGlycopeptide + "," + KeepOGlycopeptide;
        }

        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            var terminus = newTerminus ?? FragmentationTerminus;
            if (SearchModeType == CleavageSpecificity.None)
                return new DigestionParams(SpecificProtease.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                    GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide);
            return new DigestionParams(Protease.Name, MaxMissedCleavages, MinLength, MaxLength,
                MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide);
        }
            

        protected override DigestionAgent GetSingleTerminusAgent(FragmentationTerminus terminus)
        {
            return terminus == FragmentationTerminus.N
                ? ProteaseDictionary.Dictionary["singleN"]
                : ProteaseDictionary.Dictionary["singleC"];
        }
    }
}
