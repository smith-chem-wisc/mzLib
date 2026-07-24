using System;
using Omics.Digestion;
using Omics.Fragmentation;

namespace Proteomics.ProteolyticDigestion 
{
    public class DigestionParams : IDigestionParams, IEquatable<DigestionParams>
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public DigestionParams() : this("trypsin")
        {
        }

        public DigestionParams(string protease = "trypsin", int maxMissedCleavages = 2, int minPeptideLength = 7, int maxPeptideLength = int.MaxValue,
            int maxModificationIsoforms = 1024, InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
            int maxModsForPeptides = 2, CleavageSpecificity searchModeType = CleavageSpecificity.Full, FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            bool generateUnlabeledProteinsForSilac = true, bool keepNGlycopeptide = false, bool keepOGlycopeptide = false,
            bool respectCleavageBlockingModifications = false)
        {
            Protease = ProteaseDictionary.Dictionary[protease];
            MaxMissedCleavages = maxMissedCleavages;
            MinLength = minPeptideLength;
            MaxLength = maxPeptideLength;
            MaxMods = maxModsForPeptides;
            MaxModificationIsoforms = maxModificationIsoforms;
            InitiatorMethionineBehavior = initiatorMethionineBehavior;
            SearchModeType = searchModeType;
            FragmentationTerminus = fragmentationTerminus;
            RecordSpecificProtease();
            GeneratehUnlabeledProteinsForSilac = generateUnlabeledProteinsForSilac;
            KeepNGlycopeptide = keepNGlycopeptide;
            KeepOGlycopeptide = keepOGlycopeptide;
            RespectCleavageBlockingModifications = respectCleavageBlockingModifications;
        }

        /// <summary>
        /// Extra missed cleavages generated, beyond <see cref="MaxMissedCleavages"/>, when
        /// <see cref="RespectCleavageBlockingModifications"/> is on. A peptidoform whose C-terminal
        /// residue carries a cleavage-blocking modification is impossible, and its real counterpart
        /// reads THROUGH that residue to the next site -- which costs a missed cleavage under the
        /// ordinary (modification-blind) count. Without this slack that read-through form would never
        /// be generated at low limits, and at MaxMissedCleavages = 0 the real peptide would be lost
        /// entirely rather than merely mis-reported. The surplus is trimmed again by the open-site
        /// filter in ProteolyticPeptide.GetModifiedPeptides, so only genuinely-reachable peptidoforms
        /// survive; the cost is enumeration, not correctness.
        /// </summary>
        public const int CleavageBlockingReadThroughSlack = 2;

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MaxMissedCleavages { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxMods { get; set; }
        public DigestionAgent DigestionAgent => Protease;

        public CleavageSpecificity SearchModeType { get; private set; } //for fast semi and nonspecific searching of proteases
        public FragmentationTerminus FragmentationTerminus { get; private set; } //for fast semi searching of proteases
        public Protease SpecificProtease { get; private set; } //for fast semi and nonspecific searching of proteases
        public bool GeneratehUnlabeledProteinsForSilac { get; private set; } //used to look for unlabeled proteins (in addition to labeled proteins) for SILAC experiments
        public bool KeepNGlycopeptide { get; private set; }
        public bool KeepOGlycopeptide { get; private set; }

        /// <summary>
        /// When set, digestion treats a cleavage-blocking modification (see
        /// <see cref="Omics.Modifications.Modification.BlocksCleavage"/>) on a Lys/Arg as abolishing
        /// that cleavage site for the peptidoform carrying it. Peptidoforms whose C-terminus is such a
        /// residue are dropped -- trypsin could not have produced them -- and the blocked residue stops
        /// counting as a missed cleavage in the read-through form, so the real peptide survives even at
        /// MaxMissedCleavages = 0. Default false, which reproduces the historical (modification-blind)
        /// digestion exactly.
        /// </summary>
        public bool RespectCleavageBlockingModifications { get; private set; }

        /// <summary>
        /// The missed-cleavage budget used for GENERATION, as opposed to <see cref="MaxMissedCleavages"/>
        /// which remains the budget the caller asked for and which the open-site filter enforces.
        /// Only full-specificity digestion is inflated: the semi- and single-terminus modes do not
        /// enumerate by missed cleavage in the same way, and the C-terminal drop still applies to them.
        /// </summary>
        public int EffectiveMaxMissedCleavages =>
            RespectCleavageBlockingModifications && SearchModeType == CleavageSpecificity.Full
                ? MaxMissedCleavages + CleavageBlockingReadThroughSlack
                : MaxMissedCleavages;

        #region Properties overridden by more generic interface

        public Protease Protease { get; private set; }

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
                   && KeepOGlycopeptide == other.KeepOGlycopeptide
                   && RespectCleavageBlockingModifications == other.RespectCleavageBlockingModifications;
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
            hash.Add(RespectCleavageBlockingModifications);
            return hash.ToHashCode();
        }

        #endregion

        public override string ToString()
        {
            return MaxMissedCleavages + "," + InitiatorMethionineBehavior + "," + MinLength + "," + MaxLength + ","
                   + MaxModificationIsoforms + "," + MaxMods + "," + SpecificProtease.Name + "," + SearchModeType + "," + FragmentationTerminus + ","
                   + GeneratehUnlabeledProteinsForSilac + "," + KeepNGlycopeptide + "," + KeepOGlycopeptide + ","
                   + RespectCleavageBlockingModifications;
        }

        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            var terminus = newTerminus ?? FragmentationTerminus;
            if (SearchModeType == CleavageSpecificity.None)
                return new DigestionParams(SpecificProtease.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                    GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide,
                    RespectCleavageBlockingModifications);
            return new DigestionParams(Protease.Name, MaxMissedCleavages, MinLength, MaxLength,
                MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide);
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