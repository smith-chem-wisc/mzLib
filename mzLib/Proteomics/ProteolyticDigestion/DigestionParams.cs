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
        /// <remarks>
        /// <see cref="MaxMissedCleavages"/> keeps its meaning and its guarantee: no peptide leaves
        /// digestion reporting more missed cleavages than were asked for. A blocked residue is not a
        /// cleavage site for the peptidoform carrying it, so it is not a missed cleavage either -- the
        /// count reports cleavages that could have happened and did not, not Lys/Arg residues. Digestion
        /// does enumerate a wider span internally to reach the read-through forms, but that slack is
        /// discounted again before a peptide is emitted, so it is never observable in the result.
        ///
        /// Scope: this flag applies to full-specificity SEARCHES only
        /// (<see cref="SearchModeType"/> == <see cref="CleavageSpecificity.Full"/>). A semi or
        /// nonspecific search is left exactly as it was, deliberately and entirely. The drop and the
        /// wider generation span are two halves of one exchange -- the impossible peptidoform leaves and
        /// the read-through form that replaces it arrives -- and applying only the first half would make
        /// a peptide unidentifiable rather than correctly identified whenever the budget is too small to
        /// reach the read-through. Half a correction is worse than none, so semi gets none.
        ///
        /// Making semi searches benefit properly is follow-up work, and it is not just a matter of
        /// widening this gate: a semi peptide's C-terminus may be a genuine protease cut or a
        /// length-driven truncation, and telling them apart needs the protease's site list, which full
        /// digestion gets free from its own enumeration and semi digestion does not.
        /// </remarks>
        public bool RespectCleavageBlockingModifications { get; private set; }

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
                GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide,
                RespectCleavageBlockingModifications);
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