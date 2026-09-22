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
            bool generateUnlabeledProteinsForSilac = true, bool keepNGlycopeptide = false, bool keepOGlycopeptide = false)

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
        }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; private set; }
        public int MaxMissedCleavages { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxMods { get; set; }
        public DigestionAgent DigestionAgent => Protease;
        public DigestionAgent SpecificDigestionAgent => SpecificProtease;

        /// <summary>
        /// The kind of search: <see cref="CleavageSpecificity.Full"/> (the default), <see cref="CleavageSpecificity.Semi"/>
        /// or <see cref="CleavageSpecificity.None"/> (non-specific). Together with <see cref="FragmentationTerminus"/> it
        /// decides whether <c>Protein.Digest</c> returns peptides or seeds.
        /// </summary>
        /// <remarks>
        /// <para>For a fully specific protease such as trypsin:</para>
        /// <code>
        ///   SearchModeType  FragmentationTerminus  Protein.Digest returns
        ///   Full            any                    fully specific peptides
        ///   Semi            Both (the default)     semi-specific peptides (at least one terminus made by the protease)
        ///   Semi            N or C                 seeds fixed at that terminus, NOT peptides
        ///   None            N or C                 singleN or singleC seeds (non-specific), NOT peptides
        ///   None            Both                   singleC seeds, the same as None + C
        /// </code>
        /// <para><b>Peptides</b> can be scored as they are, so every search engine can use them. <b>Seeds</b> are long
        /// stretches fixed at one terminus whose other end is not decided yet; only an engine that decides it after the
        /// search, from the precursor mass, can use them. In MetaMorpheus that is the non-specific search engine alone, which
        /// clones these parameters to N and to C (<see cref="Clone"/>) and runs one pass for each. Classic, Modern, Glyco
        /// and crosslink searches must leave <see cref="FragmentationTerminus"/> at Both. There is no request for a list of
        /// non-specific peptides: None always gives seeds.</para>
        /// <para>A protease whose own <see cref="DigestionAgent.CleavageSpecificity"/> is Semi gives semi-specific peptides
        /// with SearchModeType Full. Test/ProteomicsTests/ProteolyticDigestion/SearchModeTypeDigestionTests.cs pins every
        /// row of the table, and SemiSpecificDigestionTests.cs checks the Semi rows peptide by peptide.</para>
        /// </remarks>
        public CleavageSpecificity SearchModeType { get; private set; }

        /// <summary>
        /// In digestion, whether a Semi or None <see cref="SearchModeType"/> asks for peptides (<c>Both</c>, the default)
        /// or for seeds fixed at the N- or C-terminus (<c>N</c>, <c>C</c>) for a search engine that trims them afterwards.
        /// It has no effect on a Full search. See the table on <see cref="SearchModeType"/>.
        /// </summary>
        public FragmentationTerminus FragmentationTerminus { get; private set; }

        /// <summary>
        /// The protease the caller named. It is the same object as <see cref="Protease"/> except when
        /// <see cref="SearchModeType"/> is None: then <see cref="Protease"/> becomes singleN (FragmentationTerminus N) or
        /// singleC (anything else) to make the non-specific seeds, and this keeps the named protease, whose sites still
        /// limit missed cleavages and are written to settings and output.
        /// </summary>
        public Protease SpecificProtease { get; private set; }
        public bool GeneratehUnlabeledProteinsForSilac { get; private set; } //used to look for unlabeled proteins (in addition to labeled proteins) for SILAC experiments
        public bool KeepNGlycopeptide { get; private set; }
        public bool KeepOGlycopeptide { get; private set; }

        #region Properties overridden by more generic interface

        /// <summary>
        /// The protease digestion runs: the one the caller named, or singleN/singleC when <see cref="SearchModeType"/> is
        /// None (see <see cref="SpecificProtease"/>).
        /// </summary>
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
