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
            bool respectCleavageBlockingModifications = false,
            bool respectCleavagePromotingModifications = false)
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
            RespectCleavagePromotingModifications = respectCleavagePromotingModifications;
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

        /// <summary>
        /// When set, digestion honours a protease whose motif REQUIRES a modification at one of its
        /// subsites -- a glycoprotease. A peptidoform whose cut is not justified by the required
        /// glycan is dropped, because the protease could not have made that cut. Default false, which
        /// reproduces the historical (modification-blind) digestion exactly.
        /// </summary>
        /// <remarks>
        /// <para><b>PRECONDITION: the glycan must be somewhere digestion can see it.</b> The requirement
        /// is evaluated against the modifications available at digestion time -- those the database
        /// annotates on the protein, and those the search configures as fixed or variable. Either source
        /// will do, and they are both consulted. But if NEITHER carries a modification that can satisfy
        /// the requirement, then no cut anywhere is justified and the protein comes back essentially
        /// undigested.</para>
        ///
        /// <para>That is the correct answer, not a defect: StcE, OpeRATOR and IMPa demonstrably do not
        /// cleave unglycosylated substrate, and returning the intact substrate is what the published
        /// unglycosylated controls describe. It does mean the flag must not be set by a workflow that
        /// resolves the glycan AFTER identification rather than placing it at digestion -- a glyco search,
        /// where the glycan enters as a precursor-mass offset and is localized against fragment ions, has
        /// no glycan in the digestion's view at all and would simply lose its peptides. For that workflow
        /// the requirement has to constrain localization instead, which this flag does not do.</para>
        ///
        /// <para>The mirror of <see cref="RespectCleavageBlockingModifications"/>, and gated the same way:
        /// full-specificity searches only. The two corrections move in opposite directions. Blocking
        /// removes a site the bare sequence had, so it needs generation SLACK to reach the read-through
        /// peptide that replaces what it drops. Promoting works in two stages instead. The first removes
        /// sites at which the required modification could never be present, before enumeration, and needs
        /// no slack: the read-through across a site that is not a site is just the ordinary peptide
        /// between the sites that remain. The second refines OCCUPANCY per peptidoform, and it does NOT
        /// yet have a read-through -- a peptidoform whose cut is feasible but unoccupied is dropped with
        /// nothing generated to replace it (truth-set case STCE-08). Note that the two corrections are
        /// therefore NOT symmetric in one further respect: an unconfigured blocking modification cannot
        /// remove a site the sequence really has, whereas an unsatisfiable promoting requirement means
        /// there is genuinely no site -- which is why this flag has no "nothing configured, go inert"
        /// escape and must not be given one.</para>
        ///
        /// Scope: like its sibling, this applies to full-specificity searches only. A semi or
        /// nonspecific search is left exactly as it was. There the peptide's termini are not all
        /// protease cuts, and dropping a peptidoform whose length-driven terminus happens to sit at an
        /// unglycosylated motif would remove a peptide the search should still see.
        /// </remarks>
        public bool RespectCleavagePromotingModifications { get; private set; }

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
                   && KeepOGlycopeptide == other.KeepOGlycopeptide
                   && RespectCleavageBlockingModifications == other.RespectCleavageBlockingModifications
                   && RespectCleavagePromotingModifications == other.RespectCleavagePromotingModifications;
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
            hash.Add(RespectCleavagePromotingModifications);
            return hash.ToHashCode();
        }

        #endregion

        public override string ToString()
        {
            return MaxMissedCleavages + "," + InitiatorMethionineBehavior + "," + MinLength + "," + MaxLength + ","
                   + MaxModificationIsoforms + "," + MaxMods + "," + SpecificProtease.Name + "," + SearchModeType + "," + FragmentationTerminus + ","
                   + GeneratehUnlabeledProteinsForSilac + "," + KeepNGlycopeptide + "," + KeepOGlycopeptide + ","
                   + RespectCleavageBlockingModifications + ","
                   + RespectCleavagePromotingModifications;
        }

        public IDigestionParams Clone(FragmentationTerminus? newTerminus = null)
        {
            var terminus = newTerminus ?? FragmentationTerminus;
            if (SearchModeType == CleavageSpecificity.None)
                return new DigestionParams(SpecificProtease.Name, MaxMissedCleavages, MinLength, MaxLength,
                    MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                    GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide,
                    RespectCleavageBlockingModifications, RespectCleavagePromotingModifications);
            return new DigestionParams(Protease.Name, MaxMissedCleavages, MinLength, MaxLength,
                MaxModificationIsoforms, InitiatorMethionineBehavior, MaxMods, SearchModeType, terminus,
                GeneratehUnlabeledProteinsForSilac, KeepNGlycopeptide, KeepOGlycopeptide,
                RespectCleavageBlockingModifications, RespectCleavagePromotingModifications);
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
