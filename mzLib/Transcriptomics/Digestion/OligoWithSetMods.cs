using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics;
using Easy.Common.Extensions;
using Omics.Fragmentation.Oligo;
using System.Text;

namespace Transcriptomics.Digestion
{
    /// <summary>
    /// Represents an oligonucleotide with set modifications, providing properties and methods for
    /// accessing and manipulating its chemical characteristics.
    /// </summary>
    /// <remarks>
    /// The monoisotopic mass, most abundant mass, and chemical formula are calculated on the fly if the corresponding properties
    /// (_monoisotopicMass, _thisChemicalFormula, _mostAbundantMonoisotopicMass) are null. This ensures that the most up-to-date values are
    /// always available based on the current state of the oligonucleotide and its modifications. Therefor, it is important to set those
    /// properties to null whenever a termini or modification is changed.
    /// </remarks>
    public class OligoWithSetMods : NucleolyticOligo, IBioPolymerWithSetMods, INucleicAcid, IEquatable<OligoWithSetMods>
    {
        public OligoWithSetMods(NucleicAcid nucleicAcid, RnaDigestionParams digestionParams, int oneBaseStartResidue,
            int oneBasedEndResidue, int missedCleavages, CleavageSpecificity cleavageSpecificity,
            Dictionary<int, Modification> allModsOneIsNTerminus, int numFixedMods, IHasChemicalFormula? fivePrimeTerminus = null,
            IHasChemicalFormula? threePrimeTerminus = null)
            : base(nucleicAcid, oneBaseStartResidue, oneBasedEndResidue, missedCleavages,
            cleavageSpecificity, fivePrimeTerminus, threePrimeTerminus)
        {
            _digestionParams = digestionParams;
            _allModsOneIsNterminus = allModsOneIsNTerminus;
            NumFixedMods = numFixedMods;
            FullSequence = this.DetermineFullSequence();
        }

        public OligoWithSetMods(string sequence, Dictionary<string, Modification>? allKnownMods = null, int numFixedMods = 0,
            RnaDigestionParams digestionParams = null, NucleicAcid n = null, int oneBaseStartResidue = 1, int oneBasedEndResidue = 0,
             int missedCleavages = 0, CleavageSpecificity cleavageSpecificity = CleavageSpecificity.Full, string description = null,
            IHasChemicalFormula? fivePrimeTerminus = null, IHasChemicalFormula? threePrimeTerminus = null)
            : base(n, oneBaseStartResidue, oneBasedEndResidue, missedCleavages,
                cleavageSpecificity, fivePrimeTerminus, threePrimeTerminus)
        {
            if (sequence.Contains("|"))
            {
                throw new MzLibUtil.MzLibException("Ambiguous oligo cannot be parsed from string: " + sequence);
            }

            _baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(sequence);
            _allModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(sequence, allKnownMods ?? Mods.AllKnownRnaModsDictionary);
            FullSequence = _allModsOneIsNterminus.ContainsKey(_baseSequence.Length + 2) 
                ? this.DetermineFullSequence() 
                : sequence;
            NumFixedMods = numFixedMods;
            _digestionParams = digestionParams;
            Description = description;

            if (n != null)
                Parent = n;
        }

        private RnaDigestionParams _digestionParams;
        private Dictionary<int, Modification> _allModsOneIsNterminus;
        private double? _monoisotopicMass;
        private ChemicalFormula? _thisChemicalFormula;
        private double? _mostAbundantMonoisotopicMass;
        private IDictionary<int, List<Modification>>? _oneBasedPossibleLocalizedModifications;
        private string? _sequenceWithChemicalFormula;

        public string FullSequence { get; private set; }
        public IDigestionParams DigestionParams => _digestionParams;
        public IHasChemicalFormula FivePrimeTerminus
        {
            get => _fivePrimeTerminus;
            set
            {
                _fivePrimeTerminus = value;
                _monoisotopicMass = null;
                _thisChemicalFormula = null;
                _mostAbundantMonoisotopicMass = null;
            }
        }

        public IHasChemicalFormula ThreePrimeTerminus
        {
            get => _threePrimeTerminus;
            set
            {
                _threePrimeTerminus = value;
                _monoisotopicMass = null;
                _thisChemicalFormula = null;
                _mostAbundantMonoisotopicMass = null;
            }
        }

        public double MonoisotopicMass
        {
            get
            {
                _monoisotopicMass ??= BaseSequence.Sum(nuc => Nucleotide.GetResidue(nuc).MonoisotopicMass) +
                                      AllModsOneIsNterminus.Values.Sum(mod => mod.MonoisotopicMass!.Value) +
                                      FivePrimeTerminus.MonoisotopicMass +
                                      ThreePrimeTerminus.MonoisotopicMass;
                return _monoisotopicMass.Value;
            }
        }

        public ChemicalFormula ThisChemicalFormula
        {
            get
            {
                if (_thisChemicalFormula is not null) return _thisChemicalFormula!;

                var fullFormula = new RNA(BaseSequence, fivePrimeTerm: FivePrimeTerminus, threePrimeTerm: ThreePrimeTerminus).GetChemicalFormula();
                foreach (var mod in AllModsOneIsNterminus.Values)
                {
                    if (mod.ChemicalFormula is null)
                    {
                        fullFormula = null;
                        break;
                    }
                    fullFormula.Add(mod.ChemicalFormula);
                }
                _thisChemicalFormula = fullFormula;
                return _thisChemicalFormula!;
            }
        }

        public double MostAbundantMonoisotopicMass
        {
            get
            {
                if (_mostAbundantMonoisotopicMass is not null) return _mostAbundantMonoisotopicMass.Value;

                var distribution = IsotopicDistribution.GetDistribution(ThisChemicalFormula);
                double maxIntensity = distribution.Intensities.Max();
                _mostAbundantMonoisotopicMass = distribution.Masses[distribution.Intensities.IndexOf(maxIntensity)].RoundedDouble();
                return _mostAbundantMonoisotopicMass!.Value;
            }
        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (_sequenceWithChemicalFormula is not null) return _sequenceWithChemicalFormula;

                var subsequence = new StringBuilder();
                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification? pepNTermVariableMod))
                {
                    if (pepNTermVariableMod is { } mod)
                        subsequence.Append('[' + mod.ChemicalFormula.Formula + ']');
                }

                for (int r = 0; r < Length; r++)
                {
                    subsequence.Append(this[r]);
                    // variable modification on this residue
                    if (!AllModsOneIsNterminus.TryGetValue(r + 2, out Modification? residueVariableMod)) continue;
                    if (residueVariableMod is { } mod)
                        subsequence.Append('[' + mod.ChemicalFormula.Formula + ']');
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification? pepCTermVariableMod))
                {
                    if (pepCTermVariableMod is { } mod)
                        subsequence.Append("-[" + mod.ChemicalFormula.Formula + ']');
                }

                _sequenceWithChemicalFormula = subsequence.ToString();
                return _sequenceWithChemicalFormula;
            }
        }

        public Dictionary<int, Modification> AllModsOneIsNterminus => _allModsOneIsNterminus;

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications => _oneBasedPossibleLocalizedModifications ??=
            _allModsOneIsNterminus.ToDictionary(p => p.Key, p => new List<Modification>() { p.Value });
        public int NumMods => AllModsOneIsNterminus.Count;
        public int NumFixedMods { get; }
        public int NumVariableMods => NumMods - NumFixedMods;

        /// <summary>
        /// Generates theoretical fragments for given dissociation type for this peptide. 
        /// The "products" parameter is filled with these fragments.
        /// </summary>
        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
            List<Product> products, IFragmentationParams? fragmentationParams = null)
        {
            products.Clear();
            fragmentationParams ??= RnaFragmentationParams.Default;
            bool modsCanSuppressBaseLossIons = fragmentationParams is RnaFragmentationParams
            {
                ModificationsCanSuppressBaseLossIons: true
            };

            List<ProductType> fivePrimeProductTypes =
                dissociationType.GetRnaTerminusSpecificProductTypesFromDissociation(FragmentationTerminus.FivePrime);
            List<ProductType> threePrimeProductTypes =
                dissociationType.GetRnaTerminusSpecificProductTypesFromDissociation(FragmentationTerminus.ThreePrime);

            bool calculateFivePrime =
                fragmentationTerminus is FragmentationTerminus.FivePrime or FragmentationTerminus.Both;
            bool calculateThreePrime =
                fragmentationTerminus is FragmentationTerminus.ThreePrime or FragmentationTerminus.Both;

            Nucleotide[] sequence;
            if (NucleicAcid is null) // If no parent, construct the nucleotide array ourselves
            {
                sequence = BaseSequence.Select(p => Nucleotide.TryGetResidue(p, out Nucleotide? nuc) ? nuc : throw new MzLibUtil.MzLibException($"Invalid nucleotide '{p}' in sequence.")).ToArray();
            }
            else
            {
                sequence = NucleicAcid.NucleicAcidArray[(OneBasedStartResidue - 1)..OneBasedEndResidue];
            }

            // intact product ion
            if (fragmentationParams.GenerateMIon && fragmentationTerminus is FragmentationTerminus.Both or FragmentationTerminus.None)
                products.AddRange(this.GetMIons(fragmentationParams));

            if (calculateFivePrime)
                foreach (var type in fivePrimeProductTypes)
                    products.AddRange(GetNeutralFragments(type, sequence, modsCanSuppressBaseLossIons));

            if (calculateThreePrime)
                foreach (var type in threePrimeProductTypes)
                    products.AddRange(GetNeutralFragments(type, sequence, modsCanSuppressBaseLossIons));
        }

        #region IEquatable

        /// <summary>
        /// Oligos are equal if they have the same full sequence, parent, and digestion agent, and terminal caps
        /// </summary>
        public override bool Equals(object? obj)
        {
            if (obj is OligoWithSetMods oligo)
            {
                return Equals(oligo);
            }
            return false;
        }

        /// <summary>
        /// Oligos are equal if they have the same full sequence, parent, and digestion agent, and terminal caps
        /// </summary>
        public bool Equals(IBioPolymerWithSetMods? other) => Equals(other as OligoWithSetMods);

        /// <summary>
        /// Oligos are equal if they have the same full sequence, parent, and digestion agent, and terminal caps
        /// </summary>
        public bool Equals(OligoWithSetMods? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (other.GetType() != GetType()) return false;

            // for those constructed from sequence and mods only
            if (Parent is null && other.Parent is null)
                return FullSequence.Equals(other.FullSequence);

            return FullSequence == other.FullSequence
                   && Equals(DigestionParams?.DigestionAgent, other.DigestionParams?.DigestionAgent)
                   && _fivePrimeTerminus.Equals(other._fivePrimeTerminus)
                   && _threePrimeTerminus.Equals(other._threePrimeTerminus)
                   // These last two are important for parsimony in MetaMorpheus
                   && OneBasedStartResidue == other!.OneBasedStartResidue
                   && Equals(Parent?.Accession, other.Parent?.Accession);
        }

        public override int GetHashCode()
        {
            var hash = new HashCode();
            hash.Add(FullSequence);
            hash.Add(OneBasedStartResidue);
            if (Parent?.Accession != null)
            {
                hash.Add(Parent.Accession);
            }
            if (DigestionParams?.DigestionAgent != null)
            {
                hash.Add(DigestionParams.DigestionAgent);
            }
            hash.Add(FivePrimeTerminus);
            hash.Add(ThreePrimeTerminus);
            return hash.ToHashCode();
        }

        #endregion

        /// <summary>
        /// Generates theoretical internal fragments for given dissociation type for this peptide. 
        /// The "products" parameter is filled with these fragments.
        /// The "minLengthOfFragments" parameter is the minimum number of nucleic acids for an internal fragment to be included
        /// </summary>
        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<Product> products, IFragmentationParams? fragmentationParams = null)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates all fragments of the specified type that can be generated from this oligonucleotide, including any modifications that should be included in those fragments.
        /// </summary>
        /// <param name="type">product type to get neutral fragments from</param>
        /// <param name="sequence">Sequence to generate fragments from, will be calculated from the parent if left null</param>
        /// <returns></returns>
        public IEnumerable<Product> GetNeutralFragments(ProductType type, Nucleotide[]? sequence = null, bool modsCanSuppressBaseLossIons = false)
        {
            sequence ??= (Parent as NucleicAcid)!.NucleicAcidArray[(OneBasedStartResidue - 1)..OneBasedEndResidue];

            // determine mass of piece remaining after fragmentation
            double monoMass = type.GetRnaMassShiftFromProductType();

            // determine mass of terminal cap and add to fragment
            bool isThreePrimeTerminal = type.GetRnaTerminusType() == FragmentationTerminus.ThreePrime;
            IHasChemicalFormula terminus = isThreePrimeTerminal ? ThreePrimeTerminus : FivePrimeTerminus;
            monoMass += terminus.MonoisotopicMass;

            // determine mass of each polymer component that is contained within the fragment and add to fragment iteratively, starting from the terminus and moving inward until the end of the oligo. 
            bool first = true; //set first to true to hand the terminus sideChainMod first
            for (int fragmentNumber = 0; fragmentNumber <= BaseSequence.Length - 1; fragmentNumber++)
            {
                // 0-based array index of nucleotide being added
                int nucleicAcidIndex = isThreePrimeTerminal ? BaseSequence.Length - fragmentNumber : fragmentNumber - 1;

                // 1-based sequence position (fragment boundary)
                int residuePosition = isThreePrimeTerminal ? BaseSequence.Length - fragmentNumber : fragmentNumber;

                // Mod at side chain being added. 2-based index for modifications on the current residue in the AllModsOneIsNterminus dictionary. 
                int sideChainModIndex = nucleicAcidIndex + 2;

                // Mod at the phosphate linkage after (5') or before (3') the current residue being added.
                int phosphateModIndex = residuePosition + 1;

                //For a2(5' fragment containing A-U):
                //    •	fragmentNumber = 2
                //    •	nucleicAcidIndex = 1 → Points to U(0 - based)
                //    •	residuePosition = 2 → Position 2 in sequence(1 - based)
                //    •	Side - chain mod lookup: AllModsOneIsNterminus[3] → Mod on U
                //    •	Backbone mod lookup: AllModsOneIsNterminus[3] → Phosphate between U - G
                //For w2(3' fragment containing G-C):
                //    •	fragmentNumber = 2
                //    •	nucleicAcidIndex = 2 → Points to G(0 - based)
                //    •	residuePosition = 2 → Position 2 from 3' end
                //    •	Side - chain mod lookup: AllModsOneIsNterminus[4] → Mod on G
                //    •	Backbone mod lookup: AllModsOneIsNterminus[3] → Phosphate between U - G

                if (first)
                {
                    // Add 3'-term mod if present
                    if (isThreePrimeTerminal && AllModsOneIsNterminus.TryGetValue(BaseSequence.Length + 2, out Modification? threePrimeMod))
                    {
                            monoMass += threePrimeMod.MonoisotopicMass ?? 0;

                    }
                    // Add 5'-term mod if present
                    else if (!isThreePrimeTerminal && AllModsOneIsNterminus.TryGetValue(1, out Modification? fivePrimeMod))
                    {
                        monoMass += fivePrimeMod.MonoisotopicMass ?? 0;
                    }

                    first = false; //set to false so only handled once
                    continue;
                }
                monoMass += sequence[nucleicAcidIndex].MonoisotopicMass;

                if (fragmentNumber < 1)
                    continue;

                // add side-chain sideChainMod only (at current position)
                if (AllModsOneIsNterminus.TryGetValue(sideChainModIndex, out Modification? sideChainMod) && sideChainMod is not BackboneModification)
                {
                    monoMass += sideChainMod.MonoisotopicMass ?? 0;
                }

                // Add backbone modifications if they are included in the fragment, otherwise add mod mass after this fragment. 
                double? backboneMassShift = null;
                if (AllModsOneIsNterminus.TryGetValue(phosphateModIndex, out Modification? mod) && mod is BackboneModification bm)
                {
                    if (Array.BinarySearch(bm.ProductsContainingModMass, type) >= 0)
                        monoMass += mod.MonoisotopicMass ?? 0;
                    else
                        backboneMassShift = mod.MonoisotopicMass;
                }

                // Handle Base Loss fragment series mass correction. 
                double neutralLoss = 0;
                var previousNucleotide = sequence[nucleicAcidIndex];
                if (type.IsBaseLoss())
                {
                    neutralLoss = previousNucleotide.BaseChemicalFormula.MonoisotopicMass;
                    var generateBaseLossIon = true;

                    if (sideChainMod is BaseModification baseMod)
                    {
                        switch (baseMod.BaseLossType)
                        {
                            case BaseLossBehavior.Suppressed when modsCanSuppressBaseLossIons:
                                generateBaseLossIon = false; // Don't generate base-loss ion
                                break;

                            case BaseLossBehavior.Suppressed:
                            case BaseLossBehavior.Modified:
                                // Add modification mass to base loss
                                if (baseMod.BaseLossModification != null)
                                {
                                    neutralLoss += baseMod.BaseLossModification?.MonoisotopicMass ?? 0;
                                }
                                break;

                            case BaseLossBehavior.Default:
                            default:
                                // Normal base loss
                                break;
                        }
                    }

                    if (!generateBaseLossIon)
                        continue;
                }

                yield return new Product(type,
                    isThreePrimeTerminal ? FragmentationTerminus.ThreePrime : FragmentationTerminus.FivePrime,
                    monoMass - neutralLoss, fragmentNumber,
                    residuePosition, 0, null, 0);

                if (backboneMassShift != null)
                    monoMass += backboneMassShift.Value; // add the backbone mass shift back for the next iteration if it was not added to this fragment 
            }
        }

        /// <summary>
        /// Outputs a duplicate IBioPolymerWithSetMods with a localized mass shift, replacing a modification when present
        /// <remarks>
        /// Used to localize an unknown mass shift in the MetaMorpheus Localization Engine
        /// </remarks>
        /// </summary>
        /// <param name="indexOfMass">The index of the modification in the AllModOneIsNTerminus Dictionary - 2 (idk why -2)</param>
        /// <param name="massToLocalize">The mass to add to the BioPolymer</param>
        public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize)
        {
            var dictWithLocalizedMass = new Dictionary<int, Modification>(AllModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (dictWithLocalizedMass.TryGetValue(indexOfMass + 2, out Modification modToReplace))
            {
                massOfExistingMod = (double)modToReplace.MonoisotopicMass;
                dictWithLocalizedMass.Remove(indexOfMass + 2);
            }

            dictWithLocalizedMass.Add(indexOfMass + 2, new Modification(_locationRestriction: "Anywhere.", _monoisotopicMass: massToLocalize + massOfExistingMod));

            var peptideWithLocalizedMass = new OligoWithSetMods(NucleicAcid, _digestionParams, OneBasedStartResidue, OneBasedEndResidue, MissedCleavages,
                CleavageSpecificityForFdrCategory, dictWithLocalizedMass, NumFixedMods, FivePrimeTerminus, ThreePrimeTerminus);

            return peptideWithLocalizedMass;
        }
    }
}
