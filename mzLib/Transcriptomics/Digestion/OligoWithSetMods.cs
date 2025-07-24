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

        public OligoWithSetMods(string sequence, Dictionary<string, Modification> allKnownMods, int numFixedMods = 0,
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

            FullSequence = sequence;
            _baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(sequence);
            _allModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(sequence, allKnownMods);
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
                        subsequence.Append('[' + mod.ChemicalFormula.Formula + ']');
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
            List<Product> products)
        {
            products.Clear();

            List<ProductType> fivePrimeProductTypes =
                dissociationType.GetRnaTerminusSpecificProductTypesFromDissociation(FragmentationTerminus.FivePrime);
            List<ProductType> threePrimeProductTypes =
                dissociationType.GetRnaTerminusSpecificProductTypesFromDissociation(FragmentationTerminus.ThreePrime);

            bool calculateFivePrime =
                fragmentationTerminus is FragmentationTerminus.FivePrime or FragmentationTerminus.Both;
            bool calculateThreePrime =
                fragmentationTerminus is FragmentationTerminus.ThreePrime or FragmentationTerminus.Both;

            var sequence = (Parent as NucleicAcid)!.NucleicAcidArray[(OneBasedStartResidue - 1)..OneBasedEndResidue];

            // intact product ion
            if (fragmentationTerminus is FragmentationTerminus.Both or FragmentationTerminus.None)
                products.AddRange(GetNeutralFragments(ProductType.M, sequence));

            if (calculateFivePrime)
                foreach (var type in fivePrimeProductTypes)
                    products.AddRange(GetNeutralFragments(type, sequence));

            if (calculateThreePrime)
                foreach (var type in threePrimeProductTypes)
                    products.AddRange(GetNeutralFragments(type, sequence));
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
            List<Product> products)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates all the fragments of the types you specify
        /// </summary>
        /// <param name="type">product type to get neutral fragments from</param>
        /// <param name="sequence">Sequence to generate fragments from, will be calculated from the parent if left null</param>
        /// <returns></returns>
        public IEnumerable<Product> GetNeutralFragments(ProductType type, Nucleotide[]? sequence = null)
        {
            sequence ??= (Parent as NucleicAcid)!.NucleicAcidArray[(OneBasedStartResidue - 1)..OneBasedEndResidue];

            if (type is ProductType.M)
            {
                yield return new Product(type, FragmentationTerminus.None, MonoisotopicMass, 0, 0, 0);
                yield break;
            }

            // determine mass of piece remaining after fragmentation
            double monoMass = type.GetRnaMassShiftFromProductType();

            // determine mass of terminal cap and add to fragment
            bool isThreePrimeTerminal = type.GetRnaTerminusType() == FragmentationTerminus.ThreePrime;
            IHasChemicalFormula terminus = isThreePrimeTerminal ? ThreePrimeTerminus : FivePrimeTerminus;
            monoMass += terminus.MonoisotopicMass;

            // determine mass of each polymer component that is contained within the fragment and add to fragment
            bool first = true; //set first to true to hand the terminus mod first
            for (int i = 0; i <= BaseSequence.Length - 1; i++)
            {
                int naIndex = isThreePrimeTerminal ? Length - i : i - 1;
                if (first)
                {
                    first = false; //set to false so only handled once
                    continue;
                }
                monoMass += sequence[naIndex].MonoisotopicMass;

                if (i < 1)
                    continue;

                // add side-chain mod
                if (AllModsOneIsNterminus.TryGetValue(naIndex + 2, out Modification mod))
                {
                    monoMass += mod.MonoisotopicMass ?? 0;
                }

                var previousNucleotide = sequence[naIndex];

                double neutralLoss = 0;
                if (type.ToString().Contains("Base"))
                {
                    neutralLoss = previousNucleotide.BaseChemicalFormula.MonoisotopicMass;
                }

                yield return new Product(type,
                    isThreePrimeTerminal ? FragmentationTerminus.ThreePrime : FragmentationTerminus.FivePrime,
                    monoMass - neutralLoss, i,
                    isThreePrimeTerminal ? BaseSequence.Length - i : i, 0, null, 0);
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
