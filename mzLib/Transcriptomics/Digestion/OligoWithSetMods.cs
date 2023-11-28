using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics;

namespace Transcriptomics
{
    public class OligoWithSetMods : NucleolyticOligo, IPrecursor, INucleicAcid
    {
        public OligoWithSetMods(NucleicAcid nucleicAcid, RnaDigestionParams digestionParams, int oneBaseStartResidue,
            int oneBasedEndResidue, int missedCleavages, CleavageSpecificity cleavageSpecificity,
            Dictionary<int, Modification> allModsOneIsNTerminus, int numFixedMods, IHasChemicalFormula? fivePrimeTerminus = null, 
            IHasChemicalFormula? threePrimeTerminus = null ) 
            : base(nucleicAcid, oneBaseStartResidue, oneBasedEndResidue, missedCleavages,
            cleavageSpecificity, fivePrimeTerminus, threePrimeTerminus)
        {
            _digestionParams = digestionParams;
            _allModsOneIsNterminus = allModsOneIsNTerminus;
            NumFixedMods = numFixedMods;
            FullSequence = (this as IPrecursor).DetermineFullSequence();
        }

        private RnaDigestionParams _digestionParams;
        private Dictionary<int, Modification> _allModsOneIsNterminus;
        private double? _monoisotopicMass;
        private ChemicalFormula? _thisChemicalFormula;
        private double? _mostAbundantMonoisotopicMass;
        private IDictionary<int, List<Modification>>? _oneBasedPossibleLocalizedModifications;

        public string FullSequence { get; private set; }
        public IDigestionParams DigestionParams => _digestionParams;
        public IBioPolymer Parent => NucleicAcid;
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
                if (_monoisotopicMass is null)
                {
                    _monoisotopicMass = BaseSequence.Sum(nuc => Nucleotide.GetResidue(nuc).MonoisotopicMass) +
                                        AllModsOneIsNterminus.Values.Sum(mod => mod.MonoisotopicMass.Value) +
                                        FivePrimeTerminus.MonoisotopicMass +
                                        ThreePrimeTerminus.MonoisotopicMass;
                }
                return _monoisotopicMass.Value;
            }
        }

        public ChemicalFormula ThisChemicalFormula
        {
            get
            {
                if (_thisChemicalFormula is null)
                {
                    var fullFormula = new RNA(BaseSequence, FivePrimeTerminus, ThreePrimeTerminus).GetChemicalFormula();
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
                }
                return _thisChemicalFormula!;
            }
        }

        public double MostAbundantMonoisotopicMass
        {
            get
            {
                if (_mostAbundantMonoisotopicMass is null)
                {
                    var distribution = IsotopicDistribution.GetDistribution(ThisChemicalFormula);
                    double maxIntensity = distribution.Intensities.Max();
                    _mostAbundantMonoisotopicMass =
                        distribution.Masses[distribution.Intensities.IndexOf(maxIntensity)].Round(9);
                }
                return _mostAbundantMonoisotopicMass.Value;
            }
        }

        public string SequenceWithChemicalFormulas => throw new NotImplementedException();

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
            List<IProduct> products)
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

        /// <summary>
        /// Generates theoretical internal fragments for given dissociation type for this peptide. 
        /// The "products" parameter is filled with these fragments.
        /// The "minLengthOfFragments" parameter is the minimum number of nucleic acids for an internal fragment to be included
        /// </summary>
        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<IProduct> products)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates all the fragments of the types you specify
        /// </summary>
        /// <param name="type">product type to get neutral fragments from</param>
        /// <param name="sequence">Sequence to generate fragments from, will be calculated from the parent if left null</param>
        /// <returns></returns>
        public IEnumerable<IProduct> GetNeutralFragments(ProductType type, Nucleotide[]? sequence = null)
        {
            sequence ??= (Parent as NucleicAcid)!.NucleicAcidArray[(OneBasedStartResidue - 1)..OneBasedEndResidue];

            if (type is ProductType.M)
            {
                yield return new RnaProduct(type, FragmentationTerminus.None, MonoisotopicMass, 0, 0, 0);
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

                yield return new RnaProduct(type,
                    isThreePrimeTerminal ? FragmentationTerminus.ThreePrime : FragmentationTerminus.FivePrime,
                    monoMass - neutralLoss, i,
                    isThreePrimeTerminal ? BaseSequence.Length - i : i, 0, null, 0);
            }
        }
    }
}
