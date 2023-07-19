using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MassSpectrometry;
using static System.Net.Mime.MediaTypeNames;

namespace Transcriptomics
{
    /// <summary>
    /// A linear polymer of Nucleic acids
    /// </summary>
    public abstract class NucleicAcid : INucleicAcid, IEquatable<NucleicAcid>, IHasMass
    {

        #region Static Properties

        /// <summary>
        /// The default chemical formula of the five prime (hydroxyl group)
        /// </summary>
        public static readonly ChemicalFormula DefaultFivePrimeTerminus = ChemicalFormula.ParseFormula("O-3P-1");

        /// <summary>
        /// The default chemical formula of the three prime terminus (hydroxyl group)
        /// </summary>
        public static readonly ChemicalFormula DefaultThreePrimeTerminus = ChemicalFormula.ParseFormula("OH");

        #endregion

        #region Constuctors

        //protected NucleicAcid(IHasMass[] modifications)
        //{
        //    Modifications = modifications;
        //}

        protected NucleicAcid(string sequence) 
            : this(sequence, DefaultFivePrimeTerminus, DefaultThreePrimeTerminus)
        {
        }

        protected NucleicAcid(string sequence, IHasChemicalFormula fivePrimeTerm, IHasChemicalFormula threePrimeTerm)
        {
            MonoisotopicMass = 0;
            Length = sequence.Length;
            _nucleicAcids = new Nucleotide[Length];
            ThreePrimeTerminus = threePrimeTerm;
            FivePrimeTerminus = fivePrimeTerm;
            ParseSequence(sequence);
        }

        #endregion

        #region Private Properties

        /// <summary>
        /// The 5-Prime chemical formula cap. This is different from the 5-prime terminus modification.
        /// </summary>
        private IHasChemicalFormula _5PrimeTerminus;

        /// <summary>
        /// The 3-Prime chemical formula cap. This is different from the 3-prime terminus modification.
        /// </summary>
        private IHasChemicalFormula _3PrimeTerminus;

        /// <summary>
        /// All of the modifications indexed by position from 5- to 3-prime termini. This array is 2 bigger than the nucleic acid array
        /// as index 0 and Count - 1 represent the 5- and 3-prime terminus, respectively
        /// </summary>
        private IHasMass[] _modifications;

        /// <summary>
        /// All of the nucleic acid residues indexed by position from 5- to 3-prime.
        /// </summary>
        private Nucleotide[] _nucleicAcids;

        /// <summary>
        /// The nucleic acid sequence with modification names interspersed. Is ignored if 'StoreSequenceString' is false
        /// </summary>
        private string _sequenceWithMods;

        /// <summary>
        /// The nucleic acid sequence. Is ignored if 'StoreSequenceString' is false
        /// </summary>
        private string _sequence;

        /// <summary>
        /// The internal flag to represent that the sequence with modifications have been changed and need to be updated
        /// </summary>
        internal bool IsDirty { get; set; }

        #endregion

        #region Internal Properties

        /// <summary>
        /// The internal data store for the modifications (2 larger than the length to handle the 5' and 3' termini)
        /// </summary>
        internal IHasMass[] Modifications { get; }

        /// <summary>
        /// The internal data store for the nucleic acids
        /// </summary>
        internal Nucleotide[] NucleicAcids => _nucleicAcids;

        #endregion

        #region Public Properties

        /// <summary>
        /// Gets or sets the 5' terminus of this nucleic acid polymer
        /// </summary>
        public IHasChemicalFormula FivePrimeTerminus
        {
            get => _5PrimeTerminus;
            set => ReplaceTerminus(ref _5PrimeTerminus, value);
        }

        /// <summary>
        /// Gets or sets the 3' terminus of this nucleic acid polymer
        /// </summary>
        public IHasChemicalFormula ThreePrimeTerminus
        {
            get => _3PrimeTerminus;
            set => ReplaceTerminus(ref _3PrimeTerminus, value);
        }

        /// <summary>
        /// Gets the number of nucleic acids in this nucleic acid polymer
        /// </summary>
        public int Length { get; private set; }

        /// <summary>
        /// The total monoisotopic mass of this peptide and all of its modifications
        /// </summary>
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Returns a copy of the nucleic acid array, used for -base mass calculations.
        /// </summary>
        public Nucleotide[] NucleicAcidArray => _nucleicAcids;

        #endregion

        #region Nucleic Acid Sequence

        /// <summary>
        /// Gets the base nucleic acid sequence
        /// </summary>
        public string BaseSequence
        {
            get
            {
                // Generate the sequence if the stored version is null or empty
                if (string.IsNullOrEmpty(_sequence))
                {
                    _sequence = new string(_nucleicAcids.Select(na => na.Letter).ToArray());
                }

                return _sequence;
            }
        }

        ///// <summary>
        ///// Gets the nucleic acid sequence with modifications
        ///// </summary>
        //public string FullSequence
        //{
        //    get
        //    {
        //        if (!IsDirty && !string.IsNullOrEmpty(_sequenceWithMods))
        //            return _sequenceWithMods;

        //        _sequenceWithMods = GetSequenceWithModifications();
        //        IsDirty = false;
        //        return _sequenceWithMods;
        //    }
        //}

        #endregion

        #region Fragmentation

        

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

            if (calculateFivePrime)
            {
                foreach (var type in fivePrimeProductTypes)
                {
                    products.AddRange(GetNeutralFragments(type));
                }
            }

            if (calculateThreePrime)
            {
                foreach (var type in threePrimeProductTypes)
                {
                    products.AddRange(GetNeutralFragments(type));
                }
            }
        }

        /// <summary>
        /// Calculates all the fragments of the types you specify
        /// </summary>
        /// <param name="type"></param>
        /// <returns></returns>
        internal IEnumerable<IProduct> GetNeutralFragments(ProductType type)
        {

            // determine mass of piece remaining after fragmentation
            double monoMass = type.GetRnaMassShiftFromProductType();

            // determine mass of terminal cap and add to fragment
            bool isThreePrimeTerminal = type.GetRnaTerminusType() == FragmentationTerminus.ThreePrime;
            IHasChemicalFormula terminus = isThreePrimeTerminal ? ThreePrimeTerminus : FivePrimeTerminus;
            monoMass += terminus.MonoisotopicMass;



            // determine mass of each polymer component that is contained within the fragment and add to fragment
            bool first = true; //set first to true to hand the terminus mod first
            bool hasMod = _modifications != null;

            for (int i = 0; i <= BaseSequence.Length - 1; i++)
            {

                int naIndex = isThreePrimeTerminal ? Length - i : i - 1;

                // Handle the terminus mods first in a special case
                IHasMass mod;
                if (first)
                {
                    first = false; //set to false so only handled once
                    if (hasMod) //checks for if there are modifications in the array (if there are mods on the RNA)
                    {
                        mod = _modifications[naIndex + 1];
                        if (mod != null)
                        {
                            monoMass += mod.MonoisotopicMass;
                        }
                    }
                    continue;
                }

                monoMass += _nucleicAcids[naIndex].MonoisotopicMass;

                if (hasMod)
                {
                    mod = _modifications[naIndex + 1];
                    if (mod != null)
                    {
                        monoMass += mod.MonoisotopicMass;
                    }
                }

                if (i < 1)
                    continue;

                var previousNucleotide = _nucleicAcids[naIndex];

                double neutralLoss = 0;
                if (type.ToString().Contains("Base"))
                {
                    neutralLoss = previousNucleotide.BaseChemicalFormula.MonoisotopicMass;
                }

                yield return new RnaProduct(type,
                    isThreePrimeTerminal ? FragmentationTerminus.ThreePrime : FragmentationTerminus.FivePrime,
                    monoMass - neutralLoss, i,
                    i, 0, null, 0);
            }
        }

        #endregion

        #region Modifications

        // TODO:

        #endregion

        #region Digestion

        public IEnumerable<NucleicAcid> Digest()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region Electrospray

        public IEnumerable<double> GetElectrospraySeries(int minCharge, int maxCharge)
        {
            for (int i = minCharge; i < maxCharge; i++)
            {
                yield return this.ToMz(i);
            }
        }

        

        #endregion

        #region Chemical Formula

        public ChemicalFormula GetChemicalFormula()
        {
            var formula = new ChemicalFormula();

            // Handle Modifications
            //if (ContainsModifications())
            //{
            //    for (int i = 0; i < Length + 2; i++)
            //    {
            //        IChemicalFormula chemMod = _modifications[i] as IChemicalFormula;

            //        if (chemMod == null)
            //            continue;

            //        formula.Add(chemMod.ChemicalFormula);
            //    }
            //}

            // Handle 5'-Terminus
            formula.Add(FivePrimeTerminus.ThisChemicalFormula);

            // Handle 3'-Terminus
            formula.Add(ThreePrimeTerminus.ThisChemicalFormula);

            // Handle Nucleic Acid Residues
            for (int i = 0; i < Length; i++)
            {
                formula.Add(_nucleicAcids[i].ThisChemicalFormula);
            }

            return formula;
        }

        #endregion

        #region Private Methods

        private bool ReplaceTerminus(ref IHasChemicalFormula terminus, IHasChemicalFormula value)
        {
            if (Equals(value, terminus))
                return false;

            if (terminus != null)
                MonoisotopicMass -= terminus.MonoisotopicMass;

            terminus = value;

            if (value != null)
                MonoisotopicMass += value.MonoisotopicMass;

            return true;
        }

        /// <summary>
        /// Parses a string sequence of nucleic acids characters into a peptide object
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        private bool ParseSequence(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return false;

            bool inMod = false;
            bool threePrimeTerminalMod = false; // 5' or 3' terminal modification
            int index = 0;

            double monoMass = 0;
            ChemicalFormula chemFormula = new();

            StringBuilder sb = null;
            sb = new StringBuilder(sequence.Length);

            StringBuilder modSb = new StringBuilder(10);
            foreach (char letter in sequence)
            {
                if (inMod)
                {
                    if (letter == ']')
                    {
                        inMod = false; // end the modification phase

                        string modString = modSb.ToString();
                        modSb.Clear();
                        IHasMass modification = null;
                        switch (modString)
                        {
                            default:
                                double mass;
                                //Modification mod;
                                ////if (ModificationDictionary.TryGetModification(modString, out mod))
                                ////{
                                ////    modification = mod;
                                ////}
                                ////else if (ChemicalFormula.IsValidChemicalFormula(modString))
                                ////{
                                ////    modification = new ChemicalFormula(modString);
                                ////}
                                ////else if (double.TryParse(modString, out mass))
                                ////{
                                ////    modification = new Mass(mass);
                                ////}
                                ////else
                                ////{
                                ////    throw new ArgumentException("Unable to correctly parse the following modification: " + modString);
                                ////}
                                break;
                        }

                        monoMass += modification.MonoisotopicMass;

                        if (_modifications == null)
                            _modifications = new IHasMass[Length + 2];

                        if (threePrimeTerminalMod)
                        {
                            _modifications[index + 1] = modification;
                        }
                        else
                        {
                            _modifications[index] = modification;
                        }

                        threePrimeTerminalMod = false;

                        //TODO: Account for mods in chemical formula
                    }
                    else
                    {
                        modSb.Append(letter);
                    }
                }
                else
                {
                    Nucleotide residue;
                    //char upperletter = char.ToUpper(letter); // moved to nucleic acid dictionary
                    if (Nucleotide.TryGetResidue(letter, out residue))
                    {
                        _nucleicAcids[index++] = residue;
                        sb.Append(residue.Letter);
                        monoMass += residue.MonoisotopicMass;
                    }
                    else
                    {
                        switch (letter)
                        {
                            case '[': // start of a modification
                                inMod = true;
                                break;

                            case '-': // start of a 3'-terminal modification
                                threePrimeTerminalMod = (index > 0);
                                break;

                            case ' ': // ignore spaces
                                break;

                            case '*': // ignore *
                                break;

                            default:
                                throw new ArgumentException(string.Format("Nucleic Acid Letter {0} does not exist in the Nucleic Acid Dictionary. {0} is also not a valid character", letter));
                        }
                    }
                }
            }

            if (inMod)
            {
                throw new ArgumentException("Couldn't find the closing ] for a modification in this sequence: " + sequence);
            }


            _sequence = sb.ToString();
            Length = index;
            MonoisotopicMass += monoMass;
            Array.Resize(ref _nucleicAcids, Length);
            if (_modifications != null)
                Array.Resize(ref _modifications, Length + 2);

            IsDirty = true;

            return true;
        }

        #endregion

        #region Interface Implemntations and Overrides

        public bool Equals(NucleicAcid? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _5PrimeTerminus.Equals(other._5PrimeTerminus)
                   && _3PrimeTerminus.Equals(other._3PrimeTerminus)
                   && _sequenceWithMods == other._sequenceWithMods;
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((NucleicAcid)obj);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(_5PrimeTerminus, _3PrimeTerminus, _sequenceWithMods);
        }

        #endregion
    }
}
