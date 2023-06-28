using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using static System.Net.Mime.MediaTypeNames;

namespace Transcriptomics
{
    /// <summary>
    /// A linear polymer of Nucleic acids
    /// </summary>
    public class NucleicAcid : INucleicAcid, IHasMass, IEquatable<NucleicAcid>
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

        protected NucleicAcid(IHasMass[] modifications)
        {
            Modifications = modifications;
        }

        protected NucleicAcid(string sequence, IHasChemicalFormula fivePrimeTerm, IHasChemicalFormula threePrimeTerm, IHasMass[] modifications)
        {
            Modifications = modifications;
            MonoisotopicMass = 0;
            Length = sequence.Length;
            ThreePrimeTerminus = threePrimeTerm;
            FivePrimeTerminus = fivePrimeTerm;
        }

        #endregion

        #region Private Properties *DONE*

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

        #region Internal Properties *DONE*

        /// <summary>
        /// The internal data store for the modifications (2 larger than the length to handle the 5' and 3' termini)
        /// </summary>
        internal IHasMass[] Modifications { get; }

        /// <summary>
        /// The internal data store for the nucleic acids
        /// </summary>
        internal Nucleotide[] NucleicAcids
        {
            get { return _nucleicAcids; }
        }

        #endregion

        #region Public Properties *DONE*

        /// <summary>
        /// Gets or sets the 5' terminus of this nucleic acid polymer
        /// </summary>
        public IHasChemicalFormula FivePrimeTerminus
        {
            get { return _5PrimeTerminus; }
            set { ReplaceTerminus(ref _5PrimeTerminus, value); }
        }

        /// <summary>
        /// Gets or sets the 3' terminus of this nucleic acid polymer
        /// </summary>
        public IHasChemicalFormula ThreePrimeTerminus
        {
            get { return _3PrimeTerminus; }
            set { ReplaceTerminus(ref _3PrimeTerminus, value); }
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
        public Nucleotide[] NucleicAcidArray
        {
            get { return _nucleicAcids; }
        }

        #endregion

        #region Nucleic Acid Sequence

        /// <summary>
        /// Gets the base nucleic acid sequence
        /// </summary>
        public string Sequence
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
        //public string SequenceWithModifications
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

        // TODO:

        #endregion

        #region Modifications

        // TODO:

        #endregion

        #region Chemical Formula

        // TODO:

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
                                Modification mod;
                                //if (ModificationDictionary.TryGetModification(modString, out mod))
                                //{
                                //    modification = mod;
                                //}
                                //else if (ChemicalFormula.IsValidChemicalFormula(modString))
                                //{
                                //    modification = new ChemicalFormula(modString);
                                //}
                                //else if (double.TryParse(modString, out mass))
                                //{
                                //    modification = new Mass(mass);
                                //}
                                //else
                                //{
                                //    throw new ArgumentException("Unable to correctly parse the following modification: " + modString);
                                //}
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
            throw new NotImplementedException();
        }

        #endregion


     
        
    }
}
