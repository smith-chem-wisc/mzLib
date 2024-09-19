using Chemistry;
using Omics.Digestion;
using Omics.Modifications;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics.Digestion;

namespace Transcriptomics
{
    /// <summary>
    /// A linear polymer of Nucleic acids
    /// </summary>
    public abstract class NucleicAcid : INucleicAcid, IBioPolymer, IEquatable<NucleicAcid>
    {
        #region Static Properties

        /// <summary>
        /// The default chemical formula of the five prime (hydroxyl group)
        /// </summary>
        /// <remarks>
        /// This means that the five prime cap will remove the excess components of first nucleotides
        /// phospho group, leaving only the hydroxyl. This formula will be used for the five prime cap, unless
        /// the nucleic acid is constructed with a different chemical formula
        /// </remarks>
        public static readonly ChemicalFormula DefaultFivePrimeTerminus = ChemicalFormula.ParseFormula("O-3P-1");

        /// <summary>
        /// The default chemical formula of the three prime terminus (hydroxyl group)
        /// </summary>
        /// <remarks>
        /// This is used to account for the mass of the additional hydroxyl group at the three end of most oligonucleotides.
        /// This formula will be used for the three prime cap, unless the nucleic acid is constructed with a different
        /// chemical formula
        /// </remarks>
        public static readonly ChemicalFormula DefaultThreePrimeTerminus = ChemicalFormula.ParseFormula("OH");

        #endregion

        #region Constuctors

        /// <summary>
        /// For creating an RNA programatically
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="fivePrimeTerm"></param>
        /// <param name="threePrimeTerm"></param>
        /// <param name="oneBasedPossibleLocalizedModifications"></param>
        protected NucleicAcid(string sequence, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null)
        {
            MonoisotopicMass = 0;
            _nucleicAcids = new Nucleotide[sequence.Length];
            ThreePrimeTerminus = threePrimeTerm ??= DefaultThreePrimeTerminus;
            FivePrimeTerminus = fivePrimeTerm ??= DefaultFivePrimeTerminus;
            _oneBasedPossibleLocalizedModifications = oneBasedPossibleLocalizedModifications ?? new Dictionary<int, List<Modification>>();
            GeneNames = new List<Tuple<string, string>>();

            ParseSequence(sequence);
        }

        /// <summary>
        /// For Reading in from rna database
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="name"></param>
        /// <param name="identifier"></param>
        /// <param name="organism"></param>
        /// <param name="databaseFilePath"></param>
        /// <param name="fivePrimeTerm"></param>
        /// <param name="threePrimeTerm"></param>
        /// <param name="oneBasedPossibleLocalizedModifications"></param>
        /// <param name="isContaminant"></param>
        /// <param name="isDecoy"></param>
        /// <param name="additionalDatabaseFields"></param>
        protected NucleicAcid(string sequence, string name, string identifier, string organism, string databaseFilePath,
            IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null,
            bool isContaminant = false, bool isDecoy = false,
            Dictionary<string, string>? additionalDatabaseFields = null)
            : this(sequence, fivePrimeTerm, threePrimeTerm, oneBasedPossibleLocalizedModifications)
        {
            Name = name;
            DatabaseFilePath = databaseFilePath;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            Organism = organism;
            Accession = identifier;
            AdditionalDatabaseFields = additionalDatabaseFields;
        }

        #endregion

        #region Private Properties

        /// <summary>
        /// The 5-Prime chemical formula cap
        /// </summary>
        private IHasChemicalFormula _5PrimeTerminus;

        /// <summary>
        /// The 3-Prime chemical formula cap
        /// </summary>
        private IHasChemicalFormula _3PrimeTerminus;

        /// <summary>
        /// All of the nucleic acid residues indexed by position from 5- to 3-prime.
        /// </summary>
        private Nucleotide[] _nucleicAcids;

        /// <summary>
        /// The nucleic acid sequence. Is ignored if 'StoreSequenceString' is false
        /// </summary>
        private string _sequence;

        private IDictionary<int, List<Modification>> _oneBasedPossibleLocalizedModifications;

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
        public int Length => BaseSequence.Length;

        public string Name { get; }
        public string FullName => Name; // TODO: Consider if this needs to be different from the name
        public string DatabaseFilePath { get; }
        public bool IsDecoy { get; }
        public bool IsContaminant { get; }
        public string Accession { get; }

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications => _oneBasedPossibleLocalizedModifications;
        public string Organism { get; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public IEnumerable<Tuple<string, string>> GeneNames { get; }
        public Dictionary<string, string>? AdditionalDatabaseFields { get; }

        /// <summary>
        /// The total monoisotopic mass of this peptide and all of its modifications
        /// </summary>
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Returns a copy of the nucleic acid array, used for -base mass calculations.
        /// </summary>
        public Nucleotide[] NucleicAcidArray => _nucleicAcids;

        public ChemicalFormula ThisChemicalFormula => GetChemicalFormula();

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

        public char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        #endregion

        #region Digestion

        public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParameters, List<Modification> allKnownFixedMods,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
            bool topDownTruncationSearch = false)
        {
            if (digestionParameters is not RnaDigestionParams digestionParams)
                throw new ArgumentException(
                    "DigestionParameters must be of type DigestionParams for protein digestion");
            allKnownFixedMods ??= new();
            variableModifications ??= new();

            // digest based upon base sequence
            foreach (var unmodifiedOligo in digestionParams.Rnase.GetUnmodifiedOligos(this,
                         digestionParams.MaxMissedCleavages, digestionParams.MinLength, digestionParams.MaxLength))
            {
                // add fixed and variable mods to base sequence digestion products
                foreach (var modifiedOligo in unmodifiedOligo.GetModifiedOligos(allKnownFixedMods, digestionParams,
                             variableModifications))
                {
                    yield return modifiedOligo;
                }
            }
        }

        public IEnumerable<OligoWithSetMods> Digest(RnaDigestionParams digestionParameters,
            List<Modification> allKnownFixedMods,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null,
            (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
            bool topDownTruncationSearch = false)
        {
            return Digest((IDigestionParams)digestionParameters, allKnownFixedMods, variableModifications, silacLabels, turnoverLabels, topDownTruncationSearch)
                .Cast<OligoWithSetMods>();
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

            int index = 0;

            double monoMass = 0;
            ChemicalFormula chemFormula = new();

            StringBuilder sb = null;
            sb = new StringBuilder(sequence.Length);

            foreach (char letter in sequence)
            {
                Nucleotide residue;
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
                        case ' ': // ignore spaces
                            break;

                        case '*': // ignore *
                            break;

                        default:
                            throw new ArgumentException(string.Format(
                                "Nucleic Acid Letter {0} does not exist in the Nucleic Acid Dictionary. {0} is also not a valid character",
                                letter));
                    }
                }
            }

            _sequence = sb.ToString();
            MonoisotopicMass += monoMass;
            Array.Resize(ref _nucleicAcids, Length);

            return true;
        }

        #endregion

        #region Interface Implemntations and Overrides

        public bool Equals(NucleicAcid? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return _5PrimeTerminus.Equals(other._5PrimeTerminus)
                   && _3PrimeTerminus.Equals(other._3PrimeTerminus);
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
            return HashCode.Combine(_5PrimeTerminus, _3PrimeTerminus, _sequence);
        }

        #endregion
    }
}
