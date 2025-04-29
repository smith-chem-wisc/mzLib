using Chemistry;
using Omics.Digestion;
using Omics.Modifications;
using Omics;
using System.Text;
using MzLibUtil;
using Transcriptomics.Digestion;
using Omics.BioPolymer;

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
        /// For creating an RNA programatically. Filters out modifications that do not match their nucleotide target site.
        /// </summary>
        protected NucleicAcid(string sequence,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null,
            IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null)
        {
            ConsensusVariant = this;
            MonoisotopicMass = 0;
            _nucleicAcids = new Nucleotide[sequence.Length];
            ThreePrimeTerminus = threePrimeTerm ?? DefaultThreePrimeTerminus;
            FivePrimeTerminus = fivePrimeTerm ?? DefaultFivePrimeTerminus;
            ParseSequenceString(sequence);

            SequenceVariations = [];
            AppliedSequenceVariations = [];
            TruncationProducts = [];
            OriginalNonVariantModifications = oneBasedPossibleLocalizedModifications ?? new Dictionary<int, List<Modification>>();
            OneBasedPossibleLocalizedModifications = oneBasedPossibleLocalizedModifications != null 
                ? ((IBioPolymer)this).SelectValidOneBaseMods(oneBasedPossibleLocalizedModifications) 
                : new Dictionary<int, List<Modification>>();
        }

        /// <summary>
        /// For Reading in from rna database. Filters out modifications that do not match their nucleotide target site.
        /// </summary>
        protected NucleicAcid(string sequence, string accession,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null,
            IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            string? name = null, string? organism = null,
            string? databaseFilePath = null,
            bool isContaminant = false, bool isDecoy = false, List<Tuple<string, string>>? geneNames = null,
            Dictionary<string, string>? additionalDatabaseFields = null,
            List<TruncationProduct>? truncationProducts = null,
            List<SequenceVariation>? sequenceVariations = null,
            List<SequenceVariation>? appliedSequenceVariations = null,
            string? sampleNameForVariants = null, string? fullName = null)
            : this(sequence, oneBasedPossibleLocalizedModifications, fivePrimeTerm, threePrimeTerm)
        {
            Name = name ?? "";
            DatabaseFilePath = databaseFilePath ?? "";
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            Organism = organism ?? "";
            Accession = accession;
            AdditionalDatabaseFields = additionalDatabaseFields;

            GeneNames = geneNames ?? [];
            TruncationProducts = truncationProducts ?? [];
            SequenceVariations = sequenceVariations ?? [];
            AppliedSequenceVariations = appliedSequenceVariations ?? [];
            SampleNameForVariants = sampleNameForVariants ?? "";
            FullName = fullName ?? name;
        }

        public abstract IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods = null);

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
        public string FullName { get; } 
        public string DatabaseFilePath { get; }
        public bool IsDecoy { get; }
        public bool IsContaminant { get; }
        public string Accession { get; }

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; protected set; }

        public string Organism { get; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public List<Tuple<string, string>> GeneNames { get; }
        public Dictionary<string, string>? AdditionalDatabaseFields { get; }

        /// <summary>
        /// The total monoisotopic mass of this nucleic acid and all of its modifications
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
                throw new MzLibException(
                    "DigestionParameters must be of type DigestionParams for protein digestion", new ArgumentException());
            allKnownFixedMods ??= new();
            variableModifications ??= new();

            // digest based upon base sequence
            foreach (var unmodifiedOligo in digestionParams.Rnase.GetUnmodifiedOligos(this,
                         digestionParams.MaxMissedCleavages, digestionParams.MinLength, digestionParams.MaxLength))
            {
                // add fixed and variable mods to base sequence digestion products
                foreach (var modifiedOligo in unmodifiedOligo.GenerateModifiedOligos(allKnownFixedMods, digestionParams,
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

        #region Sequence Variants

        public IBioPolymer ConsensusVariant { get; init; }

        /// <summary>
        /// Sequence Variants as defined in the parsed XML database
        /// </summary>
        public List<SequenceVariation> SequenceVariations { get; protected set; }

        /// <summary>
        /// Truncation products as defined in the parsed XML Database
        /// </summary>
        public List<TruncationProduct> TruncationProducts { get; protected set; }

        /// <summary>
        /// Sequence variations that have been applied to the base sequence.
        /// </summary>
        public List<SequenceVariation> AppliedSequenceVariations { get; protected set; }

        /// <summary>
        /// Sample name from which applied variants came, e.g. tumor or normal.
        /// </summary>
        public string SampleNameForVariants { get; protected set; }

        /// <summary>
        /// Original modifications as defined in the Parsed XML database
        /// </summary>
        public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

        // Abstract so we can do this construction in the appropriate derived class
        public abstract TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
            IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
            where TBioPolymerType : IHasSequenceVariants;

        #endregion

        #region Electrospray

        public IEnumerable<double> GetElectrospraySeries(int minCharge, int maxCharge)
        {
            if (minCharge > maxCharge)
                (minCharge, maxCharge) = (maxCharge, minCharge);
            
            for (int i = maxCharge; i > minCharge - 1; i--)
                yield return this.ToMz(i); 
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

        private void ReplaceTerminus(ref IHasChemicalFormula? terminus, IHasChemicalFormula? value)
        {
            if (Equals(value, terminus))
                return;

            if (terminus != null)
                MonoisotopicMass -= terminus.MonoisotopicMass;

            terminus = value;

            if (value != null)
                MonoisotopicMass += value.MonoisotopicMass;
        }

        /// <summary>
        /// Parses a string sequence of nucleic acid characters into an array of Nucleotide objects,
        /// updates the sequence string, and calculates the monoisotopic mass.
        /// </summary>
        /// <param name="sequence">The string sequence of nucleic acid characters to parse.</param>
        private void ParseSequenceString(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return;

            int index = 0;
            double monoMass = 0;

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
        }

        #endregion

        #region Interface Implemntations and Overrides

        public bool Equals(NucleicAcid? other)
        {
            // interface equals first because it does null and reference checks
            return (this as IBioPolymer).Equals(other)
                   && _5PrimeTerminus.Equals(other._5PrimeTerminus)
                   && _3PrimeTerminus.Equals(other._3PrimeTerminus);
        }

        public override bool Equals(object? obj)
        {
            if (obj is NucleicAcid oligo)
            {
                return Equals(oligo);
            }
            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(_5PrimeTerminus, _3PrimeTerminus, _sequence);
        }

        #endregion
    }
}
