using Chemistry;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Transcriptomics
{
    public class RNA : NucleicAcid
    {
        /// <summary>
        /// For constructing RNA from a string
        /// </summary>
        public RNA(string sequence, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null)
            : base(sequence, fivePrimeTerm, threePrimeTerm, oneBasedPossibleLocalizedModifications)
        {
        }

        /// <summary>
        /// For use with RNA loaded from a database
        /// </summary>
        public RNA(string sequence, string name, string accession, string organism, string databaseFilePath,
            IHasChemicalFormula? fivePrimeTerminus = null, IHasChemicalFormula? threePrimeTerminus = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleModifications = null,
            bool isContaminant = false, bool isDecoy = false, List<Tuple<string, string>> geneNames = null,
            Dictionary<string, string>? databaseAdditionalFields = null, List<TruncationProduct>? truncationProducts = null,
            List<SequenceVariation>? sequenceVariations = null, List<SequenceVariation>? appliedSequenceVariations = null,
            string? sampleNameForVariants = null, string? fullName = null)
            : base(sequence, name, accession, organism, databaseFilePath, fivePrimeTerminus, threePrimeTerminus,
                oneBasedPossibleModifications, isContaminant, isDecoy, geneNames, databaseAdditionalFields, truncationProducts,
                sequenceVariations, appliedSequenceVariations, sampleNameForVariants, fullName)
        {
        }
        
        /// <summary>
        /// For creating a variant of an existing nucleic acid. Filters out modifications that do not match their nucleotide target site.
        /// </summary>
        public RNA(string variantBaseSequence, NucleicAcid original, IEnumerable<SequenceVariation> appliedSequenceVariants,
            IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)

            : this(variantBaseSequence, original.Name, VariantApplication.GetAccession(original, appliedSequenceVariants), original.Organism, original.DatabaseFilePath,
                original.FivePrimeTerminus, original.ThreePrimeTerminus, oneBasedModifications, original.IsContaminant, original.IsDecoy, original.GeneNames.ToList(),
                original.AdditionalDatabaseFields, applicableProteolysisProducts.ToList(), original.SequenceVariations.ToList(), appliedSequenceVariants.ToList(), sampleNameForVariants, original.FullName)
        {
            NonVariant = original.NonVariant;
        }

        public override TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
            IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
        {
            var variantRNA = new RNA(variantBaseSequence, original as RNA, appliedSequenceVariants, applicableProteolysisProducts, oneBasedModifications, sampleNameForVariants);
            return (TBioPolymerType)(IHasSequenceVariants)variantRNA;
        }
    }
}
