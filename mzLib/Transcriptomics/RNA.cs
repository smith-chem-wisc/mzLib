using Chemistry;
using Easy.Common.Extensions;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Transcriptomics
{
    public class RNA : NucleicAcid, IEquatable<RNA>
    {
        /// <summary>
        /// For constructing RNA from a string
        /// </summary>
        public RNA(string sequence,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null,
            IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null)
            : base(sequence, oneBasedPossibleLocalizedModifications, fivePrimeTerm, threePrimeTerm)
        {
        }

        /// <summary>
        /// For use with RNA loaded from a database
        /// </summary>
        public RNA(string sequence, string accession,
            IDictionary<int, List<Modification>>? oneBasedPossibleModifications = null,
            IHasChemicalFormula? fivePrimeTerminus = null, IHasChemicalFormula? threePrimeTerminus = null,
            string? name = null, string? organism = null, string? databaseFilePath = null,
            bool isContaminant = false, bool isDecoy = false, List<Tuple<string, string>> geneNames = null,
            Dictionary<string, string>? databaseAdditionalFields = null,
            List<TruncationProduct>? truncationProducts = null,
            List<SequenceVariation>? sequenceVariations = null,
            List<SequenceVariation>? appliedSequenceVariations = null,
            string? sampleNameForVariants = null, string? fullName = null)
            : base(sequence, accession, oneBasedPossibleModifications, fivePrimeTerminus, threePrimeTerminus,
                name, organism, databaseFilePath, isContaminant, isDecoy, geneNames, databaseAdditionalFields,
                truncationProducts, sequenceVariations, appliedSequenceVariations, sampleNameForVariants, fullName)
        {
        }
        
        /// <summary>
        /// For creating a variant of an existing nucleic acid. Filters out modifications that do not match their nucleotide target site.
        /// </summary>
        public RNA(string variantBaseSequence, NucleicAcid original, IEnumerable<SequenceVariation>? appliedSequenceVariants,
            IEnumerable<TruncationProduct>? applicableTruncationProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
            : this(variantBaseSequence, VariantApplication.GetAccession(original, appliedSequenceVariants), 
                  oneBasedModifications,
                  original.FivePrimeTerminus, original.ThreePrimeTerminus,
                  VariantApplication.GetVariantName(original.Name, appliedSequenceVariants), 
                  original.Organism, original.DatabaseFilePath, original.IsContaminant, 
                  original.IsDecoy, original.GeneNames, original.AdditionalDatabaseFields,
                  [..applicableTruncationProducts ?? new List<TruncationProduct>()], original.SequenceVariations, 
                  [..appliedSequenceVariants ?? new List<SequenceVariation>()], sampleNameForVariants, 
                  VariantApplication.GetVariantName(original.FullName, appliedSequenceVariants))
        {
            ConsensusVariant = original.ConsensusVariant;
            OriginalNonVariantModifications = ConsensusVariant.OriginalNonVariantModifications;
            AppliedSequenceVariations = (appliedSequenceVariants ?? new List<SequenceVariation>()).ToList();
            SampleNameForVariants = sampleNameForVariants;
        }

        public override IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods = null)
        {
            // Create a new rna with the new base sequence and modifications
            RNA newRna = this.CreateNew(newBaseSequence, newMods);
            return newRna;
        }

        public override TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
            IEnumerable<TruncationProduct> applicableTruncationProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
        {
            if (original is not RNA rna)
                throw new ArgumentException("The original nucleic acid must be RNA to create an RNA variant.");

            var variantRNA = new RNA(variantBaseSequence, rna, appliedSequenceVariants, 
                applicableTruncationProducts, oneBasedModifications, sampleNameForVariants);
            return (TBioPolymerType)(IHasSequenceVariants)variantRNA;
        }

        public bool Equals(RNA? other)
        {
            // interface equals first because it does null and reference checks
            return (this as NucleicAcid).Equals(other);
        }
    }
}
