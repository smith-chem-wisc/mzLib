using Omics.BioPolymer;
using Omics.Modifications;
using System.Text;
using Transcriptomics.Digestion;

namespace Transcriptomics
{
    public static class ClassExtensions
    {
        /// <summary>
        /// Creates a new instance of a nucleic acid or oligo with set modifications, optionally updating its sequence, modifications, and decoy status.
        /// </summary>
        /// <typeparam name="T">The type of the nucleic acid, which must implement <see cref="INucleicAcid"/>.</typeparam>
        /// <param name="target">The target nucleic acid or oligo with set modifications to base the new instance on.</param>
        /// <param name="sequence">The new sequence string, if any. If null, the original sequence is used.</param>
        /// <param name="modifications">A dictionary of modifications to apply, if any. If null, the original modifications are used.</param>
        /// <param name="isDecoy">A flag indicating whether the sequence is a decoy, if any. If null, the original decoy status is used.</param>
        /// <returns>A new instance of the specified nucleic acid type with the provided or existing properties.</returns>
        /// <remarks>
        /// This method facilitates the generation of new sequences for both nucleic acids and oligos with set modifications by allowing
        /// optional updates to the sequence string, modifications, and decoy status. It ensures that the new instances are properly
        /// initialized with the provided or existing properties, enabling further analysis of modified sequences and future generation of decoys on the fly.
        /// </remarks>
        public static T CreateNew<T>(this T target, string? sequence = null, IDictionary<int, List<Modification>>? modifications = null,
        bool? isDecoy = null, List<TruncationProduct>? truncationProducts = null, List<SequenceVariation>? sequenceVariations = null,
        List<SequenceVariation>? appliedSequenceVariations = null, string decoyIdentifier = "DECOY")
            where T : INucleicAcid
        {
            // set new object parameters where not null
            object? returnObj = null;
            string newSequence = sequence ?? target.BaseSequence;
            IDictionary<int, List<Modification>> newModifications = modifications ?? target.OneBasedPossibleLocalizedModifications;
            
            switch (target)
            {
                case RNA rna:
                {
                    bool newIsDecoy = isDecoy ?? rna.IsDecoy;
                    string accession = newIsDecoy ? $"{decoyIdentifier}_{rna.Accession}" : rna.Accession;
                    List<TruncationProduct> newTruncs = truncationProducts ?? rna.TruncationProducts;
                    List<SequenceVariation> newVariations = sequenceVariations ?? rna.SequenceVariations;
                    List<SequenceVariation> newAppliedVariations = appliedSequenceVariations ?? rna.AppliedSequenceVariations;

                        returnObj = new RNA(newSequence, accession, newModifications, rna.FivePrimeTerminus,
                        rna.ThreePrimeTerminus, rna.Name, rna.Organism, rna.DatabaseFilePath, rna.IsContaminant,
                        newIsDecoy, rna.GeneNames, rna.AdditionalDatabaseFields, newTruncs,
                        newVariations, newAppliedVariations, rna.SampleNameForVariants, rna.FullName);
                    break;
                }
                case OligoWithSetMods oligo:
                {
                    var oldParent = oligo.Parent as RNA ?? throw new NullReferenceException();
                    bool newIsDecoy = isDecoy ?? oldParent.IsDecoy;
                    string accession = newIsDecoy ? $"{decoyIdentifier}_{oldParent.Accession}" : oldParent.Accession;
                    List<TruncationProduct> newTruncs = truncationProducts ?? oldParent.TruncationProducts;
                    List<SequenceVariation> newVariations = sequenceVariations ?? oldParent.SequenceVariations;
                    List<SequenceVariation> newAppliedVariations = appliedSequenceVariations ?? oldParent.AppliedSequenceVariations;

                    var newParent = new RNA(newSequence, accession, newModifications,oldParent.FivePrimeTerminus, oldParent.ThreePrimeTerminus, 
                    oldParent.Name, oldParent.Organism, oldParent.DatabaseFilePath, oldParent.IsContaminant, newIsDecoy, oldParent.GeneNames, oldParent.AdditionalDatabaseFields,
                    newTruncs, newVariations, newAppliedVariations, oldParent.SampleNameForVariants, oldParent.FullName);


                    returnObj = new OligoWithSetMods(
                        newParent,
                        (oligo.DigestionParams as RnaDigestionParams)!,
                        oligo.OneBasedStartResidue,
                        oligo.OneBasedEndResidue,
                        oligo.MissedCleavages,
                        oligo.CleavageSpecificityForFdrCategory,
                        newModifications.ToDictionary(p => p.Key, p => p.Value.First()),
                        oligo.NumFixedMods,
                        oligo.FivePrimeTerminus,
                        oligo.ThreePrimeTerminus);
                    break;
                }
                default:
                    throw new ArgumentException("INucleicAcid type not yet implemented");
            }

            return (T)returnObj ?? throw new NullReferenceException("Error creating new INucleicAcid");
        }

        /// <summary>
        /// Transcribes a DNA sequence into an RNA sequence
        /// </summary>
        /// <param name="dna">The input dna sequence</param>
        /// <param name="isCodingStrand">True if the input sequence is the coding strand, False if the input sequence is the template strand</param>
        /// <returns></returns>
        public static string Transcribe(this string dna, bool isCodingStrand = true)
        {
            var sb = new StringBuilder();
            foreach (var residue in dna)
            {
                if (isCodingStrand)
                {
                    sb.Append(residue == 'T' ? 'U' : residue);
                }
                else
                {
                    switch (residue)
                    {
                        case 'A':
                            sb.Append('U');
                            break;
                        case 'T':
                            sb.Append('A');
                            break;
                        case 'C':
                            sb.Append('G');
                            break;
                        case 'G':
                            sb.Append('C');
                            break;
                        default:
                            sb.Append(residue);
                            break;
                    }
                }
            }
            return sb.ToString();
        }
    }
}
