using System.Collections.Generic;

namespace Proteomics
{
    public class ProteinWithAppliedVariants
    {
        public string VariantBaseSequence { get; }
        public Protein Protein { get; }
        public List<SequenceVariation> AppliedSequenceVariations { get; }
        public int Length { get { return Protein.Length; } }
        public string Individual { get; }

        public ProteinWithAppliedVariants(string variantBaseSequence, Protein protein, List<SequenceVariation> appliedSequenceVariations, string individual)
        {
            VariantBaseSequence = variantBaseSequence;
            Protein = protein;
            AppliedSequenceVariations = appliedSequenceVariations ?? new List<SequenceVariation>();
            Individual = individual;
        }

        public char this[int zeroBasedIndex]
        {
            get
            {
                return VariantBaseSequence[zeroBasedIndex];
            }
        }
    }
}