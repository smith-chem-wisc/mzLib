using Omics.Modifications;

namespace Omics.BioPolymer;

public interface IHasSequenceVariants
{
    string BaseSequence { get; }
    IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
    IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
    
    IBioPolymer NonVariant { get; }
    List<SequenceVariation> AppliedSequenceVariations { get; }
    IEnumerable<TruncationProduct> ProteolysisProducts { get; }


    IHasSequenceVariants CreateVariant(string variantBaseSequence, IHasSequenceVariants original, IEnumerable<SequenceVariation> appliedSequenceVariants,
        IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants);
}