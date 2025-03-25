using Omics.Modifications;
namespace Omics.BioPolymer;

/// <summary>
/// Represents a BioPolymer capable of having many sequence variants. 
/// Sequence variant parsing and construction is handled in <see cref="VariantApplication"/>
/// </summary>
public interface IHasSequenceVariants
{
    string BaseSequence { get; }
    IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
    IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
    IBioPolymer NonVariant { get; }
    List<SequenceVariation> AppliedSequenceVariations { get; }
    IEnumerable<TruncationProduct> TruncationProducts { get; }

    /// <summary>
    /// Used to construct a new variant of the same type as the original and is called in <see cref="VariantApplication"/>
    /// </summary>
    /// <remarks>The generic structure enables proteins to produce proteins and RNA to produce RNA</remarks>
    /// <returns></returns>
    TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
        IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
        where TBioPolymerType : IHasSequenceVariants;
}