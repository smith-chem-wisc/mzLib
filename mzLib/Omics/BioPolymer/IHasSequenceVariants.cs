using Omics.Modifications;
namespace Omics.BioPolymer;

/// <summary>
/// Represents a BioPolymer capable of having many sequence variants. 
/// Sequence variant parsing and construction is handled in <see cref="VariantApplication"/>
/// </summary>
public interface IHasSequenceVariants
{
    /// <summary>
    /// Primary Sequence
    /// </summary>
    string BaseSequence { get; }

    /// <summary>
    /// Sample name from which applied variants came, e.g. tumor or normal.
    /// </summary>
    public string SampleNameForVariants { get; }

    /// <summary>
    /// Modifications (values) located at one-based protein positions (keys)
    /// </summary>
    IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }

    /// <summary>
    /// Original modifications as defined in the Parsed XML database
    /// </summary>
    IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }

    /// <summary>
    /// Original BioPolymer used to construct all sequence variants. 
    /// </summary>
    IBioPolymer ConsensusVariant { get; }

    /// <summary>
    /// Sequence Variants from the database that are represented in this IHasSequenceVariants Base Sequence
    /// </summary>
    List<SequenceVariation> AppliedSequenceVariations { get; }

    /// <summary>
    /// Sequence Variants as defined in the parsed XML database
    /// </summary>
    public List<SequenceVariation> SequenceVariations { get; }

    /// <summary>
    /// Truncation products as defined in the parsed XML Database
    /// </summary>
    List<TruncationProduct> TruncationProducts { get; }

    
    /// <summary>
    /// Used to construct a new variant of the same type as the original and is called in <see cref="VariantApplication"/>
    /// </summary>
    /// <remarks>The generic structure enables proteins to produce proteins and RNA to produce RNA</remarks>
    /// <returns></returns>
    TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original, IEnumerable<SequenceVariation> appliedSequenceVariants,
        IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string sampleNameForVariants)
        where TBioPolymerType : IHasSequenceVariants;
}