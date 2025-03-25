using Omics.Modifications;

namespace Omics.BioPolymer;

public interface IHasSequenceVariants
{
    string BaseSequence { get; }
    public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
    
    IBioPolymer NonVariant { get; }
    List<SequenceVariation> AppliedSequenceVariations { get; }
    IEnumerable<TruncationProduct> ProteolysisProducts { get; }



}