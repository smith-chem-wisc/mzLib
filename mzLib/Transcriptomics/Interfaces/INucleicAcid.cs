using Chemistry;
using Omics;
using Omics.Modifications;

namespace Transcriptomics
{
    public interface INucleicAcid : IHasChemicalFormula, IBioPolymer
    {
        /// <summary>
        /// The amino acid sequence
        /// </summary>
        string BaseSequence { get; }

        /// <summary>
        /// The length of the amino acid sequence
        /// </summary>
        int Length { get; }

        /// <summary>
        /// Modifications 
        /// </summary>
        IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }


        IHasChemicalFormula FivePrimeTerminus { get; set; }

        IHasChemicalFormula ThreePrimeTerminus { get; set; }
    }
}
