using Chemistry;

namespace Transcriptomics
{
    public interface INucleotide : IHasChemicalFormula
    {
        char Letter { get; }
        string Symbol { get; }
    }
}
