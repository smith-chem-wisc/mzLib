using Omics.Digestion;
using Omics.Modifications;

namespace Omics
{
    public interface IBioPolymer 
    {
        string Name { get; }
        string FullName { get; }
        string BaseSequence { get; }
        int Length { get; }
        string DatabaseFilePath { get; }
        bool IsDecoy { get; }
        bool IsContaminant { get; }
        string Organism { get; }
        string Accession { get; }
        /// <summary>
        /// This tuple originates from reading UniProt XML protein entries. It consists of two parts.
        /// The second part is the gene name, which is also often referred to as the gene symbol.
        /// The first part is a descriptor for the second part. For example, "primary" is a word used for the
        /// first part to indicate that the gene symbol/name that follows is the name most widely accepted
        /// by the community. The word "synonym" in the first part would indicate that the gene name in the
        /// second part is an alternative. The word "ORF", which stands for open reading frame in the first part
        /// is also acceptable. Here, the second part will include a range, for example "AO055_05240"
        /// </summary>
        IEnumerable<Tuple<string, string>> GeneNames { get; }
        IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
        char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false);
    }
}
