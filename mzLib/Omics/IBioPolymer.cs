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
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the nameof the specific molecule.
        /// There may be many genes and names of a certain type produced when reading an XML protein database.
        /// For RNA there may be a gene name and a specific transcript name due to splice variants.
        /// </summary>
        IEnumerable<Tuple<string, string>> GeneNames { get; }
        IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }
        char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false);
    }
}
