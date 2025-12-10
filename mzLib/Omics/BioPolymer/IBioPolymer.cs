using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;

namespace Omics
{
    public interface IBioPolymer : IEquatable<IBioPolymer>, IHasSequenceVariants
    {
        string Name { get; }
        string FullName { get; }
        new string BaseSequence { get; } // new keyword is to define inheritance from IHasSequenceVariants to IBioPolymer
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
        List<Tuple<string, string>> GeneNames { get; }
        new IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; } // new keyword is to define inheritance from IHasSequenceVariants to IBioPolymer
        char this[int zeroBasedIndex] => BaseSequence[zeroBasedIndex];

        IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications,
            List<Modification> variableModifications, List<SilacLabel>? silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false);

        IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>>? newMods);

        bool IEquatable<IBioPolymer>.Equals(IBioPolymer? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (other.GetType() != GetType()) return false;
            return Accession == other.Accession
                && BaseSequence == other.BaseSequence;
        }

        /// <summary>
        /// Filters modifications that do not match their target amino acid.
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        IDictionary<int, List<Modification>> SelectValidOneBaseMods(IDictionary<int, List<Modification>> dict)
        {
            Dictionary<int, List<Modification>> validModDictionary = new Dictionary<int, List<Modification>>();
            foreach (KeyValuePair<int, List<Modification>> entry in dict)
            {
                List<Modification> validMods = new List<Modification>();
                foreach (Modification m in entry.Value)
                {
                    //mod must be valid mod and the motif of the mod must be present in the protein at the specified location
                    if (m.ValidModification && ModificationLocalization.ModFits(m, BaseSequence, 0, BaseSequence.Length, entry.Key))
                    {
                        validMods.Add(m);
                    }
                }

                if (validMods.Any())
                {
                    if (validModDictionary.Keys.Contains(entry.Key))
                    {
                        validModDictionary[entry.Key].AddRange(validMods);
                    }
                    else
                    {
                        validModDictionary.Add(entry.Key, validMods);
                    }
                }
            }
            return validModDictionary;
        }
    }
}
