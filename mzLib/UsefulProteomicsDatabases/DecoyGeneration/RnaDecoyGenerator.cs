using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Omics.Modifications;
using Transcriptomics;
using Omics.BioPolymer;
using Omics;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Provides methods for generating decoy nucleic acids from any implementor of <see cref="INucleicAcid"/>.
    /// </summary>
    /// <remarks>
    /// This class supports various types of decoy generation, including reversing, sliding, and shuffling sequences.
    /// It allows for the creation of decoy sequences while preserving certain characteristics such as modification sites and termini.
    /// The <c>GenerateDecoys</c> method serves as the main entry point, delegating to specific decoy generation methods based on the specified <see cref="DecoyType"/>.
    /// TODO: Implement Shuffle and Slide Decoys
    /// TODO: Consider passing digestion motif as optional parameter to leave digestion sites intact. Currently leaving the 3' intact as it is the predominant cleavage motif.
    /// TODO: Consider palindromic sequences and the result they have on fragment ions (d/z are identical, c/y are identical). This will be particularly important for slided decoys
    /// </remarks>
    public static class RnaDecoyGenerator
    {
        public static List<T> GenerateDecoys<T>(List<T> nucleicAcids, DecoyType decoyType, int maxThreads = -1, string decoyIdentifier = "DECOY") where T : INucleicAcid
        {
            switch (decoyType)
            {
                case DecoyType.None:
                    return new List<T>();
                case DecoyType.Reverse:
                    return GenerateReverseDecoys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Slide:
                    return GenerateSlidedDecoys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Shuffle:
                    return GenerateShuffledDeocys(nucleicAcids, maxThreads, decoyIdentifier);
                case DecoyType.Random:
                default:
                    throw new ArgumentOutOfRangeException(nameof(decoyType), decoyType, null);
            }
        }

        /// <summary>
        /// Generated decoys in which the sequence is reversed,
        /// leaving modification on their nucleic acid of origin,
        /// and 3' termini intact as it is the most likely cleavage site. 
        /// </summary>
        /// <param name="nucleicAcids"></param>
        /// <param name="maxThreads"></param>
        /// <returns></returns>
        private static List<T> GenerateReverseDecoys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            List<T> decoyNucleicAcids = new List<T>();
            Parallel.ForEach(nucleicAcids, new ParallelOptions() { MaxDegreeOfParallelism = maxThreads }, nucleicAcid =>
            {
                // reverse sequence
                var reverseSequence =
                    new string(nucleicAcid.BaseSequence.Reverse().ToArray());

                // create a mapping of original to reversed indices
                var indexMapping = new Dictionary<int, int>();
                for (int i = 0; i < nucleicAcid.BaseSequence.Length; i++)
                {
                    indexMapping[i + 1] = nucleicAcid.BaseSequence.Length - i;
                }

                // reverse modifications; terminus-restricted mods stay on their terminus, otherwise a 5' mod
                // lands on the decoy's 3' residue (decoy digestion products skip location checks)
                var reverseModifications = new Dictionary<int, List<Modification>>();
                foreach (var kvp in nucleicAcid.OneBasedPossibleLocalizedModifications)
                {
                    foreach (var mod in kvp.Value)
                    {
                        var key = DecoyModificationKey(mod, kvp.Key, indexMapping);
                        if (!reverseModifications.TryGetValue(key, out var modsAtKey))
                            reverseModifications.Add(key, modsAtKey = new List<Modification>());
                        modsAtKey.Add(mod);
                    }
                }

                var reverseFixedModifications = new Dictionary<int, Modification>();
                if (nucleicAcid is IBioPolymer bioPol)
                {
                    foreach (var kvp in bioPol.OneBasedFixedModifications)
                    {
                        reverseFixedModifications[DecoyModificationKey(kvp.Value, kvp.Key, indexMapping)] = kvp.Value;
                    }
                }

                List<TruncationProduct> reverseTruncs = new List<TruncationProduct>();
                List<SequenceVariation> reverseVariations = new List<SequenceVariation>();
                List<SequenceVariation> reverseAppliedVariations = new List<SequenceVariation>();
                if (nucleicAcid is IHasSequenceVariants variantContaining)
                {
                    // Reverse Applied Variants
                    foreach (SequenceVariation variation in variantContaining.AppliedSequenceVariations)
                    {
                        var reverseBegin = indexMapping[variation.OneBasedBeginPosition];
                        var reverseEnd = indexMapping[variation.OneBasedEndPosition];
                        var reverseModificationsForVariation = new Dictionary<int, List<Modification>>();
                        foreach (var modKvp in variation.OneBasedModifications)
                        {
                            var reverseModKey = indexMapping[modKvp.Key];
                            reverseModificationsForVariation.Add(reverseModKey, modKvp.Value);
                        }
                        reverseAppliedVariations.Add(new SequenceVariation(reverseBegin, reverseEnd, variation.OriginalSequence, variation.VariantSequence, variation.VariantCallFormatDataString.Description, reverseModificationsForVariation));
                    }

                    // Reverse Applied Variants
                    foreach (SequenceVariation variation in variantContaining.SequenceVariations)
                    {
                        var reverseBegin = indexMapping[variation.OneBasedBeginPosition];
                        var reverseEnd = indexMapping[variation.OneBasedEndPosition];
                        var reverseModificationsForVariation = new Dictionary<int, List<Modification>>();
                        foreach (var modKvp in variation.OneBasedModifications)
                        {
                            var reverseModKey = indexMapping[modKvp.Key];
                            reverseModificationsForVariation.Add(reverseModKey, modKvp.Value);
                        }
                        reverseVariations.Add(new SequenceVariation(reverseBegin, reverseEnd, variation.OriginalSequence, variation.VariantSequence, variation.VariantCallFormatDataString.Description, reverseModificationsForVariation));
                    }

                    // Reverse Truncations
                    foreach (TruncationProduct truncation in variantContaining.TruncationProducts)
                    {
                        var reverseBegin = indexMapping[truncation.OneBasedEndPosition!.Value];
                        var reverseEnd = indexMapping[truncation.OneBasedBeginPosition!.Value];

                        reverseTruncs.Add(new(reverseBegin, reverseEnd, $"{decoyIdentifier} {truncation.Type}"));
                    }
                }

                T newNucleicAcid = nucleicAcid.CreateNew(reverseSequence, reverseModifications, true, reverseTruncs, reverseVariations, reverseAppliedVariations, decoyIdentifier, reverseFixedModifications);
                lock (decoyNucleicAcids)
                {
                    decoyNucleicAcids.Add(newNucleicAcid);
                }
            });
            return decoyNucleicAcids;
        }

        /// <summary>
        /// Decoy position of a modification: 5'- and 3'-terminal mods keep their own terminus, all others are reversed.
        /// </summary>
        private static int DecoyModificationKey(Modification mod, int oneBasedPosition, Dictionary<int, int> indexMapping)
        {
            return mod.LocationRestriction switch
            {
                "5'-terminal." or "Oligo 5'-terminal." or "N-terminal." or "Peptide N-terminal." when oneBasedPosition == 1 => 1,
                "3'-terminal." or "Oligo 3'-terminal." or "C-terminal." or "Peptide C-terminal." when oneBasedPosition == indexMapping.Count => indexMapping.Count,
                _ => indexMapping[oneBasedPosition]
            };
        }

        private static List<T> GenerateSlidedDecoys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }

        private static List<T> GenerateShuffledDeocys<T>(List<T> nucleicAcids, int maxThreads, string decoyIdentifier) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }

    }
}
