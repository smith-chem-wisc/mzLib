using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics.Modifications;
using Transcriptomics;

namespace UsefulProteomicsDatabases.Transcriptomics
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
        public static List<T> GenerateDecoys<T>(List<T> nucleicAcids, DecoyType decoyType, int maxThreads = -1) where T : INucleicAcid
        {
            switch (decoyType)
            {
                case DecoyType.None:
                    return new List<T>();
                case DecoyType.Reverse:
                    return GenerateReverseDecoys(nucleicAcids, maxThreads);
                case DecoyType.Slide:
                    return GenerateSlidedDecoys(nucleicAcids, maxThreads);
                case DecoyType.Shuffle:
                    return GenerateShuffledDeocys(nucleicAcids, maxThreads);
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
        private static List<T> GenerateReverseDecoys<T>(List<T> nucleicAcids, int maxThreads = -1) where T : INucleicAcid
        {
            List<T> decoyNucleicAcids = new List<T>();
            Parallel.ForEach(nucleicAcids, new ParallelOptions() { MaxDegreeOfParallelism = maxThreads }, nucleicAcid =>
            {
                // reverse sequence
                var reverseSequence =
                    new string(nucleicAcid.BaseSequence[..^1].Reverse().Append(nucleicAcid.BaseSequence.Last()).ToArray());

                // reverse modifications
                var reverseModifications = new Dictionary<int, List<Modification>>();
                foreach (var kvp in nucleicAcid.OneBasedPossibleLocalizedModifications)
                {
                    var reverseKey = kvp.Key == reverseSequence.Length ? kvp.Key : reverseSequence.Length - kvp.Key;
                    reverseModifications.Add(reverseKey, kvp.Value);
                }

                T newNucleicAcid = nucleicAcid.CreateNew(reverseSequence, reverseModifications, true);
                lock (decoyNucleicAcids)
                {
                    decoyNucleicAcids.Add(newNucleicAcid);
                }
            });
            return decoyNucleicAcids;
        }

        private static List<T> GenerateSlidedDecoys<T>(List<T> nucleicAcids, int maxThreads = -1) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }

        private static List<T> GenerateShuffledDeocys<T>(List<T> nucleicAcids, int maxThreads = -1) where T : INucleicAcid
        {
            throw new NotImplementedException();
        }

    }
}
