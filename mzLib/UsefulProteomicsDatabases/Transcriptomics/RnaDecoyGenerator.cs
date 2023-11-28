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
    /// Generates Decoy Nucleic Acids from any implementor of INucleicAcid
    /// TODO: Implement Shuffle and Slide Decoys
    /// TODO: Consider passing digestion motif as optional parameter to leave digestion sites intact
    /// </summary>
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
        /// leaving modification on their nucleic acid of origin, and 3' termini intact
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
