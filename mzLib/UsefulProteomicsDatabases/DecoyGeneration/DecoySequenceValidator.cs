#nullable enable
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Transcriptomics;

namespace UsefulProteomicsDatabases;

public static class DecoySequenceValidator
{
    /// <summary>
    /// This function takes in a decoy protein and a list of forbidden sequences that the decoy
    /// protein should not contain. Optionally, a list of the peptides within the base sequence
    /// of the decoy protein that need to be scrambled can be passed as well. It will scramble the required sequences,
    /// leaving cleavage sites intact. 
    /// </summary>
    /// <param name="originalDecoy"> A Decoy protein to be cloned </param>
    /// <param name="digestionParams"> Digestion parameters </param>
    /// <param name="forbiddenSequences"> A HashSet of forbidden sequences that the decoy protein should not contain. Typically, a set of target base sequences </param>
    /// <param name="sequencesToScramble"> Optional IEnumberable of sequences within the decoy protein that need to be replaced.
    ///                                     If this is passed, only sequences within the IEnumerable will be replaced!!! </param>
    /// <returns> A cloned copy of the decoy protein with a scrambled sequence </returns>
    public static TBioPolymer ScrambleDecoyBioPolymer<TBioPolymer>(TBioPolymer originalDecoy, IDigestionParams digestionParams,
        HashSet<string> forbiddenSequences, IEnumerable<string>? sequencesToScramble = null) where TBioPolymer : IBioPolymer
    {
        var random = new Random(42);

        // If no sequencesToScramble are passed in, we check to see if any 
        // peptides in the decoy are forbidden sequences
        var digested = originalDecoy
            .Digest(digestionParams, new List<Modification>(), new List<Modification>())
            .Select(pep => pep.FullSequence)
            .ToHashSet();

        HashSet<string> toScramble = sequencesToScramble?.ToHashSet() ?? digested.Where(forbiddenSequences.Contains).ToHashSet();

        // if it's a nucleic acid, we need to check for palindromic sequences and scramble those
        if (originalDecoy is NucleicAcid)
        {
            // If the sequence is palindromic, and not in the sequencesToScramble, we need to scramble it
            foreach (var sequenceToInvestigate in digested.Except(toScramble))
            {
                // Check if the sequence is palindromic
                if (IsPalindromic(sequenceToInvestigate, out int degreeOfPalindromicity, 3)) // 3 chosen as a magic number to ensure the first 3 terminal ions have unique mass, this should be refined with experimental evidence
                {
                    // If it is, we need to scramble it
                    toScramble.Add(sequenceToInvestigate);
                }
            }
        }

        if (!toScramble.Any())
            return originalDecoy;

        string scrambledSequence = originalDecoy.BaseSequence;

        // Clone the original protein's modifications
        var scrambledModificationDictionary = originalDecoy.OriginalNonVariantModifications.ToDictionary(kvp => kvp.Key, kvp => kvp.Value);

        // Start small and then go big. If we scramble a zero-missed cleavage peptide, but the missed cleavage peptide contains the previously scrambled peptide
        // Then we can avoid unnecessary operations as the scrambledSequence will no longer contain the longer sequence of the missed cleavage peptide
        foreach (string peptideSequence in toScramble.OrderBy(seq => seq.Length))
        {
            if (scrambledSequence.Contains(peptideSequence))
            {
                string scrambledPeptideSequence = ScrambleSequence(peptideSequence, digestionParams.DigestionAgent.DigestionMotifs, random,
                    out var swappedArray);
                int scrambleAttempts = 1;

                // Try five times to scramble the peptide sequence without creating a forbidden sequence
                while ((forbiddenSequences.Contains(scrambledPeptideSequence) || (originalDecoy is NucleicAcid && IsPalindromic(scrambledPeptideSequence, out _, 3))) 
                    & scrambleAttempts <= 5)
                {
                    scrambledPeptideSequence = ScrambleSequence(peptideSequence, digestionParams.DigestionAgent.DigestionMotifs, random,
                        out swappedArray);
                    scrambleAttempts++;
                }

                scrambledSequence = scrambledSequence.Replace(peptideSequence, scrambledPeptideSequence);

                if (!scrambledModificationDictionary.Any()) continue;

                // rearrange the modifications 
                foreach (int index in scrambledSequence.IndexOfAll(scrambledPeptideSequence))
                {
                    // Get mods that were affected by the scramble
                    var relevantMods = scrambledModificationDictionary.Where(kvp =>
                        kvp.Key >= index + 1 && kvp.Key < index + peptideSequence.Length + 1).ToList();

                    // Modify the dictionary to reflect the new positions of the modifications
                    foreach (var kvp in relevantMods)
                    {
                        int newKey = swappedArray[kvp.Key - 1 - index] + 1 + index;
                        // To prevent collisions, we have to check if mods already exist at the new idx.
                        if (scrambledModificationDictionary.TryGetValue(newKey, out var modsToSwap))
                        {
                            // If there are mods at the new idx, we swap the mods
                            scrambledModificationDictionary[newKey] = kvp.Value;
                            scrambledModificationDictionary[kvp.Key] = modsToSwap;
                        }
                        else
                        {
                            scrambledModificationDictionary.Add(newKey, kvp.Value);
                            scrambledModificationDictionary.Remove(kvp.Key);
                        }
                    }
                }
            }
        }

        IBioPolymer newDecoy = originalDecoy.CloneWithNewSequenceAndMods(scrambledSequence, scrambledModificationDictionary);
        return (TBioPolymer)newDecoy;
    }

    /// <summary>
    /// Scrambles a peptide sequence, preserving the position of any cleavage sites.
    /// </summary>
    /// <param name="swappedPositionArray">An array that maps the previous position (index) to the new position (value)</param>
    public static string ScrambleSequence(string sequence, List<DigestionMotif> motifs, Random rng, out int[] swappedPositionArray)
    {
        // First, find the location of every cleavage motif. These sites shouldn't be scrambled.
        HashSet<int> zeroBasedCleavageSitesLocations = new();
        foreach (var motif in motifs)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                (bool fits, bool prevents) = motif.Fits(sequence, i);
                if (fits && !prevents)
                {
                    zeroBasedCleavageSitesLocations.Add(i);
                }
            }
        }

        // Next, scramble the sequence using the Fisher-Yates shuffle algorithm.
        char[] sequenceArray = sequence.ToCharArray();
        // We're going to keep track of the positions of the characters in the original sequence,
        // This will enable us to adjust the location of modifications that are present in the original sequence
        // to the new scrambled sequence.
        int[] tempPositionArray = Enumerable.Range(0, sequenceArray.Length).ToArray();
        int n = sequenceArray.Length;
        while (n > 1)
        {
            n--;
            if (zeroBasedCleavageSitesLocations.Contains(n))
            {
                // Leave the cleavage site in place
                continue;
            }
            int k = rng.Next(n + 1);
            // don't swap the position of a cleavage site
            while (zeroBasedCleavageSitesLocations.Contains(k))
            {
                k = rng.Next(n + 1);
            }

            // rearrange the sequence array
            char tempResidue = sequenceArray[k];
            sequenceArray[k] = sequenceArray[n];
            sequenceArray[n] = tempResidue;

            // update the position array to represent the swaps
            int tempPosition = tempPositionArray[k];
            tempPositionArray[k] = tempPositionArray[n];
            tempPositionArray[n] = tempPosition;
        }

        // This maps the previous position (index) to the new position (value)
        swappedPositionArray = new int[tempPositionArray.Length];
        for (int i = 0; i < tempPositionArray.Length; i++)
        {
            swappedPositionArray[tempPositionArray[i]] = i;
        }

        return new string(sequenceArray);
    }

    /// <summary>
    /// Every zero-based position that participates in a cleavage-motif match: the FULL span of each
    /// match, not only the location at which it starts.
    /// </summary>
    /// <remarks>
    /// <see cref="ScrambleSequence"/> records only the match location. That is complete for a
    /// single-residue motif ("K", "R", "E"), but not for a multi-residue one: StcE's "TX|T" matches
    /// three residues, and leaving the trailing two free lets a rearrangement destroy that cleavage
    /// site and invent others elsewhere. Rearrangements that must round-trip through digestion --
    /// entrapment generation -- hold the whole span still, so they use this.
    /// </remarks>
    public static HashSet<int> CleavageSitePositions(string sequence, List<DigestionMotif> motifs)
    {
        HashSet<int> positions = new();
        if (string.IsNullOrEmpty(sequence) || motifs is null)
        {
            return positions;
        }

        foreach (DigestionMotif motif in motifs)
        {
            int span = motif.InducingCleavage?.Length ?? 0;
            if (span == 0)
            {
                continue;
            }

            for (int i = 0; i < sequence.Length; i++)
            {
                (bool fits, bool prevents) = motif.Fits(sequence, i);
                if (!fits || prevents)
                {
                    continue;
                }

                for (int offset = 0; offset < span && i + offset < sequence.Length; offset++)
                {
                    positions.Add(i + offset);
                }
            }
        }

        return positions;
    }

    /// <summary>
    /// The number of DISTINCT rearrangements of <paramref name="sequence"/> that leave every
    /// cleavage-site residue in place: the multinomial over the free positions, in closed form.
    /// Nothing is enumerated, so this is cheap even when the answer is astronomically large.
    /// </summary>
    /// <remarks>
    /// Returns <see cref="BigInteger.One"/> when no rearrangement exists -- a homopolymeric tract
    /// such as "SSSSSSR" has exactly one arrangement once its cleavage sites are pinned. That is
    /// arithmetic rather than a collision, so no amount of searching or reseeding can improve it.
    /// Intended for peptide-length sequences; the free-position count drives the cost.
    /// </remarks>
    public static BigInteger PermutationSpaceSize(string sequence, List<DigestionMotif> motifs)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return BigInteger.One;
        }

        SortedDictionary<char, int> counts = FreeResidueCounts(sequence, motifs, out int freeCount);
        return Multinomial(counts, freeCount);
    }

    /// <summary>
    /// The rearrangement of <paramref name="sequence"/> at <paramref name="index"/> in lexicographic
    /// order over the distinct rearrangements of its free (non-cleavage-site) positions.
    /// </summary>
    /// <param name="swappedPositionArray">Maps the previous position (index) to the new position
    /// (value), the same contract as <see cref="ScrambleSequence"/>, so modification remapping is
    /// shared between the two.</param>
    /// <remarks>
    /// The result is composition-preserving, hence isomeric with the input: same residues, same
    /// mass, different order. Unlike <see cref="ScrambleSequence"/> this consumes no random state,
    /// so a given (sequence, index) always yields the same string -- across runs, threads, machines
    /// and framework versions. Callers that need an unpredictable-but-reproducible choice derive
    /// the index from a seed and the sequence rather than from a random number generator.
    /// </remarks>
    /// <exception cref="MzLibException"><paramref name="index"/> lies outside the space.</exception>
    public static string UnrankPermutation(string sequence, List<DigestionMotif> motifs, BigInteger index,
        out int[] swappedPositionArray)
    {
        swappedPositionArray = Enumerable.Range(0, sequence?.Length ?? 0).ToArray();
        if (string.IsNullOrEmpty(sequence))
        {
            return sequence;
        }

        HashSet<int> pinned = CleavageSitePositions(sequence, motifs);
        List<int> freePositions = new();
        for (int i = 0; i < sequence.Length; i++)
        {
            if (!pinned.Contains(i))
            {
                freePositions.Add(i);
            }
        }

        SortedDictionary<char, int> counts = new();
        Dictionary<char, Queue<int>> origins = new();
        foreach (int position in freePositions)
        {
            char residue = sequence[position];
            counts.TryGetValue(residue, out int seen);
            counts[residue] = seen + 1;

            if (!origins.TryGetValue(residue, out Queue<int>? queue))
            {
                queue = new Queue<int>();
                origins[residue] = queue;
            }
            queue.Enqueue(position);
        }

        BigInteger size = Multinomial(counts, freePositions.Count);
        if (index < BigInteger.Zero || index >= size)
        {
            throw new MzLibException(
                $"Permutation index {index} is outside the {size} distinct permutations of '{sequence}'.");
        }

        char[] rearranged = sequence.ToCharArray();
        List<char> residues = counts.Keys.ToList();
        BigInteger remainingPermutations = size;
        int remainingPositions = freePositions.Count;
        BigInteger rank = index;

        foreach (int slot in freePositions)
        {
            foreach (char residue in residues)
            {
                if (counts[residue] == 0)
                {
                    continue;
                }

                // Permutations whose next residue is this one. Exact: the division cancels.
                BigInteger branch = remainingPermutations * counts[residue] / remainingPositions;
                if (rank < branch)
                {
                    rearranged[slot] = residue;
                    swappedPositionArray[origins[residue].Dequeue()] = slot;
                    counts[residue]--;
                    remainingPermutations = branch;
                    remainingPositions--;
                    break;
                }

                rank -= branch;
            }
        }

        return new string(rearranged);
    }

    /// <summary>Residue counts over the positions no cleavage motif holds in place.</summary>
    private static SortedDictionary<char, int> FreeResidueCounts(string sequence, List<DigestionMotif> motifs,
        out int freeCount)
    {
        HashSet<int> pinned = CleavageSitePositions(sequence, motifs);
        SortedDictionary<char, int> counts = new();
        freeCount = 0;

        for (int i = 0; i < sequence.Length; i++)
        {
            if (pinned.Contains(i))
            {
                continue;
            }

            counts.TryGetValue(sequence[i], out int seen);
            counts[sequence[i]] = seen + 1;
            freeCount++;
        }

        return counts;
    }

    /// <summary>
    /// n! / (m1! * m2! * ...), built from successive binomials so no factorial is ever materialised.
    /// </summary>
    private static BigInteger Multinomial(SortedDictionary<char, int> counts, int total)
    {
        BigInteger permutations = BigInteger.One;
        int placed = 0;

        foreach (int count in counts.Values)
        {
            placed += count;
            permutations *= Binomial(placed, count);
        }

        return total == 0 ? BigInteger.One : permutations;
    }

    /// <summary>Binomial coefficient, computed incrementally so every intermediate stays exact.</summary>
    private static BigInteger Binomial(int n, int k)
    {
        if (k < 0 || k > n)
        {
            return BigInteger.Zero;
        }

        k = Math.Min(k, n - k);
        BigInteger result = BigInteger.One;
        for (int i = 1; i <= k; i++)
        {
            result = result * (n - k + i) / i;
        }

        return result;
    }

    /// <summary>
    /// Determines if a string is palindromic and calculates the degree of palindromicity.
    /// </summary>
    /// <param name="input">The input string to check.</param>
    /// <param name="palindromicCharacters">The number of palindromic characters.</param>
    /// <param name="characterCountCutoff">Number of palindromic residues required to return true, ABCDA counts as 1, ABCDBA counts as 2, and ABCBA counts as 3 </param>
    /// <returns>True if the string is palindromic, otherwise false.</returns>
    public static bool IsPalindromic(string input, out int palindromicCharacters, int? characterCountCutoff = null )
    {
        palindromicCharacters = 0;
        if (string.IsNullOrEmpty(input))
        {
            return false;
        }

        // if null, cutoff will ensure only fully palindromic sequences are returned as true
        characterCountCutoff ??= (input.Length + 1) / 2;

        int left = 0;
        int right = input.Length - 1;

        while (left <= right)
        {
            char leftChar = input[left];
            char rightChar = input[right];
            if (leftChar == rightChar)
            {
                palindromicCharacters++;
            }
            else
            {
                break;
            }

            left++;
            right--;
        }

        return palindromicCharacters >= characterCountCutoff;
    }
}