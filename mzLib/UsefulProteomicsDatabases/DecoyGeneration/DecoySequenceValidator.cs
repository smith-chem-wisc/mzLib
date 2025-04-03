#nullable enable
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases.DecoyGeneration;

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
    public static T ScrambleDecoyBioPolymer<T>(T originalDecoy, IDigestionParams digestionParams,
        HashSet<string> forbiddenSequences, IEnumerable<string>? sequencesToScramble = null) where T : IBioPolymer
    {
        var random = new Random(42);

        // If no sequencesToScramble are passed in, we check to see if any 
        // peptides in the decoy are forbidden sequences
        sequencesToScramble ??= originalDecoy
            .Digest(digestionParams, new List<Modification>(), new List<Modification>())
            .Select(pep => pep.FullSequence)
            .Where(forbiddenSequences.Contains);

        if (!sequencesToScramble.Any())
            return originalDecoy;

        string scrambledSequence = originalDecoy.BaseSequence;

        // Clone the original protein's modifications
        var scrambledModificationDictionary = originalDecoy.OriginalNonVariantModifications.ToDictionary(kvp => kvp.Key, kvp => kvp.Value);


        // Start small and then go big. If we scramble a zero-missed cleavage peptide, but the missed cleavage peptide contains the previously scrambled peptide
        // Then we can avoid unnecessary operations as the scrambledSequence will no longer contain the longer sequence of the missed cleavage peptide
        foreach (string peptideSequence in sequencesToScramble.OrderBy(seq => seq.Length))
        {
            if (scrambledSequence.Contains(peptideSequence))
            {
                string scrambledPeptideSequence = ScrambleSequence(peptideSequence, digestionParams.DigestionAgent.DigestionMotifs, random,
                    out var swappedArray);
                int scrambleAttempts = 1;

                // Try five times to scramble the peptide sequence without creating a forbidden sequence
                while (forbiddenSequences.Contains(scrambledPeptideSequence) & scrambleAttempts <= 5)
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
        return (T)newDecoy;
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
}