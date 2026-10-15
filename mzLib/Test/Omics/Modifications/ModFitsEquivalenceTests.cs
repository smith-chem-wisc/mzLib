using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// ModFits compares a single-residue motif with the residue at the index directly. These tests hold it to the
    /// previous implementation, kept verbatim below, over random motifs (one to four letters, with ambiguity codes and
    /// lower-case context), every location restriction, and indexes inside, at the edges of and outside the sequence.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ModFitsEquivalenceTests
    {
        private const string Letters = "ACDEFGHIKLMNPQRSTVWYXBJZU";
        private static readonly string[] Restrictions =
        {
            "Anywhere.", "Unassigned.", "N-terminal.", "5'-terminal.", "C-terminal.", "3'-terminal.",
            "Peptide N-terminal.", "Oligo 5'-terminal.", "Peptide C-terminal.", "Oligo 3'-terminal.",
        };

        [Test]
        public static void ModFitsAnswersAsBefore()
        {
            var random = new Random(15);
            var modifications = new List<Modification>();
            for (int m = 0; m < 400; m++)
            {
                modifications.Add(RandomModification(random, m));
            }
            Assert.That(modifications.Count(m => m.Target.ToString().Length == 1), Is.GreaterThan(80));

            long fits = 0, calls = 0;
            for (int trial = 0; trial < 3000; trial++)
            {
                int length = 1 + random.Next(60);
                char[] residues = new char[length];
                for (int r = 0; r < length; r++)
                {
                    residues[r] = Letters[random.Next(Letters.Length)];
                }
                string sequence = new(residues);

                foreach (Modification mod in modifications.OrderBy(_ => random.Next()).Take(40))
                {
                    for (int call = 0; call < 20; call++)
                    {
                        int bioPolymerIndex = random.Next(-2, length + 3);
                        int productLength = 1 + random.Next(length);
                        int productIndex = random.Next(-1, productLength + 2);

                        bool expected = ModFitsBefore(mod, sequence, productIndex, productLength, bioPolymerIndex);
                        bool actual = ModificationLocalization.ModFits(mod, sequence, productIndex, productLength, bioPolymerIndex);
                        Assert.That(actual, Is.EqualTo(expected),
                            $"motif {mod.Target} {mod.LocationRestriction} in {sequence} at {bioPolymerIndex} ({productIndex} of {productLength})");
                        fits += expected ? 1 : 0;
                        calls++;
                    }
                }
            }
            Assert.That(fits, Is.GreaterThan(calls / 100), "the random calls must hit both answers");
        }

        private static Modification RandomModification(Random random, int index)
        {
            int length = random.Next(3) == 0 ? 2 + random.Next(3) : 1;
            int capital = random.Next(length);
            char[] motif = new char[length];
            for (int i = 0; i < length; i++)
            {
                char letter = Letters[random.Next(Letters.Length)];
                motif[i] = i == capital ? letter : char.ToLowerInvariant(letter);
            }
            Assert.That(ModificationMotif.TryGetMotif(new string(motif), out ModificationMotif target));
            return new Modification(_originalId: "mod" + index, _modificationType: "Test", _target: target,
                _locationRestriction: Restrictions[random.Next(Restrictions.Length)], _monoisotopicMass: 1 + index);
        }

        #region the previous implementation, verbatim

        private static bool ModFitsBefore(Modification attemptToLocalize, string sequence, int digestionProductOneBasedIndex, int digestionProductLength, int bioPolymerOneBasedIndex)
        {
            // First find the capital letter...
            var motif = attemptToLocalize.Target.ToString();
            var motifStartLocation = -1;
            for (int i = 0; i < motif.Length; i++)
            {
                if (!char.IsUpper(motif[i]))
                    continue;

                motifStartLocation = i;
                break;
            }

            // Look up starting at and including the capital letter
            var proteinToMotifOffset = bioPolymerOneBasedIndex - motifStartLocation - 1;
            var indexUp = 0;
            while (indexUp < motif.Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= sequence.Length
                    || !MotifMatchesBefore(motif[indexUp], sequence[indexUp + proteinToMotifOffset]))
                {
                    return false;
                }
                indexUp++;
            }
            switch (attemptToLocalize.LocationRestriction)
            {
                    // Only the intact (undigested) terminus
                case "N-terminal." when bioPolymerOneBasedIndex > 2:
                case "5'-terminal." when bioPolymerOneBasedIndex > 2:
                case "C-terminal." when bioPolymerOneBasedIndex < sequence.Length:
                case "3'-terminal." when bioPolymerOneBasedIndex < sequence.Length:

                    // All Digested Termini AND original undigested termini
                case "Peptide N-terminal." when digestionProductOneBasedIndex > 1:
                case "Oligo 5'-terminal." when digestionProductOneBasedIndex > 1:
                case "Peptide C-terminal." when digestionProductOneBasedIndex < digestionProductLength:
                case "Oligo 3'-terminal." when digestionProductOneBasedIndex < digestionProductLength:
                    return false;

                default:
                    // I guess Anywhere. and Unassigned. are true since how do you localize anywhere or unassigned.

                    return true;
            }
        }

        private static bool MotifMatchesBefore(char motifChar, char sequenceChar)
        {
            char upperMotifChar = char.ToUpper(motifChar);
            return upperMotifChar switch
            {
                'X' => true,
                'B' => sequenceChar is 'D' or 'N',
                'J' => sequenceChar is 'I' or 'L',
                'Z' => sequenceChar is 'E' or 'Q',
                _ => upperMotifChar == sequenceChar
            };
        }

        #endregion
    }
}
