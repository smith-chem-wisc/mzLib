using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using UsefulProteomicsDatabases;
using Transcriptomics.Digestion;
using Omics.Digestion;
using Transcriptomics;
using System.Diagnostics.CodeAnalysis;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class DecoySequenceValidatorTests
{
    [TestCase("AABBAA", 3, 3, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Equal")]
    [TestCase("AABBAA", 3, null, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Null")]
    [TestCase("AABBAA", 3, 4, false, TestName = "Full Palindrome-Even Sequence Length-Cutoff Greater")]
    [TestCase("AABBAA", 3, 2, true, TestName = "Full Palindrome-Even Sequence Length-Cutoff Lesser")]
    [TestCase("ABCBA", 3, 4, false, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Greater")]
    [TestCase("ABCBA", 3, null, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Null")]
    [TestCase("ABCBA", 3, 3, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Equal")]
    [TestCase("ABCBA", 3, 2, true, TestName = "Full Palindrome-Odd Sequence Length-Cutoff Lesser")]
    [TestCase("ABCDEFCBA", 3, 2, true, TestName = "Partial Palindrome-Odd Sequence-Cutoff Lesser")]
    [TestCase("ABCDEFCBA", 3, null, false, TestName = "Partial Palindrome-Odd Sequence-Cutoff Null")]
    [TestCase("ABCDEFCBA", 3, 3, true, TestName = "Partial Palindrome-Odd Sequence-Cutoff Equal")]
    [TestCase("ABCDEFCBA", 3, 4, false, TestName = "Partial Palindrome-Odd Sequence-Cutoff Greater")]
    [TestCase("ABCDEGFCBA", 3, 2, true, TestName = "Partial Palindrome-Even Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGFCBA", 3, null, false, TestName = "Partial Palindrome-Even Sequence-Cutoff Null")]
    [TestCase("ABCDEGFCBA", 3, 3, true, TestName = "Partial Palindrome-Even Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFCBA", 3, 4, false, TestName = "Partial Palindrome-Even Sequence-Cutoff Greater")]
    [TestCase("ABCDEGF", 0, 2, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGF", 0, null, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Null")]
    [TestCase("ABCDEGF", 0, 3, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFD", 0, 4, false, TestName = "Non-Palindrome-Odd Sequence-Cutoff Greater")]
    [TestCase("ABCDEGFD", 0, 2, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Lesser")]
    [TestCase("ABCDEGFD", 0, null, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Null")]
    [TestCase("ABCDEGFD", 0, 3, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Equal")]
    [TestCase("ABCDEGFD", 0, 4, false, TestName = "Non-Palindrome-Even Sequence-Cutoff Greater")]
    [TestCase("", 0, 1, false, TestName = "Empty Sequence")]
    [TestCase(null, 0, 1, false, TestName = "Null Sequence")]
    public void IsPalindromicTests(string input, int expectedDegree, int? degreeCutoff, bool expectedResult)
    {
        bool result = DecoySequenceValidator.IsPalindromic(input, out int degreeOfPalindromicity, degreeCutoff);
        Assert.That(result, Is.EqualTo(expectedResult));
        Assert.That(degreeOfPalindromicity, Is.EqualTo(expectedDegree));
    }


    [Test]
    public static void TestDecoyScramblerModificationHandling_Rna()
    {
        IDigestionParams d = new RnaDigestionParams("RNase T1",
                    maxMissedCleavages: 1,
                    minLength: 5);

        ModificationMotif.TryGetMotif("G", out ModificationMotif motifG);
        ModificationMotif.TryGetMotif("A", out ModificationMotif motifA);
        Modification modG = new Modification("myMod", null, "myModType", null, motifG, "Anywhere.", null, 10, null, null, null, null, null, null);
        Modification modA = new Modification("myMod", null, "myModType", null, motifA, "Anywhere.", null, 10, null, null, null, null, null, null);

        IDictionary<int, List<Modification>> modDictDecoy = new Dictionary<int, List<Modification>>
            {
                {9, new List<Modification> { modG } },
                {7, new List<Modification> { modA } }
            };

        RNA target = new RNA("GUUUAUUUGUAUUUUUU", "target");
        RNA decoy = new RNA("UUUUUUAUGUUUAUUUG", "decoy", modDictDecoy);

        var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
        var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

        HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();
        var offendingDecoys = decoyPep.Where(p => targetPepSeqs.Contains(p.FullSequence)).Select(d => d.FullSequence).ToList();
        Assert.That(offendingDecoys.Count, Is.EqualTo(1));

        RNA scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs, offendingDecoys);
        offendingDecoys = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>())
            .Where(p => targetPepSeqs.Contains(p.FullSequence))
            .Select(d => d.FullSequence)
            .ToList();
        Assert.That(offendingDecoys.Count, Is.EqualTo(0));


        var aIndex = scrambledDecoy.BaseSequence.IndexOf("A");
        var gIndex = scrambledDecoy.BaseSequence.IndexOf("G"); // We modified the first residue, so we don't need all locations, just the first
        var aIndices = scrambledDecoy.BaseSequence.IndexOfAll("A");
        var gIndices = scrambledDecoy.BaseSequence.IndexOfAll("G");

        Assert.That(gIndices.Count(), Is.EqualTo(2));
        Assert.That(aIndices.Count(), Is.EqualTo(2));
        Assert.That(aIndices.First(), Is.EqualTo(aIndex));

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.ContainsKey(aIndex + 1), Is.True);
        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications[aIndex + 1].Contains(modA), Is.True);

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.ContainsKey(gIndex + 1), Is.True);
        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications[gIndex + 1].Contains(modG), Is.True);

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
    }

    [Test]
    [TestCase("GAACCAAGACGUACGUACACACG")]
    [TestCase("GAACCAAGAACCAAGGAACCAAGGAACCAAG")]
    [TestCase("GAAUAAGAAUAAGAAUAAGAAUAAG")]
    public static void ScramblesNucleicAcidPalindrome(string sequence)
    { 
        var rna = new RNA(sequence, "test");
        var digestionParams = new RnaDigestionParams("RNase T1", maxMissedCleavages: 1, minLength: 5);
        var oligos = rna.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

        bool oneIsPalindromic = false;
        int palindromeCount = oligos.Count(p => DecoySequenceValidator.IsPalindromic(p.BaseSequence, out _));  
        Assert.That(palindromeCount, Is.GreaterThan(0));


        // scramble with no forbidden sequences, so only the palindromic oligo component should be scrambled
        var decoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(rna, digestionParams, [], new List<string>());
        var decoyOligos = decoy.Digest(digestionParams, new List<Modification>(), new List<Modification>());
        palindromeCount = decoyOligos.Count(p => DecoySequenceValidator.IsPalindromic(p.BaseSequence, out _));
        Assert.That(palindromeCount, Is.EqualTo(0));
    }

    [Test]
    public static void TestDecoyScramblerModificationHandling_Protein()
    {
        DigestionParams d = new DigestionParams(
                    maxMissedCleavages: 1,
                    minPeptideLength: 5,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

        ModificationMotif.TryGetMotif("G", out ModificationMotif motifG);
        ModificationMotif.TryGetMotif("F", out ModificationMotif motifF);
        Modification modG = new Modification("myMod", null, "myModType", null, motifG, "Anywhere.", null, 10, null, null, null, null, null, null);
        Modification modF = new Modification("myMod", null, "myModType", null, motifF, "Anywhere.", null, 10, null, null, null, null, null, null);

        IDictionary<int, List<Modification>> modDictDecoy = new Dictionary<int, List<Modification>>
            {
                {8, new List<Modification> { modG } },
                {10, new List<Modification> { modF } }
            };

        Protein target = new Protein("MEDEEKFVGYKYGVFK", "target"); //, oneBasedModifications: modDictTarget);
        Protein decoy = new Protein("EEDEMKYGVFKFVGYK", "decoy", oneBasedModifications: modDictDecoy);

        var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
        var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

        HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();
        var offendingDecoys = decoyPep.Where(p => targetPepSeqs.Contains(p.FullSequence)).Select(d => d.FullSequence).ToList();
        Protein scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs, offendingDecoys);

        var fIndex = scrambledDecoy.BaseSequence.IndexOf("F");
        var gIndex = scrambledDecoy.BaseSequence.IndexOf("G"); // We modified the first residue, so we don't need all locations, just the first
        var fIndices = scrambledDecoy.BaseSequence.IndexOfAll("F");
        var gIndices = scrambledDecoy.BaseSequence.IndexOfAll("G");

        Assert.That(gIndices.Count(), Is.EqualTo(2));
        Assert.That(fIndices.Count(), Is.EqualTo(2));
        Assert.That(fIndices.First(), Is.EqualTo(fIndex));

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.ContainsKey(fIndex + 1), Is.True);
        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications[fIndex + 1].Contains(modF), Is.True);

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.ContainsKey(gIndex + 1), Is.True);
        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications[gIndex + 1].Contains(modG), Is.True);

        Assert.That(scrambledDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
    }

    [Test, Timeout(5000)]
    public static void TestDecoyScramblerNoInfiniteLoops()
    {
        DigestionParams d = new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 3,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

        Protein target = new Protein("MEK", "target");
        Protein decoy = new Protein("EMK", "decoy");

        var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
        var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

        HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();

        // We'll pretend that this is also a target sequence and can't be used as a decoy
        HashSet<string> offendingDecoys = new HashSet<string> { "EMK" };

        // You can't win in this scenario, there's no way to scramble that results in a different decoy
        Protein scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs.Union(offendingDecoys).ToHashSet(), offendingDecoys);
        var scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

        Assert.That(decoyPep.Count(), Is.EqualTo(scrambledPep.Count()));

        d = new DigestionParams(
                    maxMissedCleavages: 1,
                    minPeptideLength: 3,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

        offendingDecoys = new HashSet<string> { "KEK" };

        var impossibleDecoy = new Protein("KEK", "target"); // This guy could crash the shuffling algorithm
        scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(impossibleDecoy, d, offendingDecoys, offendingDecoys);

        Assert.That("KEK", Is.EqualTo(scrambledDecoy.BaseSequence));
    }

}