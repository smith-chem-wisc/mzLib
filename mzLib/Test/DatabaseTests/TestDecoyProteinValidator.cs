using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
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

    // [Timeout] over [CancelAfter] deliberately: CancelAfter cannot enforce a budget on a
    // non-cooperative synchronous test, so it would not catch the infinite loop this guards against.
#pragma warning disable CS0618
    [Test, Timeout(5000)]
#pragma warning restore CS0618
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

    #region Deterministic permutation (entrapment generation)

    private static List<DigestionMotif> Trypsin => DigestionMotif.ParseDigestionMotifsFromString("K|,R|");

    [Test]
    public void PermutationSpaceSize_IsTheMultinomialOverFreePositions()
    {
        // No K or R, so all nine residues are free: T x6, P x2, A x1 -> 9!/(6!*2!*1!) = 252.
        Assert.That(DecoySequenceValidator.PermutationSpaceSize("TTTPAPTTT", Trypsin), Is.EqualTo(new BigInteger(252)));

        // R is pinned and the six free residues are identical, so exactly one arrangement exists.
        // This is not a collision -- it is arithmetic, and no amount of searching can help.
        Assert.That(DecoySequenceValidator.PermutationSpaceSize("SSSSSSR", Trypsin), Is.EqualTo(BigInteger.One));

        // Seven distinct free residues with K pinned -> 7! = 5040.
        Assert.That(DecoySequenceValidator.PermutationSpaceSize("ACDEFGHK", Trypsin), Is.EqualTo(new BigInteger(5040)));
    }

    [Test]
    public void UnrankPermutation_IsDeterministicAndEnumeratesTheSpaceExactly()
    {
        const string sequence = "TTTPAPTTT";
        BigInteger size = DecoySequenceValidator.PermutationSpaceSize(sequence, Trypsin);
        var distinct = new HashSet<string>();

        for (BigInteger i = BigInteger.Zero; i < size; i++)
        {
            string once = DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, i, out _);
            string twice = DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, i, out _);

            Assert.That(twice, Is.EqualTo(once), "the same index must always give the same permutation");
            Assert.That(string.Concat(once.OrderBy(c => c)), Is.EqualTo(string.Concat(sequence.OrderBy(c => c))),
                "the permutation must be composition-preserving, hence isomeric with its target");
            distinct.Add(once);
        }

        Assert.That(distinct.Count, Is.EqualTo(252), "every index must give a different permutation");
    }

    [Test]
    public void UnrankPermutation_LeavesEveryCleavageSiteInPlace()
    {
        const string sequence = "ACDKEFGRHK";
        BigInteger size = DecoySequenceValidator.PermutationSpaceSize(sequence, Trypsin);

        for (BigInteger i = BigInteger.Zero; i < BigInteger.Min(size, 200); i++)
        {
            string permuted = DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, i, out _);
            Assert.That(permuted[3], Is.EqualTo('K'));
            Assert.That(permuted[7], Is.EqualTo('R'));
            Assert.That(permuted[9], Is.EqualTo('K'));
        }
    }

    [Test]
    public void UnrankPermutation_SwappedPositionArrayMapsOldIndexToNew()
    {
        // Same contract as ScrambleSequence's out parameter: previous position (index) -> new position (value).
        const string sequence = "ACDEFGHK";
        string permuted = DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, new BigInteger(1234), out int[] swapped);

        Assert.That(swapped.Length, Is.EqualTo(sequence.Length));
        for (int i = 0; i < sequence.Length; i++)
        {
            Assert.That(permuted[swapped[i]], Is.EqualTo(sequence[i]),
                "the residue at old index " + i + " must be found at its mapped new index");
        }
    }

    [Test]
    public void UnrankPermutation_RejectsAnIndexOutsideTheSpace()
    {
        BigInteger size = DecoySequenceValidator.PermutationSpaceSize("TTTPAPTTT", Trypsin);

        // The message names the offending index and the size of the space; a caller probing for a
        // free permutation reads it to tell "exhausted" from "asked for the wrong thing".
        var toolarge = Assert.Throws<MzLibException>(() =>
            DecoySequenceValidator.UnrankPermutation("TTTPAPTTT", Trypsin, size, out _));
        Assert.That(toolarge!.Message, Does.Contain("252").And.Contain("TTTPAPTTT"));

        var negative = Assert.Throws<MzLibException>(() =>
            DecoySequenceValidator.UnrankPermutation("TTTPAPTTT", Trypsin, BigInteger.MinusOne, out _));
        Assert.That(negative!.Message, Does.Contain("-1"));
    }

    [Test]
    public void PermutationOfAnEmptyOrUnmotifedSequenceIsWellDefined()
    {
        // An empty sequence has exactly one arrangement, not zero, and must not throw.
        Assert.That(DecoySequenceValidator.PermutationSpaceSize("", Trypsin), Is.EqualTo(BigInteger.One));
        Assert.That(DecoySequenceValidator.UnrankPermutation("", Trypsin, BigInteger.Zero, out int[] swapped),
            Is.EqualTo(""));
        Assert.That(swapped, Is.Empty);

        // No motifs at all -- a top-down "protease" -- pins nothing, so every residue is free.
        var noMotifs = new List<DigestionMotif>();
        Assert.That(DecoySequenceValidator.CleavageSitePositions("ACDEFGHK", noMotifs), Is.Empty);
        Assert.That(DecoySequenceValidator.PermutationSpaceSize("ACDEFGHK", noMotifs),
            Is.EqualTo(new BigInteger(40320)));   // 8! -- the K is free now
    }

    [Test]
    public void CleavageSitePositions_AreEmptyForAnEmptySequenceOrNoMotifs()
    {
        Assert.That(DecoySequenceValidator.CleavageSitePositions("", Trypsin), Is.Empty);
        Assert.That(DecoySequenceValidator.CleavageSitePositions(null!, Trypsin), Is.Empty);
        Assert.That(DecoySequenceValidator.CleavageSitePositions("ACDEFGHK", null!), Is.Empty);
    }

    [Test]
    public void CleavageSitePositions_NeverReturnAnIndexOutsideTheSequence()
    {
        // A multi-residue motif matching at the very end must not run off the end of the sequence.
        var stcE = DigestionMotif.ParseDigestionMotifsFromString("TX|T");
        foreach (string sequence in new[] { "ATPT", "TPT", "TTTPAPTTT", "AAATPTGGGSQSCCC" })
        {
            foreach (int position in DecoySequenceValidator.CleavageSitePositions(sequence, stcE))
            {
                Assert.That(position, Is.InRange(0, sequence.Length - 1),
                    "position " + position + " is outside '" + sequence + "'");
            }
        }
    }

    // Golden vectors. The properties asserted above -- exhaustive enumeration, preserved
    // composition, fixed cleavage sites -- all hold for ANY consistent ordering, so none of them
    // would notice if the ordering itself changed. Entrapment pairing depends on a given index
    // always yielding the same string, so the ordering is part of the contract and is pinned here.
    // Expected values come from an independent implementation, not from a snapshot of this one.
    [TestCase("ACDEFGHK", 0, "ACDEFGHK", TestName = "Unrank golden - first index is lexicographically smallest")]
    [TestCase("ACDEFGHK", 1234, "CGDEHAFK", TestName = "Unrank golden - interior index")]
    [TestCase("ACDEFGHK", 5039, "HGFEDCAK", TestName = "Unrank golden - last index is lexicographically largest")]
    [TestCase("TTTPAPTTT", 0, "APPTTTTTT", TestName = "Unrank golden - repeated residues, first index")]
    [TestCase("TTTPAPTTT", 251, "TTTTTTPPA", TestName = "Unrank golden - repeated residues, last index")]
    [TestCase("ACDKEFGRHK", 0, "ACDKEFGRHK", TestName = "Unrank golden - interior cleavage sites stay put")]
    public void UnrankPermutation_OrderingIsPartOfTheContract(string sequence, int index, string expected)
    {
        Assert.That(DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, new BigInteger(index), out _),
            Is.EqualTo(expected));
    }

    [Test]
    public void CleavageSitePositions_PinTheFullMotifSpanNotJustItsFirstResidue()
    {
        // StcE's TX|T matches three residues. ScrambleSequence pins only the match location, so a
        // scramble there can move the X and the trailing T -- destroying a real cleavage site and
        // inventing others. The entrapment path must pin the whole span.
        var stcE = DigestionMotif.ParseDigestionMotifsFromString("TX|T");
        var pinned = DecoySequenceValidator.CleavageSitePositions("ATPTA", stcE);

        Assert.That(pinned.OrderBy(p => p), Is.EqualTo(new[] { 1, 2, 3 }));
    }

    [Test]
    public void CleavageSitePositions_HonourThePreventingResidue()
    {
        // trypsin|P must not pin the K in "AKP", because that K is not a cleavage site.
        var trypsinP = DigestionMotif.ParseDigestionMotifsFromString("K[P]|,R[P]|");
        Assert.That(DecoySequenceValidator.CleavageSitePositions("AKPCR", trypsinP).OrderBy(p => p),
            Is.EqualTo(new[] { 4 }));
    }

    #endregion
}