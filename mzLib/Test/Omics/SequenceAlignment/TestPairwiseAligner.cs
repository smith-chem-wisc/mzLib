using NUnit.Framework;
using Omics.SequenceAlignment;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.Omics.SequenceAlignment;

/// <summary>
/// The pairwise aligner maps a residue of one sequence to the residue of another that it
/// corresponds to. The first use is joining a UniProt entry's residue to the same gene's Ensembl
/// protein, so most pairs are identical or nearly so, and the traps are the ends (signal peptides,
/// extended termini), internal indels, and residues the matrix does not name (U, X).
/// Expected scores are worked by hand from NCBI's BLOSUM62.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class TestPairwiseAligner
{
    [Test]
    public void Blosum62_IsNcbisMatrixVerbatim()
    {
        var m = SubstitutionMatrix.Blosum62;
        Assert.That(m.Name, Is.EqualTo("BLOSUM62"));
        Assert.That(string.Concat(m.Codes), Is.EqualTo("ARNDCQEGHILKMFPSTWYVBZX*"));
        Assert.That(m.Score('A', 'A'), Is.EqualTo(4));
        Assert.That(m.Score('W', 'W'), Is.EqualTo(11));
        Assert.That(m.Score('A', 'R'), Is.EqualTo(-1));
        Assert.That(m.Score('R', 'A'), Is.EqualTo(-1));
        Assert.That(m.Score('B', 'D'), Is.EqualTo(4));
        Assert.That(m.Score('*', '*'), Is.EqualTo(1));
        Assert.That(m.Score('X', 'X'), Is.EqualTo(-1));
    }

    [Test]
    public void Blosum62_ScoresAnUnnamedResidueAsX_AndIgnoresCase()
    {
        var m = SubstitutionMatrix.Blosum62;
        Assert.That(m.Score('U', 'C'), Is.EqualTo(m.Score('X', 'C')));
        Assert.That(m.Score('O', 'K'), Is.EqualTo(m.Score('X', 'K')));
        Assert.That(m.Score('a', 'a'), Is.EqualTo(4));
    }

    [Test]
    public void Blosum62_IdenticalSequencesAlignToThemselves()
    {
        Assert.That(SubstitutionMatrix.Blosum62.IdenticalSequencesAlignToThemselves, Is.True);
    }

    [Test]
    public void Align_IdenticalSequences_PairEveryResidueWithItself()
    {
        var a = new PairwiseAligner().Align("PEPTIDE", "PEPTIDE");
        for (int k = 1; k <= 7; k++)
        {
            Assert.That(a.TryMapSourceToTarget(k, out int t), Is.True);
            Assert.That(t, Is.EqualTo(k));
        }
        Assert.That(a.Score, Is.EqualTo(7 + 5 + 7 + 5 + 4 + 6 + 5));
        Assert.That(a.Identical, Is.EqualTo(7));
        Assert.That(a.AlignedSource, Is.EqualTo("PEPTIDE"));
    }

    [Test]
    public void Align_IdenticalSequencesWithAnX_GiveTheSameAnswerThroughTheTable()
    {
        // X scores -1 against itself, so the shortcut is not taken; the table must agree.
        var a = new PairwiseAligner().Align("PXPTIDE", "PXPTIDE");
        Assert.That(a.AlignedSource, Is.EqualTo("PXPTIDE"));
        Assert.That(a.AlignedTarget, Is.EqualTo("PXPTIDE"));
        Assert.That(a.Score, Is.EqualTo(7 - 1 + 7 + 5 + 4 + 6 + 5));
    }

    [Test]
    public void Align_ASubstitution_KeepsThePositionAndCountsIt()
    {
        var a = new PairwiseAligner().Align("MKTAYIAK", "MKTSYIAK");
        Assert.That(a.TryMapSourceToTarget(4, out int t), Is.True);
        Assert.That(t, Is.EqualTo(4));
        Assert.That(a.Substituted, Is.EqualTo(1));
        Assert.That(a.Identical, Is.EqualTo(7));
        Assert.That(a.Score, Is.EqualTo(5 + 5 + 5 + 1 + 7 + 4 + 4 + 5));
    }

    [Test]
    public void Align_AnInternalDeletion_LeavesTheDeletedResiduesFacingAGap()
    {
        var a = new PairwiseAligner().Align("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKISFVKSHFSRQ");
        Assert.That(a.AlignedTarget, Is.EqualTo("MKTAYIAK---ISFVKSHFSRQ"));
        foreach (int deleted in new[] { 9, 10, 11 })
        {
            Assert.That(a.TryMapSourceToTarget(deleted, out int none), Is.False);
            Assert.That(none, Is.EqualTo(0));
        }
        Assert.That(a.TryMapSourceToTarget(12, out int afterGap), Is.True);
        Assert.That(afterGap, Is.EqualTo(9));
        Assert.That(a.TryMapTargetToSource(9, out int back), Is.True);
        Assert.That(back, Is.EqualTo(12));
        Assert.That(a.SourceResiduesFacingGap, Is.EqualTo(3));
        // 39 before the gap, 55 after it, and 11 + 3 x 1 for the gap.
        Assert.That(a.Score, Is.EqualTo(39 + 55 - 14));
    }

    [Test]
    public void Align_AnExtraNTerminus_IsFreeByDefault_AndChargedWhenAsked()
    {
        var free = new PairwiseAligner().Align("MLLLPEPTIDE", "PEPTIDE");
        Assert.That(free.TryMapSourceToTarget(5, out int t), Is.True);
        Assert.That(t, Is.EqualTo(1));
        Assert.That(free.Score, Is.EqualTo(39));

        var charged = new PairwiseAligner(freeEndGaps: false).Align("MLLLPEPTIDE", "PEPTIDE");
        Assert.That(charged.AlignedTarget, Is.EqualTo("----PEPTIDE"));
        Assert.That(charged.Score, Is.EqualTo(39 - (11 + 4)));
    }

    [Test]
    public void Align_AnExtraCTerminusOnTheTarget_IsFreeByDefault()
    {
        var a = new PairwiseAligner().Align("PEPTIDE", "PEPTIDEKKKK");
        Assert.That(a.AlignedSource, Is.EqualTo("PEPTIDE----"));
        Assert.That(a.TargetResiduesFacingGap, Is.EqualTo(4));
        Assert.That(a.TryMapTargetToSource(8, out _), Is.False);
        Assert.That(a.Score, Is.EqualTo(39));
    }

    [Test]
    public void Align_ATie_IsBrokenTheSameWayEveryTime()
    {
        // "A" pairs equally well with either A of "AA"; the documented rule ends at the full corner.
        var a = new PairwiseAligner().Align("AA", "A");
        Assert.That(a.TryMapSourceToTarget(1, out _), Is.False);
        Assert.That(a.TryMapSourceToTarget(2, out int t), Is.True);
        Assert.That(t, Is.EqualTo(1));
    }

    [Test]
    public void Align_AnEmptySequence_FacesEveryResidueWithAGap()
    {
        var a = new PairwiseAligner().Align("", "PEP");
        Assert.That(a.AlignedSource, Is.EqualTo("---"));
        Assert.That(a.Score, Is.EqualTo(0));
    }

    [Test]
    public void Align_TheAlignedStringsAreTheInputsAndRescoreToTheScore()
    {
        var random = new Random(1381);
        const string residues = "ACDEFGHIKLMNPQRSTVWY";
        foreach (bool freeEnds in new[] { true, false })
        {
            var aligner = new PairwiseAligner(freeEndGaps: freeEnds);
            for (int trial = 0; trial < 200; trial++)
            {
                string s = RandomSequence(random, residues, random.Next(0, 40));
                string t = Mutate(random, s, residues);
                var a = aligner.Align(s, t);
                Assert.That(a.AlignedSource.Replace("-", ""), Is.EqualTo(s));
                Assert.That(a.AlignedTarget.Replace("-", ""), Is.EqualTo(t));
                Assert.That(Rescore(a.AlignedSource, a.AlignedTarget, aligner), Is.EqualTo(a.Score),
                    $"{a.AlignedSource} / {a.AlignedTarget}");
            }
        }
    }

    [Test]
    public void Align_RefusesATableLargerThanMaxCells()
    {
        var aligner = new PairwiseAligner(maxCells: 10);
        var ex = Assert.Throws<ArgumentException>(() => aligner.Align("AAAA", "AAAC"));
        Assert.That(ex!.Message, Does.Contain("MaxCells"));
        // Identical sequences need no table, so they are not refused.
        Assert.That(aligner.Align("AAAAAAAA", "AAAAAAAA").Identical, Is.EqualTo(8));
    }

    [Test]
    public void Map_RefusesAPositionOutsideTheSequence()
    {
        var a = new PairwiseAligner().Align("PEPTIDE", "PEPTIDE");
        Assert.Throws<ArgumentOutOfRangeException>(() => a.TryMapSourceToTarget(0, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.TryMapTargetToSource(8, out _));
    }

    [Test]
    public void Id_NamesTheMethodAndEveryParameter()
    {
        Assert.That(new PairwiseAligner().Id, Is.EqualTo("global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"));
        Assert.That(new PairwiseAligner(gapOpen: 10, gapExtend: 2, freeEndGaps: false).Id,
            Is.EqualTo("global-affine;BLOSUM62;open=10;extend=2;end-gaps=charged"));
        Assert.That(new PairwiseAligner().Align("PEP", "PEP").AlignerId, Is.EqualTo(new PairwiseAligner().Id));
    }

    [Test]
    public void Constructor_RefusesNegativeGapCosts()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new PairwiseAligner(gapOpen: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PairwiseAligner(gapExtend: 0));
    }

    [Test]
    public void Constructor_RefusesMaxCellsLargerThanAnArrayCanHold()
    {
        // The traceback is one array of MaxCells bytes; past Array.MaxLength it cannot be built.
        Assert.Throws<ArgumentOutOfRangeException>(() => new PairwiseAligner(maxCells: (long)Array.MaxLength + 1));
        Assert.That(new PairwiseAligner(maxCells: Array.MaxLength).MaxCells, Is.EqualTo(Array.MaxLength));
    }

    [TestCase("PEP-TIDE", "PEPTIDE")]
    [TestCase("PEPTIDE", "PEP.TIDE")]
    [TestCase("PEP-TIDE", "PEP-TIDE")]
    public void Align_RefusesAGappedRow(string source, string target)
    {
        // A gap character is not a residue; scored as X it would align silently.
        var ex = Assert.Throws<ArgumentException>(() => new PairwiseAligner().Align(source, target));
        Assert.That(ex!.Message, Does.Contain("gap"));
    }

    private static IEnumerable<TestCaseData> MalformedMatrices()
    {
        yield return new TestCaseData("  A  R\nA  4 -1\nR  0  5\n", "not symmetric").SetName("Parse_RefusesAnAsymmetricMatrix");
        yield return new TestCaseData("  A  R\nA  4 -1\n", "1 rows, expected 2").SetName("Parse_RefusesAMissingRow");
        yield return new TestCaseData("  A  R\nA  4\nR -1  5\n", "line 2: 1 scores, expected 2").SetName("Parse_RefusesAShortRow");
        yield return new TestCaseData("  A  R\nA  4 x\nR -1  5\n", "line 2: score 'x' is not an integer").SetName("Parse_RefusesANonIntegerScore");
        yield return new TestCaseData("  A  R\nR  4 -1\nA -1  5\n", "line 2: row 'R', expected 'A'").SetName("Parse_RefusesARowOutOfOrder");
        yield return new TestCaseData("  A  A\nA  4  4\nA  4  4\n", "line 1: residue code 'A' is not a single new letter").SetName("Parse_RefusesARepeatedCode");
    }

    [TestCaseSource(nameof(MalformedMatrices))]
    public void Parse_RefusesRatherThanGuesses(string text, string message)
    {
        var ex = Assert.Throws<InvalidDataException>(() => SubstitutionMatrix.Parse("TEST", new StringReader(text)));
        Assert.That(ex!.Message, Does.Contain(message));
    }

    [Test]
    public void Parse_AMatrixWithoutX_RefusesAnUnnamedResidue()
    {
        var m = SubstitutionMatrix.Parse("AR", new StringReader("# comment\n  A  R\nA  4 -1\nR -1  5\n"));
        Assert.That(m.Score('a', 'R'), Is.EqualTo(-1));
        Assert.Throws<ArgumentException>(() => m.Score('A', 'U'));
    }

    private static string RandomSequence(Random random, string residues, int length)
        => new string(Enumerable.Range(0, length).Select(_ => residues[random.Next(residues.Length)]).ToArray());

    // A related sequence: substitutions, an internal indel, and changed ends.
    private static string Mutate(Random random, string s, string residues)
    {
        var chars = s.ToList();
        for (int k = 0; k < chars.Count; k++)
        {
            if (random.NextDouble() < 0.1) chars[k] = residues[random.Next(residues.Length)];
        }
        if (chars.Count > 4 && random.NextDouble() < 0.5) chars.RemoveRange(random.Next(chars.Count - 3), random.Next(1, 4));
        if (random.NextDouble() < 0.5) chars.InsertRange(random.Next(chars.Count + 1), RandomSequence(random, residues, random.Next(1, 5)));
        if (random.NextDouble() < 0.3) chars.InsertRange(0, RandomSequence(random, residues, random.Next(1, 6)));
        if (random.NextDouble() < 0.3) chars.AddRange(RandomSequence(random, residues, random.Next(1, 6)));
        return string.Concat(chars);
    }

    // Scores an alignment from its two strings, independently of the aligner's tables.
    private static int Rescore(string alignedSource, string alignedTarget, PairwiseAligner aligner)
    {
        int score = 0, k = 0, length = alignedSource.Length;
        while (k < length)
        {
            if (alignedSource[k] != '-' && alignedTarget[k] != '-')
            {
                score += aligner.Matrix.Score(alignedSource[k], alignedTarget[k]);
                k++;
                continue;
            }
            bool sourceGap = alignedTarget[k] == '-';
            int start = k;
            while (k < length && (sourceGap ? alignedTarget[k] == '-' && alignedSource[k] != '-'
                                            : alignedSource[k] == '-' && alignedTarget[k] != '-'))
            {
                k++;
            }
            bool atEnd = start == 0 || k == length;
            if (!(aligner.FreeEndGaps && atEnd))
            {
                score -= aligner.GapOpen + (k - start) * aligner.GapExtend;
            }
        }
        return score;
    }
}
