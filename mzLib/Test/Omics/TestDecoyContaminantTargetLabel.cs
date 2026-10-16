using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;

namespace Test.Omics;

[TestFixture]
public class TestDecoyContaminantTargetLabel
{
    [TestCase(false, false, false, "T")]
    [TestCase(true, false, false, "D")]
    [TestCase(false, true, false, "C")]
    [TestCase(false, false, true, "ET")]
    [TestCase(true, false, true, "ED")]
    [TestCase(true, true, false, "D")]
    public void TheLabelFollowsEdEtDCT(bool isDecoy, bool isContaminant, bool isEntrapment, string expected)
    {
        Assert.That(DecoyContaminantTargetLabel.For(isDecoy, isContaminant, isEntrapment), Is.EqualTo(expected));
    }

    [Test]
    public void ABioPolymerIsLabelledFromItsOwnFlags()
    {
        var entrapmentDecoy = new Protein("PEPTIDEK", "DECOY_Random_P1_f0", isDecoy: true, isEntrapment: true);
        var entrapmentTarget = new Protein("PEPTIDEK", "Random_P1_f0", isEntrapment: true);

        Assert.That(DecoyContaminantTargetLabel.For(entrapmentDecoy), Is.EqualTo("ED"));
        Assert.That(DecoyContaminantTargetLabel.For(entrapmentTarget), Is.EqualTo("ET"));
        Assert.That(DecoyContaminantTargetLabel.For(new Protein("PEPTIDEK", "P1")), Is.EqualTo("T"));
    }

    /// <summary>
    /// Readers test for a letter, never for equality: a PSM mapping to several parents joins their
    /// labels with '|', and <c>== "D"</c> reads an entrapment decoy as a target. A value that is not
    /// a label at all names no parents, so the "E" in "Excel" is not entrapment.
    /// </summary>
    [TestCase("T", false, false, false)]
    [TestCase("D", true, false, false)]
    [TestCase("C", false, true, false)]
    [TestCase("ET", false, false, true)]
    [TestCase("ED", true, false, true)]
    [TestCase("T|ET", false, false, true)]
    [TestCase("ED|D", true, false, true)]
    [TestCase("T|C", false, true, false)]
    [TestCase("C|D", true, true, false)]
    [TestCase("", false, false, false)]
    [TestCase(null, false, false, false)]
    [TestCase("Output too long for Excel", false, false, false)]
    public void ReadingALabelFindsEachLetterAnywhereInIt(string? label, bool isDecoy, bool isContaminant, bool isEntrapment)
    {
        Assert.That(DecoyContaminantTargetLabel.IsDecoy(label), Is.EqualTo(isDecoy));
        Assert.That(DecoyContaminantTargetLabel.IsContaminant(label), Is.EqualTo(isContaminant));
        Assert.That(DecoyContaminantTargetLabel.IsEntrapment(label), Is.EqualTo(isEntrapment));
    }

    /// <summary>
    /// A shared PSM's label names one parent per candidate, in the writer's order. Anything that is
    /// not a label (a blank cell, MetaMorpheus's "Output too long for Excel") names no parents.
    /// </summary>
    [TestCase("T", "T")]
    [TestCase("T|T|D", "T,T,D")]
    [TestCase(" T | ET ", "T,ET")]
    [TestCase("ED|D", "ED,D")]
    [TestCase("", "")]
    [TestCase(" ", "")]
    [TestCase(null, "")]
    [TestCase("Output too long for Excel", "")]
    [TestCase("T|", "")]
    public void ALabelNamesOneParentPerCandidate(string? label, string expected)
    {
        Assert.That(string.Join(",", DecoyContaminantTargetLabel.Parents(label)), Is.EqualTo(expected));
    }

    /// <summary>
    /// MetaMorpheus counts a shared PSM as a fraction of a decoy when it computes FDR: each
    /// candidate adds 1/n, decoy or target by its own parent (FdrAnalysisEngine.CalculateQValue,
    /// <c>groupDecoyHits / totalHits</c>). So T|T|D is a third of a decoy, not a decoy.
    /// </summary>
    [TestCase("T", 0.0)]
    [TestCase("D", 1.0)]
    [TestCase("ED", 1.0)]
    [TestCase("T|T|D", 1.0 / 3)]
    [TestCase("T|D", 0.5)]
    [TestCase("ET|ED", 0.5)]
    [TestCase("D|ED", 1.0)]
    [TestCase("T|C", 0.0)]
    [TestCase("", 0.0)]
    [TestCase(null, 0.0)]
    [TestCase("Output too long for Excel", 0.0)]
    public void ASharedPsmIsTheFractionOfADecoyThatFdrCounts(string? label, double expected)
    {
        Assert.That(DecoyContaminantTargetLabel.DecoyFraction(label), Is.EqualTo(expected).Within(1e-12));
    }

    /// <summary>
    /// An entrapment candidate counts as its share of the PSM only when its full sequence is its
    /// own. A sequence that a target or contaminant candidate also carries is that real peptide,
    /// not an ambiguity, so it counts as target. A label collapsed to one value applies to every
    /// candidate, and so does a collapsed full sequence.
    /// </summary>
    [TestCase("ET", "PEPTIDEK", 1.0)]
    [TestCase("T", "PEPTIDEK", 0.0)]
    [TestCase("ED", "PEPTIDEK", 0.0)]
    [TestCase("T|ET", "PEPTIDEK|PEPTLDEK", 0.5)]
    [TestCase("T|ET", "PEPTIDEK", 0.0)]
    [TestCase("T|ET|ET", "AAK|BBK|CCK", 2.0 / 3)]
    [TestCase("C|ET", "AAK|AAK", 0.0)]
    [TestCase("ET|ED", "AAK|BBK", 0.5)]
    [TestCase("ET", "AAK|BBK", 1.0)]
    [TestCase("D|T|T|ET", "TLQLIRK|TLQLLRK|TLQLLRK|TLQLLRK", 0.0)]
    [TestCase("D|T|ET", "TLQLIRK|TLQLLRK|TLQLIRK", 1.0 / 3)]
    [TestCase("", "PEPTIDEK", 0.0)]
    [TestCase(null, null, 0.0)]
    [TestCase("Output too long for Excel", "PEPTIDEK", 0.0)]
    public void AnEntrapmentCandidateCountsOnlyWhenItsSequenceIsItsOwn(string? label, string? fullSequence, double expected)
    {
        Assert.That(DecoyContaminantTargetLabel.EntrapmentFraction(label, fullSequence), Is.EqualTo(expected).Within(1e-12));
    }

    /// <summary>
    /// Columns that cannot be lined up candidate by candidate have no answer, and saying 0 or 1
    /// would hide that.
    /// </summary>
    [TestCase("T|ET", "AAK|BBK|CCK")]
    [TestCase("T|ET|ET", "AAK|BBK")]
    public void MisalignedColumnsHaveNoEntrapmentFraction(string label, string fullSequence)
    {
        Assert.That(DecoyContaminantTargetLabel.EntrapmentFraction(label, fullSequence), Is.NaN);
    }

    /// <summary>
    /// A full sequence that names nothing -- blank, or the placeholder MetaMorpheus writes for an
    /// over-long joined cell -- cannot say whether an entrapment candidate shares a target's
    /// peptide. It used to read as one collapsed sequence that every candidate shared, so a
    /// <c>T|ET</c> row read 0, although MetaMorpheus writes that placeholder only when the sequences
    /// differ. Where the sequences cannot change the answer, the answer stands.
    /// </summary>
    [TestCase("T|ET", "Output too long for Excel", double.NaN)]
    [TestCase("C|ET", "Output too long for Excel", double.NaN)]
    [TestCase("T|ET", "", double.NaN)]
    [TestCase("T|ET", " ", double.NaN)]
    [TestCase("T|ET", null, double.NaN)]
    [TestCase("ET|D", "Output too long for Excel", 0.5)]
    [TestCase("ET", "Output too long for Excel", 1.0)]
    [TestCase("T|D", "Output too long for Excel", 0.0)]
    public void AnUnreadableFullSequenceAnswersOnlyWhatItCan(string label, string? fullSequence, double expected)
    {
        double actual = DecoyContaminantTargetLabel.EntrapmentFraction(label, fullSequence);

        if (double.IsNaN(expected))
        {
            Assert.That(actual, Is.NaN);
        }
        else
        {
            Assert.That(actual, Is.EqualTo(expected).Within(1e-12));
        }
    }

    /// <summary>
    /// The precedence is ED, ET, D, C, T, so entrapment outranks contaminant. No case set both
    /// flags. It matters for groups, whose flags are any-member: one contaminant member and one
    /// entrapment member set both, and the group is written as entrapment.
    /// </summary>
    [TestCase(false, true, true, "ET")]
    [TestCase(true, true, true, "ED")]
    public void EntrapmentOutranksContaminant(bool isDecoy, bool isContaminant, bool isEntrapment, string expected)
    {
        Assert.That(DecoyContaminantTargetLabel.For(isDecoy, isContaminant, isEntrapment), Is.EqualTo(expected));
    }

    /// <summary>
    /// A value that is not a label names no parent, even when it contains the letter asked about.
    /// The only such case was the "E" in "Output too long for Excel", so reading decoy or
    /// contaminant from the raw string, as before 3cc95941, still passed.
    /// </summary>
    [TestCase("Decoy")]
    [TestCase("Contaminant")]
    [TestCase("TD")]
    [TestCase("T|DC")]
    public void AValueThatIsNotALabelNamesNoParentEvenWhenItContainsTheLetter(string label)
    {
        Assert.That(DecoyContaminantTargetLabel.Parents(label), Is.Empty);
        Assert.That(DecoyContaminantTargetLabel.IsDecoy(label), Is.False);
        Assert.That(DecoyContaminantTargetLabel.IsContaminant(label), Is.False);
        Assert.That(DecoyContaminantTargetLabel.IsEntrapment(label), Is.False);
    }

    /// <summary>
    /// FDR counts every candidate as 1/n, and a contaminant is a target there, so a contaminant
    /// candidate stays in the denominator. No case mixed C with D.
    /// </summary>
    [TestCase("C|D", 0.5)]
    [TestCase("T|C|D", 1.0 / 3)]
    public void AContaminantCandidateCountsInTheDecoyFractionsDenominator(string label, double expected)
    {
        Assert.That(DecoyContaminantTargetLabel.DecoyFraction(label), Is.EqualTo(expected).Within(1e-12));
    }

    /// <summary>
    /// A collapsed full sequence applies to every candidate, so the denominator is the number of
    /// candidates the label names, not the one sequence written. The only collapsed-sequence case
    /// with several parents was <c>T|ET</c>, whose numerator is 0, so any denominator passed.
    /// </summary>
    [TestCase("ET|D", "AAK", 0.5)]
    [TestCase("D|ET", "AAK", 0.5)]
    public void ACollapsedFullSequenceStillDividesByEveryCandidate(string label, string fullSequence, double expected)
    {
        Assert.That(DecoyContaminantTargetLabel.EntrapmentFraction(label, fullSequence), Is.EqualTo(expected).Within(1e-12));
    }

    /// <summary>
    /// A target or contaminant claims its sequence wherever it stands in the label. Every case put
    /// the real protein first, so a claim that only looked backwards passed.
    /// </summary>
    [TestCase("ET|T", "PEPTIDEK", 0.0)]
    [TestCase("ET|D|T", "AAK|BBK|AAK", 0.0)]
    [TestCase("ET|C", "AAK", 0.0)]
    public void ARealProteinClaimsItsSequenceWhereverItStandsInTheLabel(string label, string fullSequence, double expected)
    {
        Assert.That(DecoyContaminantTargetLabel.EntrapmentFraction(label, fullSequence), Is.EqualTo(expected).Within(1e-12));
    }
}
