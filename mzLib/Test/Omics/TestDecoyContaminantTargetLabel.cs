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
    /// labels with '|', and <c>== "D"</c> reads an entrapment decoy as a target.
    /// </summary>
    [TestCase("T", false, false)]
    [TestCase("D", true, false)]
    [TestCase("C", false, false)]
    [TestCase("ET", false, true)]
    [TestCase("ED", true, true)]
    [TestCase("T|ET", false, true)]
    [TestCase("ED|D", true, true)]
    [TestCase("T|C", false, false)]
    [TestCase("", false, false)]
    [TestCase(null, false, false)]
    public void ReadingALabelFindsDecoyAndEntrapmentAnywhereInIt(string? label, bool isDecoy, bool isEntrapment)
    {
        Assert.That(DecoyContaminantTargetLabel.IsDecoy(label), Is.EqualTo(isDecoy));
        Assert.That(DecoyContaminantTargetLabel.IsEntrapment(label), Is.EqualTo(isEntrapment));
    }
}
