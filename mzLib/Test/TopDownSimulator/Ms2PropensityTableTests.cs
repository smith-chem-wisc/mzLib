using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Ms2;

namespace Test.TopDownSimulator;

[TestFixture]
public class Ms2PropensityTableTests
{
    private const string CanonicalResidues = "ACDEFGHIKLMNPQRSTVWY";

    [Test]
    public void EmbeddedTableLoadsWithTheCompletePairGrid()
    {
        var table = ResiduePairCleavagePropensityTable.Default;

        Assert.That(table.PairCount, Is.EqualTo(400),
            "The published sheet is a complete ordered 20 x 20 residue grid.");
        Assert.That(table.Residues, Is.EquivalentTo(CanonicalResidues.ToCharArray()));
    }

    [Test]
    public void EveryCanonicalPairHasAPropensityInUnitInterval()
    {
        var table = ResiduePairCleavagePropensityTable.Default;
        var missing = new List<string>();

        foreach (char n in CanonicalResidues)
        {
            foreach (char c in CanonicalResidues)
            {
                if (!table.TryGetPropensity(n, c, out double p))
                {
                    missing.Add($"{n}|{c}");
                    continue;
                }

                Assert.That(p, Is.GreaterThan(0).And.LessThanOrEqualTo(1),
                    $"Propensity for {n}|{c} should be a fraction in (0, 1].");
            }
        }

        Assert.That(missing, Is.Empty,
            "The denatured column covers all 400 pairs with no zeros, which is why it is the one modelled.");
    }

    [Test]
    public void PublishedValuesAreCarriedVerbatim()
    {
        var table = ResiduePairCleavagePropensityTable.Default;

        // Spot checks straight out of the "Data" worksheet, Propensity (dTD) column.
        Assert.That(table.GetPropensity('C', 'W'), Is.EqualTo(0.11101499423298732).Within(1e-15));
        Assert.That(table.GetPropensity('P', 'H'), Is.EqualTo(0.0072365615672109768).Within(1e-15));
        Assert.That(table.GetPropensity('A', 'A'), Is.EqualTo(0.11976339430609857).Within(1e-15));
    }

    [Test]
    public void GetPropensityThrowsForUnknownPair()
    {
        var table = ResiduePairCleavagePropensityTable.Default;
        Assert.Throws<KeyNotFoundException>(() => table.GetPropensity('X', 'A'));
    }

    [Test]
    public void ModelDefaultsToTheIdentityTransformOfThePublishedValue()
    {
        var model = new ResiduePairCleavagePropensityModel();

        Assert.That(model.PropensityScale, Is.EqualTo(1.0));
        Assert.That(model.MaximumCleavageProbability, Is.EqualTo(1.0));

        // The bond in "CWA" at index 0 is C on the N-terminal side, W on the C-terminal side:
        // the ~11% bond from the user-facing worked example.
        Assert.That(model.GetCleavageProbability("CWA", 0),
            Is.EqualTo(0.11101499423298732).Within(1e-15));
    }

    [Test]
    public void PropensityScaleIsAppliedAndClamped()
    {
        double published = ResiduePairCleavagePropensityTable.Default.GetPropensity('C', 'W');

        var doubled = new ResiduePairCleavagePropensityModel(propensityScale: 2.0);
        Assert.That(doubled.GetCleavageProbability("CWA", 0), Is.EqualTo(published * 2).Within(1e-15));

        var clamped = new ResiduePairCleavagePropensityModel(propensityScale: 100.0);
        Assert.That(clamped.GetCleavageProbability("CWA", 0), Is.EqualTo(1.0).Within(1e-15));

        var capped = new ResiduePairCleavagePropensityModel(propensityScale: 100.0, maximumCleavageProbability: 0.5);
        Assert.That(capped.GetCleavageProbability("CWA", 0), Is.EqualTo(0.5).Within(1e-15));
    }

    [Test]
    public void UnknownResidueFallsBackToMissingPairPropensity()
    {
        var zero = new ResiduePairCleavagePropensityModel();
        Assert.That(zero.GetCleavageProbability("XA", 0), Is.EqualTo(0.0));

        var fallback = new ResiduePairCleavagePropensityModel(missingPairPropensity: 0.07);
        Assert.That(fallback.GetCleavageProbability("XA", 0), Is.EqualTo(0.07).Within(1e-15));
    }

    [Test]
    public void UnmeasuredPairInACustomTableRoutesThroughTheFallback()
    {
        var table = ResiduePairCleavagePropensityTable.Parse(new[]
        {
            "NTerminalResidue\tCTerminalResidue\tDenaturedTopDownPropensity",
            "A\tA\tNA",
            "A\tC\t0.5"
        });

        var model = new ResiduePairCleavagePropensityModel(table, missingPairPropensity: 0.07);

        Assert.That(model.GetCleavageProbability("AA", 0), Is.EqualTo(0.07).Within(1e-15));
        Assert.That(model.GetCleavageProbability("AC", 0), Is.EqualTo(0.5).Within(1e-15));
    }

    [Test]
    public void BondIndexIsValidated()
    {
        var model = new ResiduePairCleavagePropensityModel();
        Assert.Throws<ArgumentOutOfRangeException>(() => model.GetCleavageProbability("PEPTIDE", -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => model.GetCleavageProbability("PEPTIDE", 6));
        Assert.DoesNotThrow(() => model.GetCleavageProbability("PEPTIDE", 5));
    }

    [Test]
    public void ParseRejectsValuesOutsideUnitInterval()
    {
        string[] percentages =
        {
            "NTerminalResidue\tCTerminalResidue\tDenaturedTopDownPropensity",
            "A\tA\t11.1"
        };

        Assert.Throws<FormatException>(() => ResiduePairCleavagePropensityTable.Parse(percentages),
            "A percentage-scaled table must be rejected rather than silently normalized.");
    }

    [Test]
    public void ParseRequiresThePropensityColumn()
    {
        string[] lines =
        {
            "NTerminalResidue\tCTerminalResidue\tSomeOtherColumn",
            "A\tA\t0.25"
        };

        Assert.Throws<FormatException>(() => ResiduePairCleavagePropensityTable.Parse(lines));
    }

    [Test]
    public void ParseHonoursCommentsAndNaMarkers()
    {
        string[] lines =
        {
            "# a comment",
            "",
            "NTerminalResidue\tCTerminalResidue\tDenaturedTopDownPropensity",
            "A\tA\tNA",
            "# another comment",
            "A\tC\t0.5"
        };

        var table = ResiduePairCleavagePropensityTable.Parse(lines);
        Assert.That(table.PairCount, Is.EqualTo(2));
        Assert.That(table.Contains('A', 'A'), Is.False, "NA means not measured.");
        Assert.That(table.GetPropensity('A', 'C'), Is.EqualTo(0.5));
    }
}
