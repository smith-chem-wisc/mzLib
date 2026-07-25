using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using TopDownSimulator.Ms2;

namespace Test.TopDownSimulator;

[TestFixture]
public class Ms2CascadeFragmentationTests
{
    /// <summary>Fixed per-bond probabilities so the cascade arithmetic can be hand-checked.</summary>
    private sealed class FixedBondCleavageModel : IBondCleavageModel
    {
        private readonly double[] _probabilities;

        public FixedBondCleavageModel(params double[] probabilities) => _probabilities = probabilities;

        public double GetCleavageProbability(string baseSequence, int zeroBasedBondIndex) =>
            _probabilities[zeroBasedBondIndex];
    }

    /// <summary>
    /// Seven residues, six bonds, and deliberately proline-free: mzLib suppresses the z-dot ion when
    /// the residue C-terminal to the cleavage is proline (the N-Cα bond there does not break under
    /// ETD), which would make the complementary-ion bookkeeping in these tests a special case.
    /// <see cref="ProlineSuppressedZDotIntensityReturnsToThePrecursor"/> covers that separately.
    /// </summary>
    private const string Sequence = "ELVISLK";

    private static PeptideWithSetModifications Peptide(string sequence = Sequence) =>
        new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());

    private static PropensityCascadeFragmentationModel HandComputedModel(
        PropensityCascadeOptions options = null) =>
        new PropensityCascadeFragmentationModel(
            new FixedBondCleavageModel(0.1, 0.2, 0.0, 0.5, 0.25, 0.4),
            options ?? new PropensityCascadeOptions { StartingIonCount = 1000.0 });

    [Test]
    public void CascadeMatchesTheHandComputedRecursion()
    {
        var cascade = HandComputedModel().ComputeCascade(Sequence);

        // R_0 = 1000; at each bond cleaved_i = p_i * R_i and R_(i+1) = R_i * (1 - p_i).
        double[] expectedRemaining = { 1000, 900, 720, 720, 360, 270 };
        double[] expectedCleaved = { 100, 180, 0, 360, 90, 108 };

        Assert.That(cascade.BondCount, Is.EqualTo(Sequence.Length - 1));
        Assert.That(cascade.RemainingBeforeBond, Is.EqualTo(expectedRemaining).Within(1e-12));
        Assert.That(cascade.CleavedIntensity, Is.EqualTo(expectedCleaved).Within(1e-12));
        Assert.That(cascade.SurvivingPrecursorIntensity, Is.EqualTo(162.0).Within(1e-12));
    }

    [Test]
    public void CascadeSatisfiesTheRecursionTermByTerm()
    {
        var cascade = HandComputedModel().ComputeCascade(Sequence);

        Assert.That(cascade.RemainingBeforeBond[0], Is.EqualTo(cascade.StartingIonCount).Within(1e-12));

        for (int i = 0; i < cascade.BondCount; i++)
        {
            double p = cascade.CleavageProbabilities[i];
            double r = cascade.RemainingBeforeBond[i];

            Assert.That(cascade.CleavedIntensity[i], Is.EqualTo(p * r).Within(1e-12),
                $"cleaved_{i} should equal p_{i} * R_{i}.");

            double next = i + 1 < cascade.BondCount
                ? cascade.RemainingBeforeBond[i + 1]
                : cascade.SurvivingPrecursorIntensity;

            Assert.That(next, Is.EqualTo(r * (1 - p)).Within(1e-12),
                $"R_{i + 1} should equal R_{i} * (1 - p_{i}).");
        }
    }

    [Test]
    public void RemainingIntensityDecaysMonotonicallyAlongTheBackbone()
    {
        // Real propensities on a long sequence: every bond has p > 0 in the denatured column, so the
        // decay is strict.
        const string longSequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK";
        var model = new PropensityCascadeFragmentationModel(
            new ResiduePairCleavagePropensityModel(),
            new PropensityCascadeOptions { StartingIonCount = 1e6 });

        var cascade = model.ComputeCascade(longSequence);

        Assert.That(cascade.BondCount, Is.EqualTo(longSequence.Length - 1));
        Assert.That(cascade.CleavageProbabilities.All(p => p > 0), Is.True);
        Assert.That(cascade.RemainingBeforeBond, Is.Ordered.Descending);

        double previous = cascade.RemainingBeforeBond[^1];
        Assert.That(cascade.SurvivingPrecursorIntensity, Is.LessThan(previous));
        Assert.That(cascade.SurvivingPrecursorIntensity, Is.GreaterThan(0));
    }

    [Test]
    public void RemainingIntensityIsFlatAcrossAZeroProbabilityBond()
    {
        var cascade = HandComputedModel().ComputeCascade(Sequence);

        // Bond index 2 has p = 0, so R_3 == R_2 - non-increasing, not strictly decreasing.
        Assert.That(cascade.RemainingBeforeBond[3], Is.EqualTo(cascade.RemainingBeforeBond[2]).Within(1e-12));
        for (int i = 1; i < cascade.BondCount; i++)
            Assert.That(cascade.RemainingBeforeBond[i], Is.LessThanOrEqualTo(cascade.RemainingBeforeBond[i - 1]));
    }

    [Test]
    public void TotalAssignedIntensityPlusSurvivorEqualsStartingIonCount()
    {
        var model = HandComputedModel();
        var result = model.Fragment(Peptide());

        double assigned = result.Fragments.Sum(f => f.Intensity);
        Assert.That(assigned + result.SurvivingPrecursorIntensity,
            Is.EqualTo(result.StartingIonCount).Within(1e-9),
            "The cascade must conserve the starting ion count.");
    }

    [Test]
    public void ConservationHoldsForEveryIonAssignmentAndChargeSet()
    {
        var assignments = new[]
        {
            FragmentIonAssignment.NTerminalOnly,
            FragmentIonAssignment.CTerminalOnly,
            FragmentIonAssignment.Both
        };

        foreach (var assignment in assignments)
        {
            foreach (int[] charges in new[] { new[] { 1 }, new[] { 1, 2, 3 } })
            {
                foreach (bool envelopes in new[] { false, true })
                {
                    var model = HandComputedModel(new PropensityCascadeOptions
                    {
                        StartingIonCount = 1000.0,
                        IonAssignment = assignment,
                        FragmentChargeStates = charges,
                        UseIsotopeEnvelopes = envelopes
                    });

                    var result = model.Fragment(Peptide());
                    double assigned = result.Fragments.Sum(f => f.Intensity);

                    Assert.That(assigned + result.SurvivingPrecursorIntensity,
                        Is.EqualTo(1000.0).Within(1e-6),
                        $"Conservation failed for {assignment}, charges [{string.Join(',', charges)}], envelopes={envelopes}.");
                }
            }
        }
    }

    [Test]
    public void EvenSplitGivesEachTerminusHalfOfEachCleavageEvent()
    {
        var result = HandComputedModel().Fragment(Peptide());
        var byBond = result.Fragments.GroupBy(f => f.OneBasedBondNumber).ToDictionary(g => g.Key, g => g.ToList());

        double[] expectedCleaved = { 100, 180, 0, 360, 90, 108 };
        for (int bond = 1; bond <= 6; bond++)
        {
            if (expectedCleaved[bond - 1] == 0)
            {
                Assert.That(byBond.ContainsKey(bond), Is.False, "A zero-probability bond should emit no ions.");
                continue;
            }

            var ions = byBond[bond];
            Assert.That(ions, Has.Count.EqualTo(2), $"Bond {bond} should yield an N- and a C-terminal ion.");

            var nIon = ions.Single(i => i.Product.Terminus == FragmentationTerminus.N);
            var cIon = ions.Single(i => i.Product.Terminus == FragmentationTerminus.C);

            Assert.That(nIon.Intensity, Is.EqualTo(expectedCleaved[bond - 1] / 2).Within(1e-9));
            Assert.That(cIon.Intensity, Is.EqualTo(expectedCleaved[bond - 1] / 2).Within(1e-9));
        }
    }

    [Test]
    public void UnevenSplitIsHonoured()
    {
        var model = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            NTerminalIonFraction = 0.8
        });

        var bondOne = model.Fragment(Peptide()).Fragments.Where(f => f.OneBasedBondNumber == 1).ToList();
        Assert.That(bondOne.Single(f => f.Product.Terminus == FragmentationTerminus.N).Intensity,
            Is.EqualTo(80.0).Within(1e-9));
        Assert.That(bondOne.Single(f => f.Product.Terminus == FragmentationTerminus.C).Intensity,
            Is.EqualTo(20.0).Within(1e-9));
    }

    [Test]
    public void DefaultsToEtdCAndZDotSeries()
    {
        var options = new PropensityCascadeOptions();
        Assert.That(options.DissociationType, Is.EqualTo(DissociationType.ETD));

        var (n, c) = options.ResolveProductTypes();
        Assert.That(n, Is.EqualTo(ProductType.c));
        Assert.That(c, Is.EqualTo(ProductType.zDot));

        var result = HandComputedModel().Fragment(Peptide());
        Assert.That(result.Fragments.Select(f => f.Product.ProductType).Distinct(),
            Is.EquivalentTo(new[] { ProductType.c, ProductType.zDot }));

        var (hcdN, hcdC) = new PropensityCascadeOptions { DissociationType = DissociationType.HCD }.ResolveProductTypes();
        Assert.That(hcdN, Is.EqualTo(ProductType.b));
        Assert.That(hcdC, Is.EqualTo(ProductType.y));
    }

    [Test]
    public void ExplicitProductTypeOverridesWin()
    {
        var model = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            DissociationType = DissociationType.EThcD,
            NTerminalProductType = ProductType.b,
            CTerminalProductType = ProductType.y
        });

        var result = model.Fragment(Peptide());
        Assert.That(result.Fragments.Select(f => f.Product.ProductType).Distinct(),
            Is.EquivalentTo(new[] { ProductType.b, ProductType.y }));
        Assert.That(result.Fragments.Sum(f => f.Intensity) + result.SurvivingPrecursorIntensity,
            Is.EqualTo(1000.0).Within(1e-9));
    }

    [Test]
    public void BondNumberMapsToComplementaryFragmentNumbers()
    {
        var result = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            IonAssignment = FragmentIonAssignment.Both
        }).Fragment(Peptide());

        foreach (var fragment in result.Fragments)
        {
            int expected = fragment.Product.Terminus == FragmentationTerminus.N
                ? fragment.OneBasedBondNumber
                : Sequence.Length - fragment.OneBasedBondNumber;

            Assert.That(fragment.Product.FragmentNumber, Is.EqualTo(expected),
                $"{fragment.Product.Annotation} is on the wrong side of bond {fragment.OneBasedBondNumber}.");
        }
    }

    [Test]
    public void MzComesFromMzLibFragmentChemistry()
    {
        var peptide = Peptide();
        var expectedProducts = new List<Product>();
        peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, expectedProducts);

        var result = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            FragmentChargeStates = new[] { 1, 2 }
        }).Fragment(peptide);

        foreach (var fragment in result.Fragments)
        {
            Assert.That(fragment.Mz, Is.EqualTo(fragment.Product.NeutralMass.ToMz(fragment.Charge)).Within(1e-12));

            var match = expectedProducts.SingleOrDefault(p =>
                p.ProductType == fragment.Product.ProductType &&
                p.FragmentNumber == fragment.Product.FragmentNumber);
            Assert.That(match, Is.Not.Null, $"{fragment.Product.Annotation} is not an mzLib product of this peptide.");
            Assert.That(fragment.Product.NeutralMass, Is.EqualTo(match.NeutralMass).Within(1e-12));
        }
    }

    [Test]
    public void ChargeStatesSplitIntensityEvenly()
    {
        var model = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            IonAssignment = FragmentIonAssignment.NTerminalOnly,
            FragmentChargeStates = new[] { 1, 2, 3, 4 }
        });

        var bondOne = model.Fragment(Peptide()).Fragments.Where(f => f.OneBasedBondNumber == 1).ToList();
        Assert.That(bondOne, Has.Count.EqualTo(4));
        Assert.That(bondOne.Select(f => f.Charge), Is.EquivalentTo(new[] { 1, 2, 3, 4 }));
        foreach (var ion in bondOne)
            Assert.That(ion.Intensity, Is.EqualTo(25.0).Within(1e-9));
    }

    [Test]
    public void IsotopeEnvelopeModeSpreadsOneFragmentOverManyPeaks()
    {
        var single = HandComputedModel().Fragment(Peptide());
        var enveloped = HandComputedModel(new PropensityCascadeOptions
        {
            StartingIonCount = 1000.0,
            UseIsotopeEnvelopes = true
        }).Fragment(Peptide());

        Assert.That(single.Fragments.All(f => f.IsotopologueIndex == 0), Is.True);
        Assert.That(enveloped.Fragments.Count, Is.GreaterThan(single.Fragments.Count));
        Assert.That(enveloped.Fragments.Any(f => f.IsotopologueIndex > 0), Is.True);

        // Per-fragment totals must still agree with the single-peak model.
        var singleByKey = single.Fragments.ToDictionary(f => (f.Product.Annotation, f.Charge), f => f.Intensity);
        var envelopedByKey = enveloped.Fragments
            .GroupBy(f => (f.Product.Annotation, f.Charge))
            .ToDictionary(g => g.Key, g => g.Sum(f => f.Intensity));

        Assert.That(envelopedByKey.Keys, Is.EquivalentTo(singleByKey.Keys));
        foreach (var key in singleByKey.Keys)
            Assert.That(envelopedByKey[key], Is.EqualTo(singleByKey[key]).Within(1e-6));
    }

    [Test]
    public void TheDefaultModelUsesThePublishedTable()
    {
        var model = new PropensityCascadeFragmentationModel();
        Assert.That(model.BondCleavageModel, Is.InstanceOf<ResiduePairCleavagePropensityModel>());
        Assert.That(model.Options.StartingIonCount, Is.EqualTo(1e6));

        var cascade = model.ComputeCascade("CWA");
        Assert.That(cascade.CleavageProbabilities[0], Is.EqualTo(0.11101499423298732).Within(1e-15));

        // The worked example: an 11% bond takes 11% of the ions and hands the rest downstream.
        Assert.That(cascade.CleavedIntensity[0], Is.EqualTo(1e6 * 0.11101499423298732).Within(1e-6));
        Assert.That(cascade.RemainingBeforeBond[1],
            Is.EqualTo(1e6 * (1 - 0.11101499423298732)).Within(1e-6));
    }

    [Test]
    public void ProlineSuppressedZDotIntensityReturnsToThePrecursor()
    {
        // PEPTIDE has proline at position 3, so mzLib emits no zDot5 - the ion complementary to
        // bond 2. That cleavage still consumes intensity from the cascade, so the C-terminal half of
        // it has nowhere to land and must be returned to the surviving precursor rather than lost.
        var result = HandComputedModel().Fragment(Peptide("PEPTIDE"));

        var bondTwo = result.Fragments.Where(f => f.OneBasedBondNumber == 2).ToList();
        Assert.That(bondTwo, Has.Count.EqualTo(1));
        Assert.That(bondTwo[0].Product.Terminus, Is.EqualTo(FragmentationTerminus.N));
        Assert.That(bondTwo[0].Intensity, Is.EqualTo(90.0).Within(1e-9));

        // Cascade survivor is 162; the orphaned 90 from bond 2 is added back.
        Assert.That(result.SurvivingPrecursorIntensity, Is.EqualTo(162.0 + 90.0).Within(1e-9));
        Assert.That(result.Fragments.Sum(f => f.Intensity) + result.SurvivingPrecursorIntensity,
            Is.EqualTo(1000.0).Within(1e-9));
    }

    [Test]
    public void SingleResidueSequenceHasNoBondsAndSurvivesIntact()
    {
        var model = new PropensityCascadeFragmentationModel(
            new FixedBondCleavageModel(), new PropensityCascadeOptions { StartingIonCount = 500 });

        var cascade = model.ComputeCascade("A");
        Assert.That(cascade.BondCount, Is.EqualTo(0));
        Assert.That(cascade.SurvivingPrecursorIntensity, Is.EqualTo(500));
    }

    [Test]
    public void InvalidOptionsAreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new PropensityCascadeFragmentationModel(
            options: new PropensityCascadeOptions { StartingIonCount = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PropensityCascadeFragmentationModel(
            options: new PropensityCascadeOptions { NTerminalIonFraction = 1.5 }));
        Assert.Throws<ArgumentException>(() => new PropensityCascadeFragmentationModel(
            options: new PropensityCascadeOptions { FragmentChargeStates = Array.Empty<int>() }));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PropensityCascadeFragmentationModel(
            options: new PropensityCascadeOptions { DissociationType = DissociationType.PQD }));
    }

    [Test]
    public void ProbabilitiesOutsideUnitIntervalAreRejected()
    {
        var model = new PropensityCascadeFragmentationModel(
            new FixedBondCleavageModel(1.5), new PropensityCascadeOptions { StartingIonCount = 10 });

        Assert.Throws<InvalidOperationException>(() => model.ComputeCascade("AA"));
    }
}
