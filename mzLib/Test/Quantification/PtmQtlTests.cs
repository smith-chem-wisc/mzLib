using NUnit.Framework;
using Quantification.PtmQtl;
using Statistics;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Quantification;

/// <summary>Hand-computed cases for intensity-based site occupancy (DEF-OCC-INT) and PTM pairs (types P and A).</summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PtmQtlTests
{
    private const string Phos = "Common Biological:Phosphorylation on S";
    private static PeptidoformObservation Obs(string run, string full, int start, double intensity, string protein = "P1") =>
        new(run, full, protein, start, start + MzLibUtil.ClassExtensions.GetBaseSequenceFromFullSequence(full).Length - 1, intensity);

    private static SiteRunOccupancy One(IReadOnlyList<SiteRunOccupancy> occ, string run, int position) =>
        occ.Single(o => o.Run == run && o.Site.Position == position);

    [Test]
    public void OccupancyIsModifiedOverCoveringIntensity()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", $"PEPS[{Phos}]K", 10, 30),
            Obs("r1", "PEPSK", 10, 70),
            Obs("r1", "AAPEPSKR", 8, 100),        // a missed-cleavage form covering the same residue, unmodified
        });
        var s = One(occ, "r1", 13);
        Assert.That(s.Site.Residue, Is.EqualTo('S'));
        Assert.That(s.Site.Key, Is.EqualTo($"P1:S13:{Phos}"));
        Assert.That(s.State, Is.EqualTo(OccupancyState.Quantified));
        Assert.That(s.Fraction, Is.EqualTo(30.0 / 200).Within(1e-15));
        Assert.That(s.UnmodifiedQuantified, Is.True);
    }

    [Test]
    public void SeenOnlyModifiedIsOccupancyOneWithTheCeilingFlag()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[] { Obs("r1", $"PEPS[{Phos}]K", 10, 30) });
        var s = One(occ, "r1", 13);
        Assert.That(s.Fraction, Is.EqualTo(1.0));
        Assert.That(s.UnmodifiedQuantified, Is.False);
    }

    [Test]
    public void TheFourStatesAreKeptApart()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("floor", $"PEPS[{Phos}]K", 10, double.NaN), Obs("floor", "PEPSK", 10, 50),
            Obs("countOnly", $"PEPS[{Phos}]K", 10, double.NaN),
            Obs("notDetected", "PEPSK", 10, 40),
            Obs("quant", $"PEPS[{Phos}]K", 10, 5), Obs("quant", "PEPSK", 10, 45),
        });
        Assert.That(One(occ, "floor", 13).State, Is.EqualTo(OccupancyState.Floor));
        Assert.That(One(occ, "floor", 13).Fraction, Is.NaN, "a floor is not a zero");
        Assert.That(One(occ, "countOnly", 13).State, Is.EqualTo(OccupancyState.CountOnly));
        Assert.That(One(occ, "notDetected", 13).State, Is.EqualTo(OccupancyState.NotDetected));
        Assert.That(One(occ, "notDetected", 13).Fraction, Is.NaN, "not detected is unknown, not 0");
        Assert.That(One(occ, "quant", 13).Fraction, Is.EqualTo(0.1).Within(1e-15));
    }

    [Test]
    public void ExclusionsAndFilterFollowDefOccPsms()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", "PEPM[Common Variable:Oxidation on M]C[Common Fixed:Carbamidomethyl on C]K", 10, 10),
            Obs("r1", "[Common Artifact:Ammonia loss on C]CPEPK", 20, 10),                   // peptide N-terminus, not protein
            Obs("r1", "[UniProt:N-acetylserine on S]SPEPK", 2, 10),                          // protein N-terminus after Met removal
            Obs("r1", $"PEPS[{Phos}]K", 30, 10),
        }, mod => !mod.StartsWith("Common Artifact"));
        var keys = occ.Select(o => o.Site.Key).Distinct().ToList();
        Assert.That(keys, Is.EquivalentTo(new[] { "P1:S2:UniProt:N-acetylserine on S", $"P1:S33:{Phos}" }));
    }

    [Test]
    public void TwoModificationsAtOnePositionShareTheDenominator()
    {
        const string Glc = "Common Biological:HexNAc on S";
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", $"PEPS[{Phos}]K", 10, 20), Obs("r1", $"PEPS[{Glc}]K", 10, 30), Obs("r1", "PEPSK", 10, 50),
        });
        var at13 = occ.Where(o => o.Site.Position == 13).ToList();
        Assert.That(at13, Has.Count.EqualTo(2));
        Assert.That(at13.Sum(o => o.Fraction), Is.EqualTo(0.5).Within(1e-15));
    }

    [Test]
    public void InconsistentCoordinatesAreRefused()
    {
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(new[]
            { new PeptidoformObservation("r1", "PEPSK", "P1", 10, 20, 5) }));
    }

    [Test]
    public void PhysicalPairsCarryCoOccupancy()
    {
        const string Acet = "Common Biological:Acetylation on K";
        var obs = new[]
        {
            Obs("r1", $"PEPS[{Phos}]K[{Acet}]R", 10, 25),   // both marks on one molecule
            Obs("r1", $"PEPS[{Phos}]KR", 10, 25),
            Obs("r1", "PEPSKR", 10, 50),
            Obs("r2", $"PEPS[{Phos}]K[{Acet}]R", 10, double.NaN),
        };
        var p = PtmPairEngine.Physical(obs).Single();
        Assert.That(p.ResultType, Is.EqualTo(PairResultType.P));
        Assert.That(string.CompareOrdinal(p.SiteA.Key, p.SiteB.Key), Is.LessThan(0), "written once, a < b");
        Assert.That(p.N, Is.EqualTo(2), "identified together in two runs");
        Assert.That(p.Statistic, Is.EqualTo(0.25).Within(1e-15), "25 of 100 covering both positions, in the one quantified run");
    }

    [Test]
    public void CoVaryingPairsAreSpearmanAcrossRunsWithOverlapExcluded()
    {
        var obs = new List<PeptidoformObservation>();
        var rng = new Random(1);
        for (int r = 0; r < 10; r++)
        {
            string run = $"r{r}";
            double occ = 0.1 + 0.08 * r;                       // rises across runs
            obs.Add(Obs(run, $"PEPS[{Phos}]K", 10, 100 * occ, "P1"));
            obs.Add(Obs(run, "PEPSK", 10, 100 * (1 - occ), "P1"));
            obs.Add(Obs(run, $"GGS[{Phos}]R", 40, 100 * occ * 0.5, "P2"));   // follows P1, on another protein
            obs.Add(Obs(run, "GGSR", 40, 100 * (1 - occ * 0.5), "P2"));
            double noise = rng.NextDouble();
            obs.Add(Obs(run, $"TTS[{Phos}]AS[{Phos}]K", 60, 50 * noise, "P1"));  // two sites on one peptide: overlapping
            obs.Add(Obs(run, $"TTS[{Phos}]ASK", 60, 30, "P1"));
            obs.Add(Obs(run, "TTSASK", 60, 20, "P1"));
        }
        var occupancy = SiteOccupancyCalculator.Calculate(obs);
        var pairs = PtmPairEngine.CoVarying(occupancy, obs);
        var cross = pairs.Single(p => p.SiteA.ProteinAccession != p.SiteB.ProteinAccession && p.SiteA.Position == 13);
        Assert.That(cross.Statistic, Is.EqualTo(1).Within(1e-12));
        Assert.That(cross.FdrFamily, Is.EqualTo("A:inter"));
        Assert.That(cross.Q, Is.Not.NaN);
        var overlapping = pairs.Single(p => p.SiteA.Position == 62 && p.SiteB.Position == 64
                                           || p.SiteA.Position == 64 && p.SiteB.Position == 62);
        Assert.That(overlapping.Overlapping, Is.True);
        Assert.That(overlapping.Q, Is.NaN, "overlapping pairs stay out of the family");
        Assert.That(cross.SpearmanMethod, Is.EqualTo(SpearmanPValueMethod.Asymptotic));
    }

    [Test]
    public void SitesQuantifiedInTooFewRunsDoNotEnterTypeA()
    {
        var obs = new List<PeptidoformObservation>();
        for (int r = 0; r < 10; r++)
        {
            obs.Add(Obs($"r{r}", $"PEPS[{Phos}]K", 10, 10 + r, "P1"));
            obs.Add(Obs($"r{r}", "PEPSK", 10, 50, "P1"));
            if (r < 6) obs.Add(Obs($"r{r}", $"GGS[{Phos}]R", 40, 5 + r, "P2"));   // 6 of 10 runs < 70%
        }
        var pairs = PtmPairEngine.CoVarying(SiteOccupancyCalculator.Calculate(obs), obs);
        Assert.That(pairs, Is.Empty);
    }
}
