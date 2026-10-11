using NUnit.Framework;
using Quantification.PtmQtl;
using StatisticalModels;
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
        Assert.That(One(occ, "notDetected", 13).Fraction, Is.NaN, "a covered-zero is a bound, not a measured 0");
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
    public void PeptideNTerminalModificationsAndNonMetStartsAreNotProteinNTermini()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", "[Common Artifact:Ammonia loss on C]CPEPK", 1, 10, "P1"),                   // peptide N-terminal mod at residue 1
            new PeptidoformObservation("r1", "[UniProt:N-acetylserine on S]SPEPK", "P2", 2, 6, 10, 'K'), // residue 2 after K1: no Met removal
            new PeptidoformObservation("r1", "[UniProt:N-acetylserine on S]SPEPK", "P3", 2, 6, 10, 'M'), // residue 2 after Met removal
            Obs("r1", "[UniProt:N-acetylserine on S]SPEPK", 2, 10, "P4"),                            // unknown: the caller guarantees Met
            Obs("r1", "[Renamed:Ammonia loss on C]CPEPK", 1, 10, "P5"),                              // unknown category: matched on IdWithMotif
            Obs("r1", "[Renamed:Unregistered on C]CPEPK", 1, 10, "P6"),                              // unknown to the registry: protein N-terminal
        });
        Assert.That(occ.Select(o => o.Site.Key).Distinct(),
            Is.EquivalentTo(new[] { "P3:S2:UniProt:N-acetylserine on S", "P4:S2:UniProt:N-acetylserine on S", "P6:C1:Renamed:Unregistered on C" }));
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

        // With no unmodified form, each modification still sees the other: neither is a ceiling.
        var noUnmodified = SiteOccupancyCalculator.Calculate(new[] { Obs("r1", $"PEPS[{Phos}]K", 10, 20), Obs("r1", $"PEPS[{Glc}]K", 10, 30) });
        Assert.That(noUnmodified.Single(o => o.Site.Modification == Phos).Fraction, Is.EqualTo(0.4).Within(1e-15));
        Assert.That(noUnmodified.All(o => o.UnmodifiedQuantified), Is.True);
    }

    [Test]
    public void InconsistentCoordinatesAreRefused()
    {
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(new[]
            { new PeptidoformObservation("r1", "PEPSK", "P1", 10, 20, 5) }));
    }

    [Test]
    public void MalformedObservationsAreRefused()
    {
        void Refused(PeptidoformObservation o) =>
            Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(new[] { o }), o?.ToString() ?? "null row");
        Assert.Throws<ArgumentNullException>(() => SiteOccupancyCalculator.Calculate(null!));
        Refused(null!);
        Refused(new PeptidoformObservation("", "PEPSK", "P1", 10, 14, 5));
        Refused(new PeptidoformObservation("r1", "PEPSK", "", 10, 14, 5));
        Refused(new PeptidoformObservation("r1", "PEPSK", "P1", 0, 4, 5));
        Refused(Obs("r1", "PEPSK", 10, -1));
        Refused(Obs("r1", "PEPSK", 10, double.PositiveInfinity));
    }

    [Test]
    public void FractionsOfOneSampleAreSummedBeforeTheRatio()
    {
        var obs = new[]
        {
            Obs("f1", $"PEPS[{Phos}]K", 10, 30), Obs("f1", "PEPSK", 10, 70),
            Obs("f2", $"PEPS[{Phos}]K", 10, 10), Obs("f2", "PEPSK", 10, 90),
            Obs("g1", "PEPSK", 10, 40),
        };
        var map = new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A", ["g1"] = "B" };
        var occ = SiteOccupancyCalculator.Calculate(obs, map);

        var a = One(occ, "A", 13);
        Assert.That(a.State, Is.EqualTo(OccupancyState.Quantified));
        Assert.That(a.Fraction, Is.EqualTo(40.0 / 200).Within(1e-15), "(30 + 10) / (100 + 100), not the mean of 0.3 and 0.1");
        Assert.That(a.CoveringIntensity, Is.EqualTo(200));
        Assert.That(One(occ, "B", 13).State, Is.EqualTo(OccupancyState.NotDetected));
        Assert.That(occ.Select(o => o.Run).Distinct(), Is.EquivalentTo(new[] { "A", "B" }));
    }

    [Test]
    public void AFormUnquantifiedInOneFractionAndCoveredInAnotherIsAFloorForTheSample()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("f1", $"PEPS[{Phos}]K", 10, double.NaN),
            Obs("f2", "PEPSK", 10, 50),
        }, new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A" });
        var a = One(occ, "A", 13);
        Assert.That(a.State, Is.EqualTo(OccupancyState.Floor));
        Assert.That(a.Fraction, Is.NaN);

        var quantifiedSomewhere = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("f1", $"PEPS[{Phos}]K", 10, double.NaN), Obs("f2", $"PEPS[{Phos}]K", 10, 20), Obs("f2", "PEPSK", 10, 60),
        }, new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A" });
        Assert.That(One(quantifiedSomewhere, "A", 13).Fraction, Is.EqualTo(20.0 / 80).Within(1e-15),
            "an unquantified fraction adds nothing; the quantified one carries the sum");
    }

    [Test]
    public void APeptidoformInTwoFractionsIsNotARepeatedRowButARepeatWithinOneFractionStillIs()
    {
        var map = new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A" };
        Assert.DoesNotThrow(() => SiteOccupancyCalculator.Calculate(new[] { Obs("f1", "PEPSK", 10, 1), Obs("f2", "PEPSK", 10, 1) }, map));
        Assert.Throws<ArgumentException>(() =>
            SiteOccupancyCalculator.Calculate(new[] { Obs("f1", "PEPSK", 10, 1), Obs("f1", "PEPSK", 10, 1) }, map));
    }

    [Test]
    public void CombinedObservationsFeedSameMoleculePairsPerSample()
    {
        const string GG = "Trypsin Digested:GG (Ubiquitination Site) on K";
        var obs = new[]
        {
            Obs("f1", $"PEPS[{Phos}]K[{GG}]R", 10, 20), Obs("f1", "PEPSKR", 10, 80),
            Obs("f2", $"PEPS[{Phos}]K[{GG}]R", 10, 10), Obs("f2", "PEPSKR", 10, 90),
            Obs("g1", $"PEPS[{Phos}]K[{GG}]R", 10, double.NaN),
        };
        var map = new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A", ["g1"] = "B" };
        var combined = SiteOccupancyCalculator.CombineObservations(obs, map);
        Assert.That(combined.Count, Is.EqualTo(3));
        Assert.That(combined.Single(o => o.Run == "A" && o.FullSequence == "PEPSKR").Intensity, Is.EqualTo(170));
        Assert.That(combined.Single(o => o.Run == "B").Intensity, Is.NaN);

        var p = PtmPairEngine.Physical(combined).Single();
        Assert.That(p.N, Is.EqualTo(2), "two samples, not three runs");
        Assert.That(p.Statistic, Is.EqualTo(30.0 / 200).Within(1e-15), "co-occupancy of sample A only; B is unquantified");
    }

    [Test]
    public void ARunWithoutASampleIsRefused()
    {
        var obs = new[] { Obs("f1", "PEPSK", 10, 1), Obs("f2", "PEPSK", 10, 1) };
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(obs, new Dictionary<string, string> { ["f1"] = "A" }));
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(obs, new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "" }));
        Assert.Throws<ArgumentNullException>(() => SiteOccupancyCalculator.Calculate(obs, (IReadOnlyDictionary<string, string>)null!));
    }

    [Test]
    public void StoredCellsCombineAsTheObservationsWouldAndKeepTheReportedFraction()
    {
        var obs = new[]
        {
            Obs("f1", $"PEPS[{Phos}]K", 10, 30), Obs("f1", "PEPSK", 10, 70),
            Obs("f2", $"PEPS[{Phos}]K", 10, double.NaN), Obs("f2", "PEPSK", 10, 50),
            Obs("f3", "PEPSK", 10, 30),
        };
        var map = new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A", ["f3"] = "A" };
        var fromObservations = One(SiteOccupancyCalculator.Calculate(obs, map), "A", 13);
        var fromCells = One(SiteOccupancyCalculator.CombineRuns(SiteOccupancyCalculator.Calculate(obs), map), "A", 13);
        Assert.That(fromCells.State, Is.EqualTo(fromObservations.State));
        Assert.That(fromCells.Fraction, Is.EqualTo(fromObservations.Fraction).Within(1e-15));
        Assert.That(fromCells.Fraction, Is.EqualTo(30.0 / 180).Within(1e-15), "the Floor and NotDetected runs add covering only");
        Assert.That(fromCells.UnmodifiedQuantified, Is.True);

        // A stored cell's exact fraction is used, not its rounded intensities.
        var site = new ModificationSite("P1", 13, 'S', Phos);
        var stored = new[]
        {
            new SiteRunOccupancy(site, "f1", OccupancyState.Quantified, 30, 100, true) { ReportedFraction = 0.25 },
            new SiteRunOccupancy(site, "f2", OccupancyState.Quantified, 1, 300, true) { ReportedFraction = 0.5 },
        };
        var c = SiteOccupancyCalculator.CombineRuns(stored, new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A" }).Single();
        Assert.That(c.Fraction, Is.EqualTo((0.25 * 100 + 0.5 * 300) / 400).Within(1e-15));
    }

    [Test]
    public void StoredCellStatesCombineByPrecedence()
    {
        var site = new ModificationSite("P1", 13, 'S', Phos);
        SiteRunOccupancy Cell(string run, OccupancyState state, double covering) => new(site, run, state, 0, covering, false);
        var map = new Dictionary<string, string> { ["a"] = "S", ["b"] = "S" };
        OccupancyState Combined(params SiteRunOccupancy[] cells) => SiteOccupancyCalculator.CombineRuns(cells, map).Single().State;

        Assert.That(Combined(Cell("a", OccupancyState.CountOnly, double.NaN), Cell("b", OccupancyState.NotDetected, 40)),
            Is.EqualTo(OccupancyState.Floor), "identified in one fraction, covered in another");
        Assert.That(Combined(Cell("a", OccupancyState.CountOnly, double.NaN), Cell("b", OccupancyState.CountOnly, double.NaN)),
            Is.EqualTo(OccupancyState.CountOnly));
        Assert.That(Combined(Cell("a", OccupancyState.NotDetected, 10), Cell("b", OccupancyState.NotDetected, 40)),
            Is.EqualTo(OccupancyState.NotDetected));
    }

    [Test]
    public void BadStoredCellsAreRefused()
    {
        var site = new ModificationSite("P1", 13, 'S', Phos);
        var map = new Dictionary<string, string> { ["a"] = "S" };
        var ok = new SiteRunOccupancy(site, "a", OccupancyState.Quantified, 1, 10, true);
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.CombineRuns(new[] { ok, ok }, map), "repeated (site, run)");
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.CombineRuns(new[] { ok with { CoveringIntensity = double.NaN } }, map));
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.CombineRuns(new[] { ok with { Run = "z" } }, map), "run without a sample");
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.CombineRuns(new SiteRunOccupancy[] { null! }, map));
    }

    [Test]
    public void CTerminalModificationsAreNotSites()
    {
        // Protein length is not supplied, so a C-terminal modification cannot be placed; the residue it follows is not its site.
        var occ = SiteOccupancyCalculator.Calculate(new[] { Obs("r1", "PEPSK-[Common Biological:Amidation on K]", 10, 30), Obs("r1", "PEPSK", 10, 70) });
        Assert.That(occ, Is.Empty);
    }

    [Test]
    public void PsmLevelRowsAreRefused()
    {
        // Two PSMs of one peptidoform in one run would count its intensity twice (occupancy 60/90, not 30/60).
        var psms = new[] { Obs("r1", $"PEPS[{Phos}]K", 10, 30), Obs("r1", $"PEPS[{Phos}]K", 10, 30), Obs("r1", "PEPSK", 10, 30) };
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.Calculate(psms));
        Assert.Throws<ArgumentException>(() => PtmPairEngine.Physical(psms), "type P validates as occupancy does");

        // The same peptidoform in another run, or on another protein, is a row of its own.
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", $"PEPS[{Phos}]K", 10, 30), Obs("r2", $"PEPS[{Phos}]K", 10, 30), Obs("r1", $"PEPS[{Phos}]K", 10, 30, "P2"),
        });
        Assert.That(occ, Has.Count.EqualTo(3));
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
        Assert.That(p.StatisticN, Is.EqualTo(1), "quantified in one of them");
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
        Assert.That(cross.StatisticN, Is.EqualTo(cross.N));
        Assert.That(cross.FdrFamily, Is.EqualTo("A:inter"));
        Assert.That(cross.Q, Is.Not.NaN);
        var overlapping = pairs.Single(p => p.SiteA.Position == 62 && p.SiteB.Position == 64
                                           || p.SiteA.Position == 64 && p.SiteB.Position == 62);
        Assert.That(overlapping.Overlapping, Is.True);
        Assert.That(overlapping.Q, Is.NaN, "overlapping pairs stay out of the family");
        Assert.That(cross.SpearmanMethod, Is.EqualTo(SpearmanPValueMethod.Asymptotic));
    }

    [Test]
    public void SharedPeptideSitesOnDifferentProteinsAreOverlapping()
    {
        // One peptide mapped to an isoform pair and to a third protein at another position (one row per protein):
        // all three sites take their occupancy from the same peptidoforms, so ρ = 1 by construction. A distinct
        // peptide on P2 is a genuine inter-protein pair.
        var obs = new List<PeptidoformObservation>();
        var rng = new Random(3);
        for (int r = 0; r < 10; r++)
        {
            string run = $"r{r}";
            double occ = 0.1 + 0.08 * r;
            foreach (var (protein, start) in new[] { ("P1", 10), ("P1-2", 10), ("P9", 50) })
            {
                obs.Add(Obs(run, $"PEPS[{Phos}]K", start, 100 * occ, protein));
                obs.Add(Obs(run, "PEPSK", start, 100 * (1 - occ), protein));
            }
            double other = rng.NextDouble();
            obs.Add(Obs(run, $"GGS[{Phos}]R", 40, 100 * other, "P2"));
            obs.Add(Obs(run, "GGSR", 40, 100 * (1 - other), "P2"));
        }
        var pairs = PtmPairEngine.CoVarying(SiteOccupancyCalculator.Calculate(obs), obs);
        Assert.That(pairs, Has.Count.EqualTo(6));
        var shared = pairs.Where(p => p.SiteA.ProteinAccession != "P2" && p.SiteB.ProteinAccession != "P2").ToList();
        Assert.That(shared, Has.Count.EqualTo(3));
        Assert.That(shared.All(p => p.Overlapping && double.IsNaN(p.Q)), Is.True, "shared-peptide pairs stay out of A:inter");
        var distinct = pairs.Except(shared).ToList();
        Assert.That(distinct.All(p => !p.Overlapping && p.FdrFamily == "A:inter" && !double.IsNaN(p.Q)), Is.True);
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

    [Test]
    public void CeilingCellsStayOutOfTypeAByDefault()
    {
        // Two sites opposed in five runs, and both seen only modified (ceilings at 1) in three low-load runs.
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 5; r++)
        {
            occupancy.Add(Cell($"r{r}", 0.1 * (r + 1), site: s));
            occupancy.Add(Cell($"r{r}", 0.1 * (5 - r), site: t));
        }
        for (int r = 5; r < 8; r++)
        {
            occupancy.Add(Cell($"r{r}", 1, unmodifiedQuantified: false, site: s));
            occupancy.Add(Cell($"r{r}", 1, unmodifiedQuantified: false, site: t));
        }
        var pair = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), 0.5).Single();
        Assert.That(pair.N, Is.EqualTo(5));
        Assert.That(pair.Statistic, Is.EqualTo(-1).Within(1e-12));

        // Kept, the shared ceilings tie at rank 7 in both sites: ρ = 20 / 40 from detection alone.
        var withCeiling = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), 0.5, excludeCeiling: false).Single();
        Assert.That(withCeiling.N, Is.EqualTo(8));
        Assert.That(withCeiling.Statistic, Is.EqualTo(0.5).Within(1e-12));
    }

    [Test]
    public void RunLevelRemovalTakesOutABatchShift()
    {
        // Twelve runs in two batches; every site reads 0.2 higher in the second. Within a batch, s and t vary
        // nearly independently (within-batch ρ = 1/35), so almost all their correlation comes from the batch. Four
        // background sites sit at the batch's centre, so each run's median is its batch level.
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var background = Enumerable.Range(0, 4).Select(i => new ModificationSite($"P{3 + i}", 7, 'S', Phos)).ToList();
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 12; r++)
        {
            double batch = r < 6 ? 0.1 : 0.3;
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseS[r % 6], site: s));
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseT[r % 6], site: t));
            foreach (var b in background) occupancy.Add(Cell($"r{r:D2}", batch + 0.035, site: b));
        }
        var plain = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>()).Single(p => p.SiteA == s && p.SiteB == t);
        Assert.That(plain.Statistic, Is.GreaterThan(0.7), "the batch alone correlates the two sites");
        Assert.That(plain.RunLevelRemoved, Is.False);

        var adjusted = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true)
            .Single(p => p.SiteA == s && p.SiteB == t);
        Assert.That(adjusted.RunLevelRemoved, Is.True);
        Assert.That(Math.Abs(adjusted.Statistic), Is.LessThan(0.1), "given the batch, the pair is unrelated");
        Assert.That(adjusted.N, Is.EqualTo(12));

        // The run level is each run's median over the six sites: the background value, batch + 0.035.
        var level = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.035).ToArray();
        var vs = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseS[r % 6]).ToArray();
        var vt = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6]).ToArray();
        var expected = SpearmanCorrelation.PartialCorrelate(vs, vt, level);
        Assert.That(adjusted.Statistic, Is.EqualTo(expected.Rho).Within(1e-12));
        Assert.That(adjusted.PValue, Is.EqualTo(expected.PValue).Within(1e-12));

        // Below MinRunsForRunLevel shared runs the partial test is not made: t with n − 3 df is too liberal there.
        var eight = occupancy.Where(o => string.CompareOrdinal(o.Run, "r08") < 0).ToList();
        var few = PtmPairEngine.CoVarying(eight, Array.Empty<PeptidoformObservation>(), removeRunLevel: true)
            .Single(p => p.SiteA == s && p.SiteB == t);
        Assert.That(few.N, Is.EqualTo(8));
        Assert.That(few.SpearmanMethod, Is.EqualTo(SpearmanPValueMethod.NotEstimable));
        Assert.That(few.Statistic, Is.NaN);
        Assert.That(few.Q, Is.NaN, "left out of the adjustment");
        Assert.That(PtmPairEngine.CoVarying(eight, Array.Empty<PeptidoformObservation>()).Single(p => p.SiteA == s && p.SiteB == t).PValue,
            Is.Not.NaN, "plain Spearman still tests it");
    }

    [Test]
    public void GroupedRunLevelTakesOutAShiftInOneChemistry()
    {
        // Two chemistries in twelve runs. Every "GG" site reads 0.2 higher in the second six runs; the "Phospho"
        // sites do not move with it. One run level for all sites is half one chemistry and half the other, so only a
        // level per chemistry removes the GG shift.
        const string Gg = "Trypsin Digested:GG (Ubiquitination Site) on K";
        var s = new ModificationSite("P1", 13, 'K', Gg);
        var t = new ModificationSite("P2", 43, 'K', Gg);
        var ggBackground = Enumerable.Range(0, 4).Select(i => new ModificationSite($"P{3 + i}", 7, 'K', Gg)).ToList();
        var phospho = Enumerable.Range(0, 7).Select(i => new ModificationSite($"Q{i}", 9, 'S', Phos)).ToList();
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        double[] phosphoLevel = { 0.30, 0.10, 0.25, 0.15, 0.20, 0.35, 0.12, 0.28, 0.18, 0.33, 0.22, 0.14 };
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 12; r++)
        {
            double batch = r < 6 ? 0.1 : 0.3;
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseS[r % 6], site: s));
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseT[r % 6], site: t));
            foreach (var b in ggBackground) occupancy.Add(Cell($"r{r:D2}", batch + 0.035, site: b));
            for (int i = 0; i < phospho.Count; i++) occupancy.Add(Cell($"r{r:D2}", phosphoLevel[r] + 0.001 * i, site: phospho[i]));
        }
        PtmPair Pair(IReadOnlyList<PtmPair> pairs) => pairs.Single(p => p.SiteA == s && p.SiteB == t);
        string Chemistry(ModificationSite site) => site.Modification.Contains("GG") ? "GG" : "phospho";

        var plain = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>()));
        var oneLevel = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true));
        var perChemistry = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true,
            runLevelGroup: Chemistry));
        Assert.That(plain.Statistic, Is.GreaterThan(0.7));
        Assert.That(Math.Abs(perChemistry.Statistic), Is.LessThan(0.1), "given the GG level, the two GG sites are unrelated");
        Assert.That(Math.Abs(perChemistry.Statistic), Is.LessThan(Math.Abs(oneLevel.Statistic)));

        // The GG level leaves the pair out: the median of the four background sites, batch + 0.035.
        var level = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.035).ToArray();
        var vs = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseS[r % 6]).ToArray();
        var vt = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6]).ToArray();
        Assert.That(perChemistry.Statistic, Is.EqualTo(SpearmanCorrelation.PartialCorrelate(vs, vt, level).Rho).Within(1e-12));

        // A cross-chemistry pair is correlated given both levels: two covariates.
        var cross = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true, runLevelGroup: Chemistry)
            .Single(p => p.SiteA == s && p.SiteB == phospho[0]);
        var phosphoOthers = Enumerable.Range(0, 12).Select(r => phosphoLevel[r] + 0.001 * 3.5).ToArray();   // median of sites 1-6
        var ggWithT = Enumerable.Range(0, 12).Select(r => new[] { (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6] }
            .Concat(Enumerable.Repeat((r < 6 ? 0.1 : 0.3) + 0.035, 4)).OrderBy(v => v).ElementAt(2)).ToArray();
        var vp = Enumerable.Range(0, 12).Select(r => phosphoLevel[r]).ToArray();
        var expected = SpearmanCorrelation.PartialCorrelate(vs, vp, new IReadOnlyList<double>[] { ggWithT, phosphoOthers });
        Assert.That(cross.Statistic, Is.EqualTo(expected.Rho).Within(1e-12));

        // A GG site seen in only two runs does not enter pairs, so it does not move the GG level either.
        var sparse = new ModificationSite("P9", 3, 'K', Gg);
        var withSparse = occupancy.Concat(new[] { Cell("r00", 0.9, site: sparse), Cell("r01", 0.95, site: sparse) }).ToList();
        var unchanged = Pair(PtmPairEngine.CoVarying(withSparse, Array.Empty<PeptidoformObservation>(), removeRunLevel: true,
            runLevelGroup: Chemistry));
        Assert.That(unchanged.Statistic, Is.EqualTo(perChemistry.Statistic).Within(1e-15));
    }

    [Test]
    public void ARunWithTooFewOtherSitesInTheGroupHasNoLevel()
    {
        // The group holds the pair and two more sites: leaving the pair out leaves 2, below MinSitesForRunLevel.
        var sites = Enumerable.Range(0, 4).Select(i => new ModificationSite($"P{i}", 5, 'S', Phos)).ToList();
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 12; r++)
            for (int i = 0; i < sites.Count; i++)
                occupancy.Add(Cell($"r{r:D2}", 0.1 + 0.01 * ((r * (i + 2)) % 7), site: sites[i]));
        var pair = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true).First();
        Assert.That(pair.SpearmanMethod, Is.EqualTo(SpearmanPValueMethod.NotEstimable));
        Assert.That(pair.N, Is.EqualTo(0), "no run has a level");
    }

    [Test]
    public void StoredOccupancyUsesTheReportedFractionAndFeedsTypeA()
    {
        // As MetaMorpheus writes it: the fraction exact, the intensities rounded to 4 significant figures.
        var s = new ModificationSite("P1", 13, 'S', "Phosphorylation on S");
        var stored = new SiteRunOccupancy(s, "r1", OccupancyState.Quantified, 1234, 5679, true) { ReportedFraction = 0.2173 };
        Assert.That(stored.Fraction, Is.EqualTo(0.2173));
        var floor = new SiteRunOccupancy(s, "r1", OccupancyState.Floor, 0, 5679, true) { ReportedFraction = 0 };
        Assert.That(floor.Fraction, Is.NaN, "a floor is censored, never a value");

        var t = new ModificationSite("P2", 40, 'S', "Phosphorylation on S");
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 8; r++)
        {
            // Rounded intensities whose ratio would order the runs differently from the reported fractions.
            occupancy.Add(new SiteRunOccupancy(s, $"r{r}", OccupancyState.Quantified, 1000, 2000, true) { ReportedFraction = 0.1 + 0.05 * r });
            occupancy.Add(new SiteRunOccupancy(t, $"r{r}", OccupancyState.Quantified, 1000 - r, 2000, true) { ReportedFraction = 0.3 + 0.02 * r });
        }
        var pair = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>()).Single();
        Assert.That(pair.Statistic, Is.EqualTo(1).Within(1e-12));
        Assert.That(pair.Overlapping, Is.False);
    }

    [Test]
    public void ACoveredZeroNeedsOnlyAnIdentificationAndCellsCountPeptidoforms()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("idOnly", "PEPSK", 10, double.NaN),                 // covered by an identification, nothing quantified
            Obs("quant", $"PEPS[{Phos}]K", 10, 5), Obs("quant", "PEPSK", 10, 45), Obs("quant", "AAPEPSKR", 8, 50),
            Obs("onlyModified", $"PEPS[{Phos}]K", 10, double.NaN),
            Obs("bothForms", $"PEPS[{Phos}]K", 10, double.NaN), Obs("bothForms", "PEPSK", 10, double.NaN),
            Obs("elsewhere", "GGGGK", 30, 20),                       // a run that never covers the site
        });
        var zero = One(occ, "idOnly", 13);
        Assert.That(zero.State, Is.EqualTo(OccupancyState.NotDetected));
        Assert.That(zero.CoveringPeptidoforms, Is.EqualTo(1));
        Assert.That(zero.ModifiedPeptidoforms, Is.EqualTo(0));
        Assert.That(One(occ, "quant", 13).CoveringPeptidoforms, Is.EqualTo(3));
        Assert.That(One(occ, "quant", 13).ModifiedPeptidoforms, Is.EqualTo(1));
        Assert.That(One(occ, "onlyModified", 13).State, Is.EqualTo(OccupancyState.CountOnly));
        Assert.That(One(occ, "onlyModified", 13).CoveringPeptidoforms, Is.EqualTo(1));
        Assert.That(One(occ, "onlyModified", 13).ModifiedPeptidoforms, Is.EqualTo(1));
        Assert.That(One(occ, "bothForms", 13).State, Is.EqualTo(OccupancyState.CountOnly));
        Assert.That(One(occ, "bothForms", 13).CoveringPeptidoforms, Is.EqualTo(2));
        Assert.That(occ.Where(o => o.Run == "elsewhere"), Is.Empty, "a site no peptidoform covered has no row");
    }

    private static SiteRunOccupancy Counted(SiteRunOccupancy cell, int covering, int modified) =>
        cell with { CoveringPeptidoforms = covering, ModifiedPeptidoforms = modified };

    [Test]
    public void BoundsEnterTypeAAsTiedRanks()
    {
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var occupancy = new List<SiteRunOccupancy>();
        double[] valuesS = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 }, valuesT = { 0.15, 0.1, 0.35, 0.3, 0.6, 0.45 };
        for (int r = 0; r < 6; r++)
        {
            occupancy.Add(Cell($"r{r}", valuesS[r], site: s));
            occupancy.Add(Cell($"r{r}", valuesT[r], site: t));
        }
        occupancy.Add(Cell("r6", double.NaN, OccupancyState.NotDetected, site: s));
        occupancy.Add(Cell("r6", 1, unmodifiedQuantified: false, site: t));                          // ceiling
        occupancy.Add(Cell("r7", double.NaN, OccupancyState.Floor, site: s));
        occupancy.Add(Cell("r7", double.NaN, OccupancyState.NotDetected, site: t));
        occupancy.Add(Cell("r8", 1, unmodifiedQuantified: false, site: s));
        occupancy.Add(Counted(Cell("r8", double.NaN, OccupancyState.CountOnly, site: t), 2, 1));     // both forms: no rank
        occupancy.Add(Counted(Cell("r9", double.NaN, OccupancyState.CountOnly, site: s), 1, 1));     // only modified: high
        occupancy.Add(Cell("r9", double.NaN, OccupancyState.Floor, site: t));

        // Without bounds, each site has a value in 6 of 10 runs: below 70 %, so no pair.
        Assert.That(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>()), Is.Empty);

        // With them, low bounds tie below every fraction and high bounds above; any such values give the same ranks.
        var pair = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), includeBounds: true).Single();
        double lo = -5, hi = 9;
        var vs = valuesS.Concat(new[] { lo, lo, hi, hi }).ToArray();
        var vt = valuesT.Concat(new[] { hi, lo, double.NaN, lo }).ToArray();
        var expected = SpearmanCorrelation.Correlate(vs, vt);
        Assert.That(pair.N, Is.EqualTo(9), "r8 is left out: t's count-only cell there has no rank");
        Assert.That(pair.Statistic, Is.EqualTo(expected.Rho).Within(1e-12));
        Assert.That(pair.PValue, Is.EqualTo(expected.PValue).Within(1e-12));
        Assert.That(pair.BoundsIncluded, Is.True);
        Assert.That(pair.TestStatus, Is.EqualTo(PairTestStatus.Tested));

        var explicitly = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), excludeCeiling: false, includeBounds: true).Single();
        Assert.That(explicitly.Statistic, Is.EqualTo(pair.Statistic));
        Assert.Throws<ArgumentException>(() => PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(),
            excludeCeiling: true, includeBounds: true));
        Assert.Throws<ArgumentException>(() => PtmPairEngine.CoVarying(
            new[] { Cell("r0", 0.1, (OccupancyState)99, site: s) }, Array.Empty<PeptidoformObservation>()));
    }

    [Test]
    public void APairWhoseSitesAreBothMostlyBoundsIsNotTested()
    {
        // Twenty runs. s and t have a value in the same two runs and are covered-zeros elsewhere; u varies everywhere.
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var u = new ModificationSite("P3", 7, 'S', Phos);
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 20; r++)
        {
            occupancy.Add(r < 2 ? Cell($"r{r:D2}", 0.1 + 0.1 * r, site: s) : Cell($"r{r:D2}", double.NaN, OccupancyState.NotDetected, site: s));
            occupancy.Add(r < 2 ? Cell($"r{r:D2}", 0.2 + 0.1 * r, site: t) : Cell($"r{r:D2}", double.NaN, OccupancyState.NotDetected, site: t));
            occupancy.Add(Cell($"r{r:D2}", 0.01 * ((7 * r) % 20 + 1), site: u));
        }
        var pairs = PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), includeBounds: true);
        var st = pairs.Single(p => p.SiteA == s && p.SiteB == t);
        Assert.That(SpearmanCorrelation.Correlate(
            Enumerable.Range(0, 20).Select(r => r < 2 ? 1.0 + r : 0).ToArray(),
            Enumerable.Range(0, 20).Select(r => r < 2 ? 1.0 + r : 0).ToArray()).PValue, Is.LessThan(1e-10),
            "the t test alone would call this pair");
        Assert.That(st.TestStatus, Is.EqualTo(PairTestStatus.BothSitesMostlyBounds));
        Assert.That(st.N, Is.EqualTo(20));
        Assert.That(st.Statistic, Is.NaN);
        Assert.That(st.PValue, Is.NaN);
        Assert.That(st.Q, Is.NaN, "left out of the adjustment");

        // One site with values in most runs is enough: the test is right or conservative there.
        var su = pairs.Single(p => p.SiteA == s && p.SiteB == u);
        Assert.That(su.TestStatus, Is.EqualTo(PairTestStatus.Tested));
        Assert.That(su.Q, Is.Not.NaN);
    }

    [Test]
    public void ARunCovariateIsPartialledOut()
    {
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var background = Enumerable.Range(0, 4).Select(i => new ModificationSite($"P{3 + i}", 7, 'S', Phos)).ToList();
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        double[] psms = { 9000, 12000, 8000, 15000, 11000, 7000, 16000, 10000, 13000, 9500, 14000, 8500 };
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 12; r++)
        {
            double batch = r < 6 ? 0.1 : 0.3;
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseS[r % 6] + psms[r] * 1e-6, site: s));
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseT[r % 6] + psms[r] * 2e-6, site: t));
            foreach (var b in background) occupancy.Add(Cell($"r{r:D2}", batch + 0.035, site: b));
        }
        var depth = Enumerable.Range(0, 12).ToDictionary(r => $"r{r:D2}", r => psms[r]);
        var vs = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseS[r % 6] + psms[r] * 1e-6).ToArray();
        var vt = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6] + psms[r] * 2e-6).ToArray();
        PtmPair Pair(IReadOnlyList<PtmPair> pairs) => pairs.Single(p => p.SiteA == s && p.SiteB == t);

        var givenDepth = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), runCovariate: depth));
        var expected = SpearmanCorrelation.PartialCorrelate(vs, vt, psms);
        Assert.That(givenDepth.Statistic, Is.EqualTo(expected.Rho).Within(1e-12));
        Assert.That(givenDepth.PValue, Is.EqualTo(expected.PValue).Within(1e-12));
        Assert.That(givenDepth.RunCovariateRemoved, Is.True);
        Assert.That(givenDepth.RunLevelRemoved, Is.False);

        // Ranked, so log2 of the PSMs gives the same answer.
        var log2 = depth.ToDictionary(kv => kv.Key, kv => Math.Log2(kv.Value));
        Assert.That(Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), runCovariate: log2)).Statistic,
            Is.EqualTo(givenDepth.Statistic).Within(1e-15));

        // With the run level too: both covariates.
        var level = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.035).ToArray();
        var both = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true, runCovariate: depth));
        Assert.That(both.Statistic, Is.EqualTo(SpearmanCorrelation.PartialCorrelate(vs, vt, new IReadOnlyList<double>[] { level, psms }).Rho).Within(1e-12));

        // A run without a covariate value is left out; below MinRunsForRunLevel shared runs the pair is not tested.
        var missingOne = depth.Where(kv => kv.Key != "r00").ToDictionary(kv => kv.Key, kv => kv.Value);
        Assert.That(Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), runCovariate: missingOne)).N, Is.EqualTo(11));
        var eight = depth.Where(kv => string.CompareOrdinal(kv.Key, "r08") < 0).ToDictionary(kv => kv.Key, kv => kv.Value);
        var few = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), runCovariate: eight));
        Assert.That(few.TestStatus, Is.EqualTo(PairTestStatus.TooFewSharedRuns));
        Assert.That(few.PValue, Is.NaN);
        Assert.Throws<ArgumentException>(() => PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(),
            runCovariate: new Dictionary<string, double> { ["r00"] = double.NaN }));
    }

    [Test]
    public void BoundsCountInTheRunLevelAtTheirRanks()
    {
        // Twelve runs. In the first six, three of the four background sites are covered-zeros: those runs read low.
        var s = new ModificationSite("P1", 13, 'S', Phos);
        var t = new ModificationSite("P2", 43, 'S', Phos);
        var background = Enumerable.Range(0, 4).Select(i => new ModificationSite($"P{3 + i}", 7, 'S', Phos)).ToList();
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        var occupancy = new List<SiteRunOccupancy>();
        for (int r = 0; r < 12; r++)
        {
            double batch = r < 6 ? 0.1 : 0.3;
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseS[r % 6], site: s));
            occupancy.Add(Cell($"r{r:D2}", batch + 0.01 * noiseT[r % 6], site: t));
            for (int i = 0; i < background.Count; i++)
                occupancy.Add(r < 6 && i > 0
                    ? Cell($"r{r:D2}", double.NaN, OccupancyState.NotDetected, site: background[i])
                    : Cell($"r{r:D2}", 0.3 + 0.01 * i + 0.001 * r, site: background[i]));
        }
        PtmPair Pair(IReadOnlyList<PtmPair> pairs) => pairs.Single(p => p.SiteA == s && p.SiteB == t);

        // The level's median over the four background sites: a bound in the first six runs (three of four are
        // covered-zeros), the middle two values in the rest.
        var withBounds = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true, includeBounds: true));
        var level = Enumerable.Range(0, 12).Select(r => r < 6 ? -5 : 0.315 + 0.001 * r).ToArray();
        var vs = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseS[r % 6]).ToArray();
        var vt = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6]).ToArray();
        Assert.That(withBounds.N, Is.EqualTo(12));
        Assert.That(withBounds.Statistic, Is.EqualTo(SpearmanCorrelation.PartialCorrelate(vs, vt, level).Rho).Within(1e-12));

        // Without them, three background sites are seen in half the runs, do not enter, and leave no run a level.
        var without = Pair(PtmPairEngine.CoVarying(occupancy, Array.Empty<PeptidoformObservation>(), removeRunLevel: true));
        Assert.That(without.TestStatus, Is.EqualTo(PairTestStatus.NotEstimable));
        Assert.That(without.N, Is.EqualTo(0));
    }

    private static PtmPair APair(string modA, int posA, string protA, int posB, string protB, double rho, double p, int n) => new()
    {
        ResultType = PairResultType.A, Overlapping = false, Statistic = rho, PValue = p, N = n,
        SiteA = new ModificationSite(protA, posA, 'S', modA), SiteB = new ModificationSite(protB, posB, 'S', Phos),
    };

    private static string Canonical(ModificationSite s) =>
        $"{s.ProteinAccession}:{s.Residue}{s.Position}:" + (s.Modification.EndsWith("Phosphoserine on S") || s.Modification.EndsWith("Phosphorylation on S") ? "UNIMOD:21" : s.Modification);

    [Test]
    public void GlobalTypeAIsSignedWeightedStoufferAcrossDatasetsWithNamesMerged()
    {
        var pairs = new[]
        {
            ("D1", APair(Phos, 13, "P1", 43, "P2", 0.8, 0.01, 16)),
            ("D2", APair("UniProt:Phosphoserine on S", 13, "P1", 43, "P2", 0.6, 0.04, 9)),   // same chemistry, other engine name
        };
        var g = GlobalPairEngine.Combine(pairs, Canonical).Single();
        Assert.That(g.Scopes, Is.EqualTo(new[] { "D1", "D2" }));
        Assert.That(g.ScopesAgreeing, Is.EqualTo(2));
        double z1 = MathNet.Numerics.Distributions.Normal.InvCDF(0, 1, 1 - 0.005), z2 = MathNet.Numerics.Distributions.Normal.InvCDF(0, 1, 1 - 0.02);
        double z = (4 * z1 + 3 * z2) / 5;
        Assert.That(g.CombinedZ, Is.EqualTo(z).Within(1e-10));
        Assert.That(g.PValue, Is.EqualTo(2 * MathNet.Numerics.Distributions.Normal.CDF(0, 1, -z)).Within(1e-12));
        Assert.That(g.Statistic, Is.EqualTo(0.7).Within(1e-12), "median rho");
        Assert.That(g.Q, Is.EqualTo(g.PValue).Within(1e-15), "one test in its family");
    }

    [Test]
    public void TwoNamesForOneSiteWithEqualRunsKeepTheFirstEngineKeyNotTheSmallerP()
    {
        // D1 has the pair under two names with equal N; "Common Biological:…" sorts before "UniProt:…" and is kept
        // although the UniProt row has the smaller p, in either input order.
        var common = ("D1", APair(Phos, 13, "P1", 43, "P2", 0.5, 0.1, 16));
        var uniProt = ("D1", APair("UniProt:Phosphoserine on S", 13, "P1", 43, "P2", 0.8, 0.01, 16));
        var d2 = ("D2", APair(Phos, 13, "P1", 43, "P2", 0.6, 0.04, 9));
        double z1 = MathNet.Numerics.Distributions.Normal.InvCDF(0, 1, 1 - 0.05), z2 = MathNet.Numerics.Distributions.Normal.InvCDF(0, 1, 1 - 0.02);
        foreach (var input in new[] { new[] { common, uniProt, d2 }, new[] { uniProt, common, d2 } })
        {
            var g = GlobalPairEngine.Combine(input, Canonical).Single();
            Assert.That(g.CombinedZ, Is.EqualTo((4 * z1 + 3 * z2) / 5).Within(1e-10));
            Assert.That(g.Statistic, Is.EqualTo(0.55).Within(1e-12), "median of 0.5 and 0.6");
        }
    }

    [Test]
    public void OppositeCorrelationsCancelAndSingleDatasetPairsAreLeftOut()
    {
        var pairs = new[]
        {
            ("D1", APair(Phos, 13, "P1", 43, "P2", 0.8, 0.01, 16)),
            ("D2", APair(Phos, 13, "P1", 43, "P2", -0.8, 0.01, 16)),
            ("D1", APair(Phos, 20, "P1", 43, "P2", 0.9, 0.001, 16)),                       // seen in one dataset only
        };
        var g = GlobalPairEngine.Combine(pairs, Canonical);
        Assert.That(g, Has.Count.EqualTo(1));
        Assert.That(g[0].CombinedZ, Is.EqualTo(0).Within(1e-12));
        Assert.That(g[0].PValue, Is.EqualTo(1).Within(1e-12));
    }

    [Test]
    public void GlobalTypePIsCountedAndUntestableARowsAreLeftOut()
    {
        PtmPair PPair(double coOccupancy, int n) => APair(Phos, 13, "P1", 20, "P1", coOccupancy, double.NaN, n) with
            { ResultType = PairResultType.P, Overlapping = true };
        var pairs = new[]
        {
            ("D1", PPair(0.25, 3)), ("D2", PPair(0.5, 2)), ("D3", PPair(double.NaN, 1)),   // D3: never quantified
            ("D1", APair(Phos, 13, "P1", 43, "P2", 0.8, 0.01, 16) with { Overlapping = true }),
            ("D2", APair(Phos, 13, "P1", 43, "P2", 0.8, 0.01, 16) with { Overlapping = true }),
            ("D1", APair(Phos, 30, "P1", 43, "P2", double.NaN, double.NaN, 2)),
            ("D2", APair(Phos, 30, "P1", 43, "P2", double.NaN, double.NaN, 2)),
            ("D1", APair("UniProt:Phosphoserine on S", 43, "P2", 43, "P2", 0.9, 0.01, 16)),   // one site under two names
            ("D2", APair("UniProt:Phosphoserine on S", 43, "P2", 43, "P2", 0.9, 0.01, 16)),
        };
        var g = GlobalPairEngine.Combine(pairs, Canonical).Single();
        Assert.That(g.ResultType, Is.EqualTo(PairResultType.P));
        Assert.That(g.Scopes, Is.EqualTo(new[] { "D1", "D2", "D3" }));
        Assert.That(g.ScopesAgreeing, Is.EqualTo(3));
        Assert.That(g.Statistic, Is.EqualTo(0.375).Within(1e-15), "median of 0.25 and 0.5; D3's NaN is not a value");
        Assert.That(g.CombinedZ, Is.NaN);
        Assert.That(g.PValue, Is.NaN);
        Assert.That(g.Q, Is.NaN);
        Assert.That(g.FdrFamily, Is.EqualTo("P:intra"));
        Assert.That(g.Method, Is.EqualTo("recurrence count"));
        Assert.That(GlobalPairEngine.Combine(pairs, Canonical, minScopes: 4), Is.Empty);
    }

    [Test]
    public void GlobalPoolingRefusesBadArguments()
    {
        var pairs = new[] { ("D1", APair(Phos, 13, "P1", 43, "P2", 0.8, 1.5, 16)), ("D2", APair(Phos, 13, "P1", 43, "P2", 0.8, 0.01, 16)) };
        Assert.Throws<ArgumentNullException>(() => GlobalPairEngine.Combine(null!, Canonical));
        Assert.Throws<ArgumentNullException>(() => GlobalPairEngine.Combine(pairs, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => GlobalPairEngine.Combine(pairs, Canonical, minScopes: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => GlobalPairEngine.Combine(pairs, Canonical), "p = 1.5 is not a p-value");
    }

    [Test]
    public void GlobalPoolingSurvivesZeroAndPerfectCorrelations()
    {
        (string, PtmPair) Scope(string scope, double rho, double p) => (scope, APair(Phos, 13, "P1", 43, "P2", rho, p, 10));
        GlobalPtmPair g = null!;

        // ρ = 0 in every scope: no direction to combine, so no test, and no exception.
        Assert.DoesNotThrow(() => g = GlobalPairEngine.Combine(new[] { Scope("D1", 0, 1), Scope("D2", 0, 1) }, Canonical).Single());
        Assert.That(g.CombinedZ, Is.NaN);
        Assert.That(g.PValue, Is.NaN);
        Assert.That(g.Q, Is.NaN);
        Assert.That(g.ScopesAgreeing, Is.EqualTo(0));

        // Spearman's asymptotic p is 0 at |ρ| = 1; opposite perfect scopes, or opposite tiny p, cancel instead of NaN.
        Assert.DoesNotThrow(() => g = GlobalPairEngine.Combine(new[] { Scope("D1", 1, 0), Scope("D2", -1, 0) }, Canonical).Single());
        Assert.That(g.CombinedZ, Is.EqualTo(0).Within(1e-12));
        Assert.That(g.PValue, Is.EqualTo(1).Within(1e-12));
        g = GlobalPairEngine.Combine(new[] { Scope("D1", 0.9, 1e-20), Scope("D2", -0.9, 1e-20) }, Canonical).Single();
        Assert.That(g.CombinedZ, Is.EqualTo(0).Within(1e-12));

        // One scope at p = 0 against a discordant one: finite Z, not an infinity that erases the other scope.
        g = GlobalPairEngine.Combine(new[] { Scope("D1", 1, 0), Scope("D2", -0.5, 0.2) }, Canonical).Single();
        Assert.That(double.IsFinite(g.CombinedZ), Is.True);
        Assert.That(g.ScopesAgreeing, Is.EqualTo(1));
    }

    private static readonly ModificationSite DeamN = new("P1", 10, 'N', "Common Artifact:Deamidation on N");

    private static SiteRunOccupancy Cell(string run, double fraction, OccupancyState state = OccupancyState.Quantified,
        bool unmodifiedQuantified = true, ModificationSite? site = null, double covering = 100) =>
        new(site ?? DeamN, run, state, fraction * covering, covering, unmodifiedQuantified) { ReportedFraction = fraction };

    private static double Expit(double x) => 1 / (1 + Math.Exp(-x));

    /// <summary>Six replicates at three trait levels, two runs each: logit = -2 + 0.5·trait + replicate offset ± 0.01.</summary>
    private static (List<SiteRunOccupancy> occ, Dictionary<string, double> trait, Dictionary<string, string> rep) TraitDesign()
    {
        var occ = new List<SiteRunOccupancy>();
        var trait = new Dictionary<string, double>();
        var rep = new Dictionary<string, string>();
        double[] levels = [0.25, 0.25, 1, 1, 2, 2];
        double[] offset = [0.1, -0.1, 0.05, -0.05, 0.08, -0.08];
        for (int r = 0; r < 6; r++)
            for (int inj = 0; inj < 2; inj++)
            {
                string run = $"R{r}_{inj}";
                trait[run] = levels[r];
                rep[run] = $"R{r}";
                // The last term keeps a residual once an injection covariate is fitted (else the fit is exact).
                occ.Add(Cell(run, Expit(-2 + 0.5 * levels[r] + offset[r] + (inj == 0 ? 0.01 : -0.01) + (inj == 1 ? 0.003 * (r % 3 - 1) : 0))));
            }
        return (occ, trait, rep);
    }

    [Test]
    public void TraitEffectIsRecoveredOnTheLogitScale()
    {
        var (occ, trait, rep) = TraitDesign();
        var e = SiteTraitEffects.Fit(occ, trait, rep, options: new SiteTraitOptions { LogitEpsilon = 1e-12 }).Single();
        Assert.That(e.Status, Is.EqualTo("Fitted"));
        Assert.That(e.Runs, Is.EqualTo(12));
        Assert.That(e.Replicates, Is.EqualTo(6));
        Assert.That(e.Effect, Is.EqualTo(0.5).Within(0.1));
        Assert.That(e.DegreesOfFreedom, Is.EqualTo(4));   // trait is constant within a replicate: G - 1 - 1
        Assert.That(e.Q, Is.EqualTo(e.PValue));
    }

    [Test]
    public void TraitFitUsesOnlyQuantifiedNonCeilingCells()
    {
        var (occ, trait, rep) = TraitDesign();
        // Replace one replicate's runs with a floor and a ceiling: both leave the fit, so 10 runs and 5 replicates remain.
        occ.RemoveAll(o => o.Run.StartsWith("R5_"));
        occ.Add(Cell("R5_0", double.NaN, OccupancyState.Floor));
        occ.Add(Cell("R5_1", 1, unmodifiedQuantified: false));
        var e = SiteTraitEffects.Fit(occ, trait, rep).Single();
        Assert.That(e.Runs, Is.EqualTo(10));
        Assert.That(e.Replicates, Is.EqualTo(5));

        var strict = SiteTraitEffects.Fit(occ, trait, rep, options: new SiteTraitOptions { MinReplicates = 6 }).Single();
        Assert.That(strict.Status, Is.EqualTo("BelowSupport"));
        Assert.That(strict.PValue, Is.NaN);
    }

    [Test]
    public void TraitFitReportsCovariateEffects()
    {
        var (occ, trait, rep) = TraitDesign();
        var cov = trait.Keys.ToDictionary(r => r, r => new[] { r.EndsWith("_0") ? 1.0 : 0.0 });
        var e = SiteTraitEffects.Fit(occ, trait, rep, cov, ["handling"], new SiteTraitOptions { LogitEpsilon = 1e-12 }).Single();
        Assert.That(e.Status, Is.EqualTo("Fitted"));
        Assert.That(e.CovariateEffects, Has.Count.EqualTo(1));
        Assert.That(e.CovariateEffects[0], Is.EqualTo(0.02).Within(0.005));   // the injection offset, +0.01 vs -0.01
        Assert.Throws<ArgumentException>(() => SiteTraitEffects.Fit(occ, trait, rep, null, ["handling"]));
    }

    [Test]
    public void TraitFitOnLog2OddsIsTheNaturalLogFitDividedByLn2()
    {
        var (occ, trait, rep) = TraitDesign();
        var cov = trait.Keys.ToDictionary(r => r, r => new[] { r.EndsWith("_0") ? 1.0 : 0.0 });
        var ln = SiteTraitEffects.Fit(occ, trait, rep, cov, ["handling"]).Single();
        var log2 = SiteTraitEffects.Fit(occ, trait, rep, cov, ["handling"],
            new SiteTraitOptions { Scale = OccupancyScale.Log2Odds }).Single();
        double ln2 = Math.Log(2);
        Assert.That(ln.Scale, Is.EqualTo(OccupancyScale.NaturalLogOdds));
        Assert.That(log2.Scale, Is.EqualTo(OccupancyScale.Log2Odds));
        Assert.That(log2.Status, Is.EqualTo("Fitted"));
        // Exact: the log2 result is the natural-log fit divided by ln 2, not a second fit.
        Assert.That(log2.Effect, Is.EqualTo(ln.Effect / ln2));
        Assert.That(log2.StandardError, Is.EqualTo(ln.StandardError / ln2));
        Assert.That(log2.CovariateEffects[0], Is.EqualTo(ln.CovariateEffects[0] / ln2));
        Assert.That(log2.ReplicateVariance, Is.EqualTo(ln.ReplicateVariance / (ln2 * ln2)));
        Assert.That(log2.ResidualVariance, Is.EqualTo(ln.ResidualVariance / (ln2 * ln2)));
        Assert.That(log2.DegreesOfFreedom, Is.EqualTo(ln.DegreesOfFreedom));
        Assert.That(log2.PValue, Is.EqualTo(ln.PValue));
        Assert.That(log2.Q, Is.EqualTo(ln.Q));
        Assert.That(log2.MedianFraction, Is.EqualTo(ln.MedianFraction));

        // A site below support records its scale too, and an unknown scale is refused.
        var below = SiteTraitEffects.Fit(occ, trait, rep,
            options: new SiteTraitOptions { Scale = OccupancyScale.Log2Odds, MinReplicates = 7 }).Single();
        Assert.That(below.Status, Is.EqualTo("BelowSupport"));
        Assert.That(below.Scale, Is.EqualTo(OccupancyScale.Log2Odds));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SiteTraitEffects.Fit(occ, trait, rep, options: new SiteTraitOptions { Scale = (OccupancyScale)2 }));
    }

    private static SiteRunOccupancy WithCounts(SiteRunOccupancy cell, int modified, int covering) =>
        cell with { CoveringPeptidoforms = covering, ModifiedPeptidoforms = modified };

    /// <summary>TraitDesign with counts on every cell, and the first three runs of the low-trait replicates covered-zeros (0 of 2).</summary>
    private static (List<SiteRunOccupancy> occ, Dictionary<string, double> trait, Dictionary<string, string> rep) SampledTraitDesign()
    {
        var (occ, trait, rep) = TraitDesign();
        var zeros = new HashSet<string> { "R0_0", "R0_1", "R1_0" };
        occ = occ.Select(o => zeros.Contains(o.Run)
            ? WithCounts(Cell(o.Run, double.NaN, OccupancyState.NotDetected), 0, 2)
            : WithCounts(o, 1, 2)).ToList();
        return (occ, trait, rep);
    }

    [Test]
    public void SampledRunsEnterTheTraitFitThroughTheLink()
    {
        var (occ, trait, rep) = SampledTraitDesign();
        var link = new DetectionLink(0, 1);
        var options = new SiteTraitOptions { BoundCells = BoundCellModel.Sampling, DetectionLink = link };
        var e = SiteTraitEffects.Fit(occ, trait, rep, options: options).Single();
        Assert.That(e.Status, Is.EqualTo("Fitted"));
        Assert.That(e.Runs, Is.EqualTo(12));
        Assert.That(e.CountedRuns, Is.EqualTo(3));
        Assert.That(e.Link, Is.EqualTo(link));
        Assert.That(e.ResidualVarianceBorrowed, Is.False);

        // The same numbers as the model fitted directly on the measured log odds and the three counts.
        var runs = trait.Keys.OrderBy(r => r, StringComparer.Ordinal).ToList();
        var y = new double[1, runs.Count];
        var counts = new BinomialCount?[1, runs.Count];
        var design = new double[runs.Count, 2];
        for (int i = 0; i < runs.Count; i++)
        {
            design[i, 0] = 1;
            design[i, 1] = trait[runs[i]];
            var cell = occ.Single(o => o.Run == runs[i]);
            double f = Math.Clamp(cell.Fraction, 5e-5, 1 - 5e-5);
            y[0, i] = cell.State == OccupancyState.Quantified ? Math.Log(f / (1 - f)) : double.NaN;
            if (cell.State == OccupancyState.NotDetected) counts[0, i] = new BinomialCount(0, 2);
        }
        var direct = LogitNormalMixedModel.Fit(y, counts, design, runs.Select(r => rep[r]).ToArray(), link);
        Assert.That(e.Effect, Is.EqualTo(direct.Coefficient(0, 1)).Within(1e-12));
        Assert.That(e.StandardError, Is.EqualTo(direct.StandardError(0, 1)).Within(1e-12));
        Assert.That(e.PValue, Is.EqualTo(direct.PValue(0, 1)).Within(1e-12));

        // By default the covered-zeros are left out, as before.
        var excluded = SiteTraitEffects.Fit(occ, trait, rep).Single();
        Assert.That(excluded.Runs, Is.EqualTo(9));
        Assert.That(excluded.CountedRuns, Is.EqualTo(0));
        Assert.That(excluded.Link, Is.Null);

        Assert.Throws<ArgumentException>(() => SiteTraitEffects.Fit(occ, trait, rep,
            options: options with { ExcludeCeiling = true }));
        Assert.Throws<ArgumentException>(() => SiteTraitEffects.Fit(occ, trait, rep,
            options: options with { DetectionLink = link with { DepthSlope = 0.2 } }));
        Assert.Throws<ArgumentException>(() => SiteTraitEffects.Fit(occ, trait, rep,
            options: options with { RunDepth = new Dictionary<string, double> { ["R0_0"] = 10 } }));
    }

    [Test]
    public void ASiteMostlyAtOneBoundIsNotFittedAndOneRarelyMeasuredBorrowsTheScatter()
    {
        var (occ, trait, rep) = TraitDesign();
        occ = occ.Select(o => WithCounts(o, 1, 2)).ToList();                     // DeamN: measured everywhere
        var runs = trait.Keys.OrderBy(r => r, StringComparer.Ordinal).ToList();
        var rare = new ModificationSite("P2", 20, 'N', DeamN.Modification);
        var bound = new ModificationSite("P3", 30, 'N', DeamN.Modification);
        for (int i = 0; i < runs.Count; i++)
        {
            // rare: two measured runs, the rest 1 of 3 (off the extremes); bound: one measured run, the rest 0 of 3.
            occ.Add(i < 2 ? WithCounts(Cell(runs[i], 0.3 + 0.01 * i, site: rare), 1, 3)
                : WithCounts(Cell(runs[i], double.NaN, OccupancyState.Floor, site: rare), 1, 3));
            occ.Add(i < 1 ? WithCounts(Cell(runs[i], 0.2, site: bound), 1, 3)
                : WithCounts(Cell(runs[i], double.NaN, OccupancyState.NotDetected, site: bound), 0, 3));
        }
        var options = new SiteTraitOptions { BoundCells = BoundCellModel.Sampling, DetectionLink = new DetectionLink(0, 1) };
        var results = SiteTraitEffects.Fit(occ, trait, rep, options: options);
        var measured = results.Single(r => r.Site == DeamN);
        var borrowed = results.Single(r => r.Site == rare);
        Assert.That(measured.Status, Is.EqualTo("Fitted"));
        Assert.That(measured.CountedRuns, Is.EqualTo(0), "no counted run: fitted as before");
        Assert.That(borrowed.ResidualVarianceBorrowed, Is.True);
        Assert.That(borrowed.ResidualVariance, Is.EqualTo(measured.ResidualVariance).Within(1e-12 * measured.ResidualVariance),
            "the median over the sites with their own: here the one");
        Assert.That(results.Single(r => r.Site == bound).Status, Is.EqualTo("MostlyOneBound"));
        Assert.That(results.Single(r => r.Site == bound).PValue, Is.NaN);
    }

    [Test]
    public void TheLinkIsFittedPerCallOrFallsBack()
    {
        var (occ, trait, rep) = SampledTraitDesign();
        var sampling = new SiteTraitOptions { BoundCells = BoundCellModel.Sampling };
        // One site's 12 cells are too few to fit a link from: without a fallback the site is not fitted.
        Assert.That(SiteTraitEffects.Fit(occ, trait, rep, options: sampling).Single().Status, Is.EqualTo("NoDetectionLink"));
        var fallback = new DetectionLink(0.1, 0.8);
        var withFallback = SiteTraitEffects.Fit(occ, trait, rep, options: sampling with { FallbackLink = fallback }).Single();
        Assert.That(withFallback.Status, Is.EqualTo("Fitted"));
        Assert.That(withFallback.Link, Is.EqualTo(fallback));

        // Nine more sites like it give 120 cells: the call fits its own link.
        var more = new List<SiteRunOccupancy>(occ);
        for (int s = 1; s <= 9; s++)
        {
            var site = new ModificationSite($"Q{s}", 5 + s, 'N', DeamN.Modification);
            more.AddRange(occ.Select((o, i) => o with
            {
                Site = site,
                ReportedFraction = o.State == OccupancyState.Quantified ? Math.Min(0.95, o.Fraction * (1 + 0.1 * s)) : o.ReportedFraction,
                ModifiedPeptidoforms = o.State == OccupancyState.Quantified ? 1 + (i + s) % 2 : 0,
                CoveringPeptidoforms = 2,
            }));
        }
        var own = SiteTraitEffects.Fit(more, trait, rep, options: sampling with { FallbackLink = fallback });
        Assert.That(own.Select(e => e.Link).Distinct().Count(), Is.EqualTo(1), "one link per call");
        Assert.That(own[0].Link, Is.Not.Null.And.Not.EqualTo(fallback));
        Assert.That(own[0].Link!.Slope, Is.GreaterThan(0), "more occupancy, more often seen");
    }

    [Test]
    public void CombiningFractionsSumsThePeptidoformCountsAndRefusesUnknownStates()
    {
        var map = new Dictionary<string, string> { ["f1"] = "s", ["f2"] = "s" };
        var cells = new[]
        {
            WithCounts(Cell("f1", 0.2), 1, 3),
            WithCounts(Cell("f2", double.NaN, OccupancyState.NotDetected), 0, 2),
        };
        var combined = SiteOccupancyCalculator.CombineRuns(cells, map).Single();
        Assert.That(combined.CoveringPeptidoforms, Is.EqualTo(5));
        Assert.That(combined.ModifiedPeptidoforms, Is.EqualTo(1));
        var missing = SiteOccupancyCalculator.CombineRuns(new[] { cells[0], Cell("f2", double.NaN, OccupancyState.NotDetected) }, map).Single();
        Assert.That(missing.CoveringPeptidoforms, Is.Null, "a run without counts leaves the sum unknown");
        Assert.Throws<ArgumentException>(() => SiteOccupancyCalculator.CombineRuns(
            new[] { Cell("f1", 0.2, (OccupancyState)99) }, map));
    }

    [Test]
    public void PooledOccupancyIsModifiedShareOfCoveringSignal()
    {
        var q1 = new ModificationSite("P1", 3, 'Q', "Common Artifact:Deamidation on Q");
        var q2 = new ModificationSite("P2", 7, 'Q', "Common Artifact:Deamidation on Q");
        var pooled = SiteTraitEffects.PooledOccupancy(new[]
        {
            Cell("r1", 0.1, site: q1, covering: 100),
            Cell("r1", 0.4, site: q2, covering: 300),
            Cell("r1", double.NaN, OccupancyState.Floor, site: q2),
            Cell("r2", 0.2, site: q1, covering: 50),
        });
        Assert.That(pooled["r1"], Is.EqualTo((10 + 120) / 400.0).Within(1e-15));
        Assert.That(pooled["r2"], Is.EqualTo(0.2).Within(1e-15));
    }
}
