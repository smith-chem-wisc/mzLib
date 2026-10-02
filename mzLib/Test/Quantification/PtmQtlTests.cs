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
    public void PeptideNTerminalModificationsAndNonMetStartsAreNotProteinNTermini()
    {
        var occ = SiteOccupancyCalculator.Calculate(new[]
        {
            Obs("r1", "[Common Artifact:Ammonia loss on C]CPEPK", 1, 10, "P1"),                   // peptide N-terminal mod at residue 1
            new PeptidoformObservation("r1", "[UniProt:N-acetylserine on S]SPEPK", "P2", 2, 6, 10, 'K'), // residue 2 after K1: no Met removal
            new PeptidoformObservation("r1", "[UniProt:N-acetylserine on S]SPEPK", "P3", 2, 6, 10, 'M'), // residue 2 after Met removal
            Obs("r1", "[UniProt:N-acetylserine on S]SPEPK", 2, 10, "P4"),                            // unknown: the caller guarantees Met
        });
        Assert.That(occ.Select(o => o.Site.Key).Distinct(),
            Is.EquivalentTo(new[] { "P3:S2:UniProt:N-acetylserine on S", "P4:S2:UniProt:N-acetylserine on S" }));
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
