using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using StatisticalModels;

namespace Quantification.PtmQtl;

/// <summary>Kind of evidence a pair row carries.</summary>
public enum PairResultType
{
    /// <summary>Physical: both modifications were identified on ONE peptidoform, i.e. on one molecule.</summary>
    P,
    /// <summary>Co-varying: the two sites' occupancies rise and fall together (or oppositely) across runs.</summary>
    A,
}

/// <summary>Whether a type A pair was tested, and if not, why.</summary>
public enum PairTestStatus
{
    /// <summary>Tested: <see cref="PtmPair.PValue"/> is set.</summary>
    Tested,
    /// <summary>
    /// The correlation could not be computed: fewer complete runs than it needs, a site constant over them, or a
    /// covariate that carries a site entirely.
    /// </summary>
    NotEstimable,
    /// <summary>A partial correlation on fewer than <see cref="PtmPairEngine.MinRunsForRunLevel"/> shared runs.</summary>
    TooFewSharedRuns,
    /// <summary>
    /// With bounds included, both sites are a bound in most shared runs: fewer than
    /// <see cref="PtmPairEngine.MinUntiedShare"/> of the runs are outside each site's largest tie.
    /// </summary>
    BothSitesMostlyBounds,
}

/// <summary>One relationship between two modification sites, within one scope (e.g. one dataset).</summary>
/// <remarks>
/// A pair is unordered and is written once, with <see cref="SiteA"/>'s key before <see cref="SiteB"/>'s by
/// ordinal comparison.
/// </remarks>
public sealed record PtmPair
{
    /// <summary>Evidence type.</summary>
    public required PairResultType ResultType { get; init; }
    /// <summary>First site (ordinal key order).</summary>
    public required ModificationSite SiteA { get; init; }
    /// <summary>Second site.</summary>
    public required ModificationSite SiteB { get; init; }
    /// <summary>Both sites on the same protein accession.</summary>
    public bool SameProtein => SiteA.ProteinAccession == SiteB.ProteinAccession;
    /// <summary>
    /// Some peptide covers both sites: one span covering both positions of one protein, or one peptide mapped to
    /// both proteins (a shared peptide, e.g. between isoforms). For type A the two occupancies then share
    /// peptides and are correlated or anti-correlated by construction, so overlapping pairs are reported but
    /// excluded from the A family's multiple-testing adjustment (q is NaN).
    /// </summary>
    public required bool Overlapping { get; init; }
    /// <summary>
    /// Type P: median over runs of the co-occupancy, the intensity of peptidoforms carrying both modifications
    /// over the intensity of peptidoforms covering both positions (NaN when never quantified).
    /// Type A: Spearman's ρ of the two occupancies across runs; the partial ρ given each run's level when
    /// <see cref="RunLevelRemoved"/>, and given the run covariate when <see cref="RunCovariateRemoved"/>. NaN when not
    /// tested (<see cref="TestStatus"/>).
    /// </summary>
    public required double Statistic { get; init; }
    /// <summary>Type A: two-sided Spearman p (method in <see cref="SpearmanMethod"/>). Type P: NaN (no test in this version).</summary>
    public required double PValue { get; init; }
    /// <summary>Benjamini-Hochberg adjusted p within <see cref="FdrFamily"/>. NaN when not tested or excluded.</summary>
    public double Q { get; internal set; } = double.NaN;
    /// <summary>The family the adjustment ran over: result type, then same- or different-protein.</summary>
    public string FdrFamily => $"{ResultType}:{(SameProtein ? "intra" : "inter")}";
    /// <summary>
    /// Type P: runs where a doubly modified peptidoform was identified. Type A: runs where both sites have a value
    /// (with <see cref="BoundsIncluded"/>, a value or a bound) and every covariate is known.
    /// </summary>
    public required int N { get; init; }
    /// <summary>
    /// Runs <see cref="Statistic"/> was computed on. Type P: runs where the doubly modified form and the forms
    /// covering both positions were quantified, at most <see cref="N"/>. Type A: equal to <see cref="N"/>.
    /// </summary>
    public int StatisticN { get; init; }
    /// <summary>Type A: how the Spearman p-value was computed. Type P: null.</summary>
    public SpearmanPValueMethod? SpearmanMethod { get; init; }
    /// <summary>
    /// Type A: whether <see cref="Statistic"/> and <see cref="PValue"/> are the partial Spearman correlation given
    /// each run's level (the median occupancy over the scope's sites in that run). Type P: false.
    /// </summary>
    public bool RunLevelRemoved { get; init; }
    /// <summary>Type A: whether the correlation is partial given the caller's per-run covariate (e.g. depth). Type P: false.</summary>
    public bool RunCovariateRemoved { get; init; }
    /// <summary>
    /// Type A: whether runs where a site was seen in only one form (or as a floor) entered as bounds, tied lowest or
    /// highest for that site. Type P: false.
    /// </summary>
    public bool BoundsIncluded { get; init; }
    /// <summary>Type A: whether the pair was tested, and if not, why. Type P: null.</summary>
    public PairTestStatus? TestStatus { get; init; }
}

/// <summary>
/// PTM-PTM relationships within one scope (one dataset): physical co-occurrence on a molecule (type P) and
/// co-variation of site occupancies across runs (type A).
/// </summary>
public static class PtmPairEngine
{
    /// <summary>A site enters type A only if its occupancy is quantified in at least this fraction of runs.</summary>
    public const double DefaultMinQuantifiedFraction = 0.7;
    /// <summary>
    /// With the run level removed, the fewest runs a pair must share to be tested. The partial test has no exact
    /// p-value, and its t approximation with n − 3 df calls a perfect ρ at n = 4 or 5 significant; plain Spearman
    /// is protected there by its exact permutation p. Below this the pair is NotEstimable and left out of the
    /// adjustment.
    /// </summary>
    public const int MinRunsForRunLevel = 10;
    /// <summary>
    /// With the run level removed, the fewest other sites of a group (the pair's own two left out) that must have a
    /// value in a run for that run to have a level. Below this the run is left out of the pair's test.
    /// </summary>
    public const int MinSitesForRunLevel = 3;
    /// <summary>
    /// With bounds included, a pair is tested only if at least one of its sites has at least this share of the
    /// shared runs outside its largest tie (its low bounds, or its high bounds, whichever are more). When both
    /// sites are mostly one bound, Spearman's t test calls too many pairs: with no relation at all, 2 to 5 times
    /// the nominal rate at p &lt; 0.001 in simulation, and two sites with values in the same 2 of 40 runs read as a
    /// perfect correlation 1 time in 780. With one site at this share or more it is right or conservative.
    /// </summary>
    public const double MinUntiedShare = 0.5;
    /// <summary>The value a low bound (a covered-zero or a floor) takes in a site's vector: below every fraction.</summary>
    private const double LowBound = -1;
    /// <summary>The value a high bound (seen only modified) takes in a site's vector: above every fraction.</summary>
    private const double HighBound = 2;

    /// <summary>
    /// Type P: every pair of sites carried together by at least one peptidoform, with the runs it was
    /// identified in and its median co-occupancy. The observations are validated as in
    /// <see cref="SiteOccupancyCalculator.Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/>.
    /// </summary>
    public static IReadOnlyList<PtmPair> Physical(IEnumerable<PeptidoformObservation> observations,
        Func<string, bool>? includeModification = null)
    {
        var obs = SiteOccupancyCalculator.Parse(observations, includeModification).Select(p => (o: p.obs, p.sites)).ToList();

        var runsTogether = new Dictionary<(ModificationSite, ModificationSite), HashSet<string>>();
        foreach (var (o, sites) in obs)
        {
            if (sites.Count < 2) continue;
            foreach (var (a, b) in Pairs(sites.Distinct().ToList()))
            {
                if (!runsTogether.TryGetValue((a, b), out var runs)) runsTogether[(a, b)] = runs = new HashSet<string>();
                runs.Add(o.Run);
            }
        }

        var byProteinRun = obs.GroupBy(x => (x.o.ProteinAccession, x.o.Run)).ToDictionary(g => g.Key, g => g.ToList());
        var result = new List<PtmPair>();
        foreach (var ((a, b), runs) in runsTogether.OrderBy(kv => kv.Key.Item1.Key, StringComparer.Ordinal)
                                                   .ThenBy(kv => kv.Key.Item2.Key, StringComparer.Ordinal))
        {
            var coOccupancy = new List<double>();
            foreach (var run in runs)
            {
                double both = 0, covering = 0;
                foreach (var (o, sites) in byProteinRun[(a.ProteinAccession, run)])
                {
                    if (double.IsNaN(o.Intensity)) continue;
                    int lo = Math.Min(a.Position, b.Position), hi = Math.Max(a.Position, b.Position);
                    if (lo < o.StartResidue || hi > o.EndResidue) continue;
                    covering += o.Intensity;
                    if (sites.Contains(a) && sites.Contains(b)) both += o.Intensity;
                }
                if (covering > 0 && both > 0) coOccupancy.Add(both / covering);
            }
            result.Add(new PtmPair
            {
                ResultType = PairResultType.P, SiteA = a, SiteB = b, Overlapping = true,
                Statistic = Median(coOccupancy), PValue = double.NaN, N = runs.Count, StatisticN = coOccupancy.Count,
            });
        }
        return result;
    }

    /// <summary>
    /// Type A: Spearman correlation across runs between the occupancies of every pair of sites quantified in
    /// at least <paramref name="minQuantifiedFraction"/> of the scope's runs, on the runs where both are
    /// quantified. By default, floors and undetected runs are not values and are left out, and so are ceilings
    /// (runs where a site was seen only modified): two unrelated sites both read 1 in low-load runs, which would
    /// correlate them through detection. With <paramref name="includeBounds"/> they enter as bounds instead, and the
    /// run's depth can be removed through <paramref name="runCovariate"/>. BH within each family, with overlapping
    /// pairs and untested pairs excluded from the family.
    /// </summary>
    /// <param name="occupancy">Output of <see cref="SiteOccupancyCalculator.Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/> for one scope.</param>
    /// <param name="observations">The same observations, used to decide which pairs overlap.</param>
    /// <param name="minQuantifiedFraction">
    /// The pair-scale guard of D4 (default 0.7): the share of the scope's runs where a site must have a value (with
    /// <paramref name="includeBounds"/>, a value or a bound) to enter.
    /// </param>
    /// <param name="excludeCeiling">
    /// Leave out a run where the site was seen only modified (occupancy 1 as a ceiling,
    /// <see cref="SiteRunOccupancy.UnmodifiedQuantified"/> false), as <see cref="SiteTraitOptions.ExcludeCeiling"/>
    /// does for traits. Null (default): true, or false with <paramref name="includeBounds"/>, where a ceiling is a
    /// bound; true together with <paramref name="includeBounds"/> is refused.
    /// </param>
    /// <param name="removeRunLevel">
    /// Correlate given each run's level: the partial Spearman correlation
    /// (<see cref="SpearmanCorrelation.PartialCorrelate(IReadOnlyList{double}, IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>)
    /// with the run's median occupancy over the other tested sites of the pair's group (see
    /// <paramref name="runLevelGroup"/>; the sites that pass <paramref name="minQuantifiedFraction"/>) that have a value in
    /// that run, the same cells the pairs use. The pair's own two sites are left out of the
    /// median, so a small group is not mostly the pair itself; a run whose group has fewer than
    /// <see cref="MinSitesForRunLevel"/> other sites there has no level. A shift that moves the group's sites in the same
    /// runs, such as an acquisition batch or sample load, then no longer correlates every pair with every other. A pair
    /// sharing fewer than <see cref="MinRunsForRunLevel"/> runs with a level is not tested. Default false.
    /// Batch only: where the run level follows the design (one condition reads higher), removing it tests the
    /// correlation within conditions instead, which is a different question.
    /// </param>
    /// <param name="runLevelGroup">
    /// With <paramref name="removeRunLevel"/>: the group whose run level applies to a site, for example its
    /// modification's chemistry. A pair within one group is correlated given that group's level; a pair across two
    /// groups, given both. Null (default): one group, every site.
    /// </param>
    /// <param name="includeBounds">
    /// Let a run where a site was seen in only one form enter as a bound, never dropped and never given an invented
    /// value: a low bound (<see cref="OccupancyState.NotDetected"/>, a covered-zero; or <see cref="OccupancyState.Floor"/>)
    /// is tied below every fraction of that site, and a high bound (a ceiling; or a
    /// <see cref="OccupancyState.CountOnly"/> cell whose covering peptidoforms all carry the modification) is tied
    /// above every fraction. A count-only cell where both forms were identified, or whose peptidoform counts were not
    /// supplied, has no rank and is left out. Bounds count toward <paramref name="minQuantifiedFraction"/> and toward a
    /// run's level. A pair whose two sites are both mostly one bound is not tested (<see cref="MinUntiedShare"/>,
    /// <see cref="PairTestStatus.BothSitesMostlyBounds"/>). Default false.
    /// </param>
    /// <param name="runCovariate">
    /// A value per run (e.g. its depth, log2 of its target PSMs) to correlate given: shallow runs make covered-zeros
    /// co-occur across unrelated sites. Ranked like every variable, so any increasing transform gives the same
    /// result. A run missing from it is left out of every pair. As with the run level, a pair sharing fewer than
    /// <see cref="MinRunsForRunLevel"/> runs is not tested. Null (default): none.
    /// </param>
    public static IReadOnlyList<PtmPair> CoVarying(IReadOnlyList<SiteRunOccupancy> occupancy,
        IEnumerable<PeptidoformObservation> observations, double minQuantifiedFraction = DefaultMinQuantifiedFraction,
        bool? excludeCeiling = null, bool removeRunLevel = false, Func<ModificationSite, string>? runLevelGroup = null,
        bool includeBounds = false, IReadOnlyDictionary<string, double>? runCovariate = null)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        ArgumentNullException.ThrowIfNull(observations);
        if (!(minQuantifiedFraction > 0 && minQuantifiedFraction <= 1)) throw new ArgumentOutOfRangeException(nameof(minQuantifiedFraction));
        if (includeBounds && excludeCeiling == true)
            throw new ArgumentException("With includeBounds a ceiling is a bound, tied highest; leave excludeCeiling unset or false.", nameof(excludeCeiling));
        if (runCovariate != null && runCovariate.Values.Any(v => !double.IsFinite(v)))
            throw new ArgumentException("Every run covariate value must be finite; leave a run out to leave it unknown.", nameof(runCovariate));
        bool dropCeiling = excludeCeiling ?? !includeBounds;

        var runs = occupancy.Select(o => o.Run).Distinct().OrderBy(r => r, StringComparer.Ordinal).ToArray();
        var runIndex = runs.Select((r, i) => (r, i)).ToDictionary(t => t.r, t => t.i, StringComparer.Ordinal);
        var vectors = new Dictionary<ModificationSite, double[]>();
        foreach (var o in occupancy)
        {
            double value = ValueOf(o, dropCeiling, includeBounds);
            if (double.IsNaN(value)) continue;
            if (!vectors.TryGetValue(o.Site, out var v))
                vectors[o.Site] = v = Enumerable.Repeat(double.NaN, runs.Length).ToArray();
            v[runIndex[o.Run]] = value;
        }
        int needed = (int)Math.Ceiling(minQuantifiedFraction * runs.Length);
        var sites = vectors.Where(kv => kv.Value.Count(double.IsFinite) >= needed).Select(kv => kv.Key)
            .OrderBy(s => s.Key, StringComparer.Ordinal).ToList();
        // Per group and run, the sorted values of the group's tested sites (the ones that enter pairs) with a value there;
        // with bounds included, a bound counts at its rank (lowest or highest), so a run where many sites fell to
        // covered-zero reads low.
        // Using only these keeps a run's level over the same sites in every run; sites seen in a few runs would make it
        // track which sites were detected.
        string GroupOf(ModificationSite s) => runLevelGroup?.Invoke(s) ?? "";
        var sortedByGroup = removeRunLevel
            ? sites.Select(st => (st, v: vectors[st])).GroupBy(t => GroupOf(t.st), StringComparer.Ordinal).ToDictionary(g => g.Key,
                g => Enumerable.Range(0, runs.Length).Select(i => g.Select(t => t.v[i]).Where(double.IsFinite).OrderBy(v => v).ToArray()).ToArray(),
                StringComparer.Ordinal)
            : null;
        double[] Level(string group, ModificationSite a, ModificationSite b)
        {
            var perRun = sortedByGroup![group];
            var level = new double[runs.Length];
            for (int i = 0; i < runs.Length; i++)
            {
                var drop = new List<double>(2);
                if (GroupOf(a) == group && double.IsFinite(vectors[a][i])) drop.Add(vectors[a][i]);
                if (GroupOf(b) == group && double.IsFinite(vectors[b][i])) drop.Add(vectors[b][i]);
                level[i] = perRun[i].Length - drop.Count >= MinSitesForRunLevel ? MedianExcluding(perRun[i], drop) : double.NaN;
            }
            return level;
        }

        // The peptides (base sequences) covering each site, on any protein they map to. Two sites overlap when one
        // peptide covers both: one span on one protein, or one shared peptide mapped to both proteins.
        var spans = observations
            .Select(o => (o.ProteinAccession, o.StartResidue, o.EndResidue, Peptide: o.FullSequence.GetBaseSequenceFromFullSequence()))
            .Distinct().GroupBy(o => o.ProteinAccession, StringComparer.Ordinal)
            .ToDictionary(g => g.Key, g => g.ToList(), StringComparer.Ordinal);
        var coveringPeptides = sites.ToDictionary(s => s, s => spans.TryGetValue(s.ProteinAccession, out var sp)
            ? sp.Where(o => o.StartResidue <= s.Position && s.Position <= o.EndResidue).Select(o => o.Peptide).ToHashSet(StringComparer.Ordinal)
            : new HashSet<string>(StringComparer.Ordinal));
        bool Overlap(ModificationSite a, ModificationSite b) => coveringPeptides[a].Overlaps(coveringPeptides[b]);

        double[]? covariate = runCovariate == null ? null
            : runs.Select(r => runCovariate.TryGetValue(r, out double v) ? v : double.NaN).ToArray();
        bool partial = removeRunLevel || covariate != null;

        var result = new List<PtmPair>();
        foreach (var (a, b) in Pairs(sites))
        {
            var covariates = new List<IReadOnlyList<double>>();
            if (removeRunLevel)
                covariates.AddRange(new[] { GroupOf(a), GroupOf(b) }.Distinct(StringComparer.Ordinal).Select(g => (IReadOnlyList<double>)Level(g, a, b)));
            if (covariate != null) covariates.Add(covariate);
            var r = partial
                ? SpearmanCorrelation.PartialCorrelate(vectors[a], vectors[b], covariates)
                : SpearmanCorrelation.Correlate(vectors[a], vectors[b]);

            var status = r.Method == SpearmanPValueMethod.NotEstimable ? PairTestStatus.NotEstimable : PairTestStatus.Tested;
            if (status == PairTestStatus.Tested && partial && r.N < MinRunsForRunLevel) status = PairTestStatus.TooFewSharedRuns;
            if (status == PairTestStatus.Tested && includeBounds && BothMostlyBounds(vectors[a], vectors[b], covariates))
                status = PairTestStatus.BothSitesMostlyBounds;
            bool tested = status == PairTestStatus.Tested;
            result.Add(new PtmPair
            {
                ResultType = PairResultType.A, SiteA = a, SiteB = b, Overlapping = Overlap(a, b),
                Statistic = tested ? r.Rho : double.NaN, PValue = tested ? r.PValue : double.NaN, N = r.N, StatisticN = r.N,
                SpearmanMethod = tested ? r.Method : SpearmanPValueMethod.NotEstimable,
                RunLevelRemoved = removeRunLevel, RunCovariateRemoved = covariate != null, BoundsIncluded = includeBounds,
                TestStatus = status,
            });
        }
        foreach (var family in result.Where(p => !p.Overlapping).GroupBy(p => p.FdrFamily))
        {
            var members = family.ToList();
            var q = MultipleTesting.BenjaminiHochberg(members.Select(p => p.PValue).ToArray());
            for (int i = 0; i < members.Count; i++) members[i].Q = q[i];
        }
        return result;
    }

    /// <summary>
    /// A cell's entry in its site's type A vector: its fraction, a bound (<see cref="LowBound"/>,
    /// <see cref="HighBound"/>), or NaN when it does not enter.
    /// </summary>
    private static double ValueOf(SiteRunOccupancy o, bool dropCeiling, bool includeBounds)
    {
        switch (o.State)
        {
            case OccupancyState.Quantified when o.UnmodifiedQuantified:
                return o.Fraction;
            case OccupancyState.Quantified:
                return includeBounds ? HighBound : dropCeiling ? double.NaN : o.Fraction;
            case OccupancyState.NotDetected:
            case OccupancyState.Floor:
                return includeBounds ? LowBound : double.NaN;
            case OccupancyState.CountOnly:
                // Seen only modified: a high bound. Both forms identified, or counts unknown: no rank.
                return includeBounds && o.CoveringPeptidoforms is int n && n > 0 && o.ModifiedPeptidoforms == n ? HighBound : double.NaN;
            default:
                throw new ArgumentException($"Unknown occupancy state {o.State} for {o.Site.Key} in {o.Run}.", "occupancy");
        }
    }

    /// <summary>
    /// True when, over the runs the test uses (both sites and every covariate known), each site has fewer than
    /// <see cref="MinUntiedShare"/> of them outside its largest tie of bounds.
    /// </summary>
    private static bool BothMostlyBounds(double[] a, double[] b, List<IReadOnlyList<double>> covariates)
    {
        int n = 0, lowA = 0, highA = 0, lowB = 0, highB = 0;
        for (int i = 0; i < a.Length; i++)
        {
            if (!double.IsFinite(a[i]) || !double.IsFinite(b[i]) || covariates.Any(c => !double.IsFinite(c[i]))) continue;
            n++;
            if (a[i] == LowBound) lowA++; else if (a[i] == HighBound) highA++;
            if (b[i] == LowBound) lowB++; else if (b[i] == HighBound) highB++;
        }
        bool Mostly(int low, int high) => n - Math.Max(low, high) < MinUntiedShare * n;
        return Mostly(lowA, highA) && Mostly(lowB, highB);
    }

    /// <summary>Unordered pairs, each written once with the ordinally smaller key first.</summary>
    private static IEnumerable<(ModificationSite, ModificationSite)> Pairs(IReadOnlyList<ModificationSite> sites)
    {
        var sorted = sites.OrderBy(s => s.Key, StringComparer.Ordinal).ToList();
        for (int i = 0; i < sorted.Count; i++)
            for (int j = i + 1; j < sorted.Count; j++)
                yield return (sorted[i], sorted[j]);
    }

    /// <summary>Median of a sorted array after removing one occurrence of each value in <paramref name="drop"/>.</summary>
    private static double MedianExcluding(double[] sorted, List<double> drop)
    {
        var removed = new List<int>(drop.Count);
        foreach (var d in drop)
        {
            int i = Array.BinarySearch(sorted, d);
            if (i < 0) continue;
            while (i > 0 && sorted[i - 1] == d) i--;
            while (removed.Contains(i)) i++;
            removed.Add(i);
        }
        removed.Sort();
        int m = sorted.Length - removed.Count;
        if (m <= 0) return double.NaN;
        double At(int k)
        {
            foreach (var r in removed) if (r <= k) k++;
            return sorted[k];
        }
        return m % 2 == 1 ? At(m / 2) : (At(m / 2 - 1) + At(m / 2)) / 2;
    }

    private static double Median(List<double> v)
    {
        if (v.Count == 0) return double.NaN;
        v.Sort();
        return v.Count % 2 == 1 ? v[v.Count / 2] : (v[v.Count / 2 - 1] + v[v.Count / 2]) / 2;
    }
}
