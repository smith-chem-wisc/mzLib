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
    /// <see cref="RunLevelRemoved"/>.
    /// </summary>
    public required double Statistic { get; init; }
    /// <summary>Type A: two-sided Spearman p (method in <see cref="SpearmanMethod"/>). Type P: NaN (no test in this version).</summary>
    public required double PValue { get; init; }
    /// <summary>Benjamini-Hochberg adjusted p within <see cref="FdrFamily"/>. NaN when not tested or excluded.</summary>
    public double Q { get; internal set; } = double.NaN;
    /// <summary>The family the adjustment ran over: result type, then same- or different-protein.</summary>
    public string FdrFamily => $"{ResultType}:{(SameProtein ? "intra" : "inter")}";
    /// <summary>Type P: runs where a doubly modified peptidoform was identified. Type A: runs where both occupancies were quantified.</summary>
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
    /// Type P: every pair of sites carried together by at least one peptidoform, with the runs it was
    /// identified in and its median co-occupancy. The observations are validated as in
    /// <see cref="SiteOccupancyCalculator.Calculate"/>.
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
    /// quantified. Floors and undetected runs are not values and are left out, and so, by default, are ceilings
    /// (runs where a site was seen only modified): two unrelated sites both read 1 in low-load runs, which would
    /// correlate them through detection. BH within each family, with overlapping pairs excluded from the family.
    /// </summary>
    /// <param name="occupancy">Output of <see cref="SiteOccupancyCalculator.Calculate"/> for one scope.</param>
    /// <param name="observations">The same observations, used to decide which pairs overlap.</param>
    /// <param name="minQuantifiedFraction">The pair-scale guard of D4 (default 0.7).</param>
    /// <param name="excludeCeiling">
    /// Leave out a run where the site was seen only modified (occupancy 1 as a ceiling,
    /// <see cref="SiteRunOccupancy.UnmodifiedQuantified"/> false), as <see cref="SiteTraitOptions.ExcludeCeiling"/>
    /// does for traits. Default true.
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
    public static IReadOnlyList<PtmPair> CoVarying(IReadOnlyList<SiteRunOccupancy> occupancy,
        IEnumerable<PeptidoformObservation> observations, double minQuantifiedFraction = DefaultMinQuantifiedFraction,
        bool excludeCeiling = true, bool removeRunLevel = false, Func<ModificationSite, string>? runLevelGroup = null)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        ArgumentNullException.ThrowIfNull(observations);
        if (!(minQuantifiedFraction > 0 && minQuantifiedFraction <= 1)) throw new ArgumentOutOfRangeException(nameof(minQuantifiedFraction));

        var runs = occupancy.Select(o => o.Run).Distinct().OrderBy(r => r, StringComparer.Ordinal).ToArray();
        var runIndex = runs.Select((r, i) => (r, i)).ToDictionary(t => t.r, t => t.i, StringComparer.Ordinal);
        var vectors = new Dictionary<ModificationSite, double[]>();
        foreach (var o in occupancy.Where(o => o.State == OccupancyState.Quantified && (!excludeCeiling || o.UnmodifiedQuantified)))
        {
            if (!vectors.TryGetValue(o.Site, out var v))
                vectors[o.Site] = v = Enumerable.Repeat(double.NaN, runs.Length).ToArray();
            v[runIndex[o.Run]] = o.Fraction;
        }
        int needed = (int)Math.Ceiling(minQuantifiedFraction * runs.Length);
        var sites = vectors.Where(kv => kv.Value.Count(double.IsFinite) >= needed).Select(kv => kv.Key)
            .OrderBy(s => s.Key, StringComparer.Ordinal).ToList();
        // Per group and run, the sorted values of the group's tested sites (the ones that enter pairs) with a value there.
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

        var result = new List<PtmPair>();
        foreach (var (a, b) in Pairs(sites))
        {
            SpearmanResult r;
            if (!removeRunLevel)
                r = SpearmanCorrelation.Correlate(vectors[a], vectors[b]);
            else
            {
                var groups = new[] { GroupOf(a), GroupOf(b) }.Distinct(StringComparer.Ordinal);
                r = SpearmanCorrelation.PartialCorrelate(vectors[a], vectors[b], groups.Select(g => (IReadOnlyList<double>)Level(g, a, b)).ToList());
            }
            bool tooFew = removeRunLevel && r.N < MinRunsForRunLevel;
            result.Add(new PtmPair
            {
                ResultType = PairResultType.A, SiteA = a, SiteB = b, Overlapping = Overlap(a, b),
                Statistic = tooFew ? double.NaN : r.Rho, PValue = tooFew ? double.NaN : r.PValue, N = r.N, StatisticN = r.N,
                SpearmanMethod = tooFew ? SpearmanPValueMethod.NotEstimable : r.Method,
                RunLevelRemoved = removeRunLevel,
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
