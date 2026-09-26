using System;
using System.Collections.Generic;
using System.Linq;
using Statistics;

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
    /// Some peptidoform covers both positions. For type A the two occupancies then share peptides and can be
    /// anti-correlated by construction, so overlapping pairs are reported but excluded from the A family's
    /// multiple-testing adjustment (q is NaN).
    /// </summary>
    public required bool Overlapping { get; init; }
    /// <summary>
    /// Type P: median over runs of the co-occupancy, the intensity of peptidoforms carrying both modifications
    /// over the intensity of peptidoforms covering both positions (NaN when never quantified).
    /// Type A: Spearman's ρ of the two occupancies across runs.
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
    /// <summary>Type A: how the Spearman p-value was computed. Type P: null.</summary>
    public SpearmanPValueMethod? SpearmanMethod { get; init; }
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
    /// Type P: every pair of sites carried together by at least one peptidoform, with the runs it was
    /// identified in and its median co-occupancy.
    /// </summary>
    public static IReadOnlyList<PtmPair> Physical(IEnumerable<PeptidoformObservation> observations,
        Func<string, bool>? includeModification = null)
    {
        ArgumentNullException.ThrowIfNull(observations);
        var obs = observations.Select(o =>
        {
            string baseSeq = MzLibUtil.ClassExtensions.GetBaseSequenceFromFullSequence(o.FullSequence);
            return (o, sites: SiteOccupancyCalculator.SitesOf(o, baseSeq, includeModification));
        }).ToList();

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
                Statistic = Median(coOccupancy), PValue = double.NaN, N = runs.Count,
            });
        }
        return result;
    }

    /// <summary>
    /// Type A: Spearman correlation across runs between the occupancies of every pair of sites quantified in
    /// at least <paramref name="minQuantifiedFraction"/> of the scope's runs, on the runs where both are
    /// quantified. Floors and undetected runs are not values and are left out. BH within each family, with
    /// overlapping pairs excluded from the family.
    /// </summary>
    /// <param name="occupancy">Output of <see cref="SiteOccupancyCalculator.Calculate"/> for one scope.</param>
    /// <param name="observations">The same observations, used to decide which pairs overlap.</param>
    /// <param name="minQuantifiedFraction">The pair-scale guard of D4 (default 0.7).</param>
    public static IReadOnlyList<PtmPair> CoVarying(IReadOnlyList<SiteRunOccupancy> occupancy,
        IEnumerable<PeptidoformObservation> observations, double minQuantifiedFraction = DefaultMinQuantifiedFraction)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        ArgumentNullException.ThrowIfNull(observations);
        if (!(minQuantifiedFraction > 0 && minQuantifiedFraction <= 1)) throw new ArgumentOutOfRangeException(nameof(minQuantifiedFraction));

        var runs = occupancy.Select(o => o.Run).Distinct().OrderBy(r => r, StringComparer.Ordinal).ToArray();
        var runIndex = runs.Select((r, i) => (r, i)).ToDictionary(t => t.r, t => t.i, StringComparer.Ordinal);
        var vectors = new Dictionary<ModificationSite, double[]>();
        foreach (var o in occupancy.Where(o => o.State == OccupancyState.Quantified))
        {
            if (!vectors.TryGetValue(o.Site, out var v))
                vectors[o.Site] = v = Enumerable.Repeat(double.NaN, runs.Length).ToArray();
            v[runIndex[o.Run]] = o.Fraction;
        }
        int needed = (int)Math.Ceiling(minQuantifiedFraction * runs.Length);
        var sites = vectors.Where(kv => kv.Value.Count(double.IsFinite) >= needed).Select(kv => kv.Key)
            .OrderBy(s => s.Key, StringComparer.Ordinal).ToList();

        // Spans covered by one peptidoform, per protein, to flag overlapping pairs.
        var spans = observations.GroupBy(o => o.ProteinAccession, StringComparer.Ordinal)
            .ToDictionary(g => g.Key, g => g.Select(o => (o.StartResidue, o.EndResidue)).Distinct().ToList(), StringComparer.Ordinal);
        bool Overlap(ModificationSite a, ModificationSite b) =>
            a.ProteinAccession == b.ProteinAccession && spans.TryGetValue(a.ProteinAccession, out var s)
            && s.Any(sp => Math.Min(a.Position, b.Position) >= sp.StartResidue && Math.Max(a.Position, b.Position) <= sp.EndResidue);

        var result = new List<PtmPair>();
        foreach (var (a, b) in Pairs(sites))
        {
            var r = SpearmanCorrelation.Correlate(vectors[a], vectors[b]);
            result.Add(new PtmPair
            {
                ResultType = PairResultType.A, SiteA = a, SiteB = b, Overlapping = Overlap(a, b),
                Statistic = r.Rho, PValue = r.PValue, N = r.N, SpearmanMethod = r.Method,
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

    private static double Median(List<double> v)
    {
        if (v.Count == 0) return double.NaN;
        v.Sort();
        return v.Count % 2 == 1 ? v[v.Count / 2] : (v[v.Count / 2 - 1] + v[v.Count / 2]) / 2;
    }
}
