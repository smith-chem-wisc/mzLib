using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using StatisticalModels;

namespace Quantification.PtmQtl;

/// <summary>One PTM pair pooled across scopes (datasets).</summary>
public sealed record GlobalPtmPair
{
    /// <summary>Evidence type.</summary>
    public required PairResultType ResultType { get; init; }
    /// <summary>Canonical key of the first site (ordinal order).</summary>
    public required string SiteA { get; init; }
    /// <summary>Canonical key of the second site.</summary>
    public required string SiteB { get; init; }
    /// <summary>Protein of the first site.</summary>
    public required string ProteinA { get; init; }
    /// <summary>Protein of the second site.</summary>
    public required string ProteinB { get; init; }
    /// <summary>Both sites on one protein accession.</summary>
    public bool SameProtein => ProteinA == ProteinB;
    /// <summary>Scopes contributing a row for this pair.</summary>
    public required IReadOnlyList<string> Scopes { get; init; }
    /// <summary>
    /// Type A: scopes whose ρ has the sign of the combined statistic. Type P: every scope (being
    /// seen on one molecule has no sign).
    /// </summary>
    public required int ScopesAgreeing { get; init; }
    /// <summary>Type A: median ρ across scopes. Type P: median across scopes of the median co-occupancy.</summary>
    public required double Statistic { get; init; }
    /// <summary>Type A: Stouffer's Z of the signed one-sided p-values, weighted by √n. Type P: NaN.</summary>
    public required double CombinedZ { get; init; }
    /// <summary>Type A: two-sided p of <see cref="CombinedZ"/>. Type P: NaN (not tested in this version).</summary>
    public required double PValue { get; init; }
    /// <summary>Benjamini-Hochberg within <see cref="FdrFamily"/>, over the pooled pairs. NaN when untested.</summary>
    public double Q { get; internal set; } = double.NaN;
    /// <summary>Result type, then same- or different-protein.</summary>
    public string FdrFamily => $"{ResultType}:{(SameProtein ? "intra" : "inter")}";
    /// <summary>How the scopes were combined.</summary>
    public string Method => ResultType == PairResultType.A
        ? "Stouffer, signed one-sided p, weights sqrt(n)"
        : "recurrence count";
}

/// <summary>
/// Pools per-scope PTM pairs (<see cref="PtmPairEngine"/>) across scopes. Sites from different scopes are
/// matched through a caller-supplied canonical key, so that two engine names for one chemistry (e.g. a
/// UniProt and a GPTMD name for phosphoserine) are one site.
/// </summary>
/// <remarks>
/// Type A: each scope's two-sided Spearman p is made one-sided in the direction of its ρ, as z = sign(ρ) Φ⁻¹(1 − p/2)
/// (the z of <see cref="PValueCombination.OneSided"/>, computed in the tail so that both signs keep their digits),
/// and the scopes are combined by weighted Stouffer with weights √n. Opposite correlations in two scopes therefore
/// cancel instead of reinforcing each other. A scope with ρ = 0 has no direction and is left out. A scope with
/// p = 0 (Spearman's asymptotic p at |ρ| = 1) is taken at p = 2 × double.Epsilon, so its z is finite (about 38)
/// and cannot make the combined Z infinite or NaN. The combined two-sided p is 2 Φ(−|Z|). BH runs within each family (same or
/// different protein) over the pooled pairs. Only non-overlapping type-A rows are pooled (see
/// <see cref="PtmPair.Overlapping"/>). Type P is counted, not tested. Scopes are assumed independent, which
/// holds for separate datasets and not for re-analyses of one sample set.
/// </remarks>
public static class GlobalPairEngine
{
    /// <summary>Pools pairs seen in at least <paramref name="minScopes"/> scopes.</summary>
    /// <param name="pairs">Per-scope pairs, each with its scope (dataset) identifier.</param>
    /// <param name="canonicalKey">Maps a site to the key that identifies it across scopes.</param>
    /// <param name="minScopes">Minimum scopes a pair must appear in (default 2).</param>
    public static IReadOnlyList<GlobalPtmPair> Combine(IEnumerable<(string Scope, PtmPair Pair)> pairs,
        Func<ModificationSite, string> canonicalKey, int minScopes = 2)
    {
        ArgumentNullException.ThrowIfNull(pairs);
        ArgumentNullException.ThrowIfNull(canonicalKey);
        if (minScopes < 1) throw new ArgumentOutOfRangeException(nameof(minScopes));

        // One row per (type, canonical pair, scope). Two engine names at one site can give one scope two
        // rows for one canonical pair; the row with the most runs is kept, then the smaller p.
        var rows = new Dictionary<(PairResultType, string, string), Dictionary<string, (PtmPair pair, string protA, string protB)>>();
        foreach (var (scope, p) in pairs)
        {
            if (p.ResultType == PairResultType.A && (p.Overlapping || double.IsNaN(p.PValue))) continue;
            string ka = canonicalKey(p.SiteA), kb = canonicalKey(p.SiteB);
            if (ka == kb) continue;
            var (a, b, pa, pb) = string.CompareOrdinal(ka, kb) < 0
                ? (ka, kb, p.SiteA.ProteinAccession, p.SiteB.ProteinAccession)
                : (kb, ka, p.SiteB.ProteinAccession, p.SiteA.ProteinAccession);
            var key = (p.ResultType, a, b);
            if (!rows.TryGetValue(key, out var byScope)) rows[key] = byScope = new(StringComparer.Ordinal);
            if (!byScope.TryGetValue(scope, out var existing)
                || p.N > existing.pair.N || (p.N == existing.pair.N && p.PValue < existing.pair.PValue))
                byScope[scope] = (p, pa, pb);
        }

        var result = new List<GlobalPtmPair>();
        foreach (var ((type, a, b), byScope) in rows.OrderBy(r => r.Key.Item1).ThenBy(r => r.Key.Item2, StringComparer.Ordinal)
                                                     .ThenBy(r => r.Key.Item3, StringComparer.Ordinal))
        {
            if (byScope.Count < minScopes) continue;
            var scopes = byScope.Keys.OrderBy(s => s, StringComparer.Ordinal).ToList();
            var members = scopes.Select(s => byScope[s]).ToList();
            var first = members[0];
            if (type == PairResultType.A)
            {
                double num = 0, den = 0;
                foreach (var (pair, _, _) in members)
                {
                    double zi = SignedZ(pair.PValue, pair.Statistic);
                    if (double.IsNaN(zi)) continue;
                    double w = Math.Sqrt(pair.N);
                    num += w * zi;
                    den += w * w;
                }
                double z = den > 0 ? num / Math.Sqrt(den) : double.NaN;
                double two = double.IsNaN(z) ? double.NaN : Math.Min(1, 2 * Normal.CDF(0, 1, -Math.Abs(z)));
                result.Add(new GlobalPtmPair
                {
                    ResultType = type, SiteA = a, SiteB = b, ProteinA = first.protA, ProteinB = first.protB, Scopes = scopes,
                    ScopesAgreeing = double.IsFinite(z) && z != 0
                        ? members.Count(m => double.IsFinite(m.pair.Statistic) && Math.Sign(m.pair.Statistic) == Math.Sign(z))
                        : 0,
                    Statistic = Median(members.Select(m => m.pair.Statistic)), CombinedZ = z, PValue = two,
                });
            }
            else
            {
                result.Add(new GlobalPtmPair
                {
                    ResultType = type, SiteA = a, SiteB = b, ProteinA = first.protA, ProteinB = first.protB, Scopes = scopes,
                    ScopesAgreeing = scopes.Count, Statistic = Median(members.Select(m => m.pair.Statistic)),
                    CombinedZ = double.NaN, PValue = double.NaN,
                });
            }
        }
        foreach (var family in result.Where(r => r.ResultType == PairResultType.A && !double.IsNaN(r.PValue)).GroupBy(r => r.FdrFamily))
        {
            var members = family.ToList();
            var q = MultipleTesting.BenjaminiHochberg(members.Select(r => r.PValue).ToArray());
            for (int i = 0; i < members.Count; i++) members[i].Q = q[i];
        }
        return result;
    }

    /// <summary>
    /// sign(effect) × Φ⁻¹(1 − p/2) for a two-sided p, computed as −Φ⁻¹(p/2) so a small p keeps its digits for either
    /// sign; p/2 is floored at double.Epsilon. NaN when p or the effect is not finite or the effect is 0.
    /// </summary>
    private static double SignedZ(double twoSidedP, double effect)
    {
        if (!double.IsFinite(twoSidedP) || !double.IsFinite(effect) || effect == 0) return double.NaN;
        if (twoSidedP < 0 || twoSidedP > 1)
            throw new ArgumentOutOfRangeException(nameof(twoSidedP), twoSidedP, "p-value is outside [0, 1].");
        return Math.Sign(effect) * -Normal.InvCDF(0, 1, Math.Max(twoSidedP / 2, double.Epsilon));
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.Where(x => !double.IsNaN(x)).OrderBy(x => x).ToList();
        if (v.Count == 0) return double.NaN;
        return v.Count % 2 == 1 ? v[v.Count / 2] : (v[v.Count / 2 - 1] + v[v.Count / 2]) / 2;
    }
}
