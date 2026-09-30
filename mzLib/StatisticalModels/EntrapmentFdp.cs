using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>A discovery candidate: its score (higher is better) and the q-value the search reported for it.</summary>
    public readonly record struct ScoredIdentification(double Score, double QValue);

    /// <summary>
    /// An entrapment discovery candidate with its paired original target.
    /// </summary>
    /// <param name="Score">Score, higher is better, on the same scale as the targets'.</param>
    /// <param name="QValue">The q-value the search reported.</param>
    /// <param name="Partner">
    /// The one original target this entrapment was made from, if the search identified it. Null when it was not
    /// identified, cannot be resolved to one target, or does not exist (foreign-species entrapment).
    /// Every null partner counts as scoring below the cutoff.
    /// </param>
    /// <param name="IsForeign">
    /// True for entrapment with no target partner by construction (foreign species). Controls only whether
    /// the paired estimate applies at all; see <see cref="EntrapmentFdp.Sweep"/>.
    /// </param>
    public readonly record struct ScoredEntrapment(double Score, double QValue, ScoredIdentification? Partner, bool IsForeign = false);

    /// <summary>The entrapment FDP estimates at one q-value threshold.</summary>
    /// <param name="Paired">Null where the paired estimator does not apply (see <see cref="EntrapmentFdp.Sweep"/>).</param>
    public sealed record EntrapmentFdpPoint(
        double QValueThreshold, int TargetCount, int EntrapmentCount, double LowerBound, double Combined, double? Paired);

    /// <summary>What an entrapment generator's exclusion means for the estimate.</summary>
    public enum EntrapmentExclusionKind
    {
        /// <summary>
        /// The peptide is really a target peptide sitting in the entrapment database, so matching it is a correct
        /// identification. Drop it from the entrapment count.
        /// </summary>
        NotReallyEntrapment,

        /// <summary>
        /// A target peptide that cannot be traced to one target. It stays a target discovery but must never serve as
        /// a pairing partner. Dropping it from the target count would bias the estimate the other way.
        /// </summary>
        Unpairable,
    }

    /// <summary>
    /// Entrapment estimates of the false discovery proportion among a search's discoveries, following Wen et al.
    /// 2025, Nat. Methods 22:1454 (the method FDRBench implements).
    /// </summary>
    /// <remarks>
    /// <para>
    /// N_T and N_E are the original-target and entrapment discoveries at a q-value threshold, and r is the
    /// effective ratio of entrapment to target database size:
    /// <list type="bullet">
    /// <item>lower bound (eq. 2): N_E / (N_T + N_E). It can only show that a search fails to control the FDR.</item>
    /// <item>combined (eq. 1): N_E (1 + 1/r) / (N_T + N_E). An upper bound; valid evidence of control.</item>
    /// <item>paired (eq. 4): (N_E + N_{E≥s&gt;T} + 2 N_{E&gt;T≥s}) / (N_T + N_E). A tighter upper bound that needs
    /// each target paired with exactly one entrapment, so r = 1.</item>
    /// </list>
    /// The paper's "sample" estimator, N_E / (r N_T), is invalid in both directions and is deliberately absent.
    /// </para>
    /// <para>
    /// Callers supply one entry per unit being estimated (precursor = peptide and charge, peptide, or protein group),
    /// with decoys removed and each unit keeping its best q-value. The estimators do not deduplicate.
    /// </para>
    /// </remarks>
    public static class EntrapmentFdp
    {
        /// <summary>0.001, 0.002, …, 0.100. Each is i/1000, so 1% is exactly 0.01.</summary>
        public static IReadOnlyList<double> DefaultQValueThresholds { get; } =
            Enumerable.Range(1, 100).Select(i => i / 1000.0).ToArray();

        /// <summary>Eq. 2: N_E / (N_T + N_E).</summary>
        public static double LowerBound(int targetCount, int entrapmentCount) =>
            (double)entrapmentCount / (targetCount + entrapmentCount);

        /// <summary>Eq. 1: N_E (1 + 1/r) / (N_T + N_E).</summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is not a positive finite number.</exception>
        public static double Combined(int targetCount, int entrapmentCount, double r)
        {
            ValidateRatio(r);
            return entrapmentCount * (1 + 1 / r) / (targetCount + entrapmentCount);
        }

        /// <summary>
        /// Eq. 4: (N_E + N_{E≥s&gt;T} + 2 N_{E&gt;T≥s}) / (N_T + N_E). Not capped at 1, as in the paper and FDRBench.
        /// </summary>
        /// <param name="partnerBelowCutoff">Discovered entrapments whose partner was not discovered.</param>
        /// <param name="partnerDiscoveredButLower">Discovered entrapments whose partner was discovered but scored lower.</param>
        public static double Paired(int targetCount, int entrapmentCount, int partnerBelowCutoff, int partnerDiscoveredButLower) =>
            (double)(entrapmentCount + partnerBelowCutoff + 2 * partnerDiscoveredButLower) / (targetCount + entrapmentCount);

        /// <summary>
        /// All three estimates at each q-value threshold. A candidate is discovered when its q-value is at or below
        /// the threshold. Thresholds with no discoveries are left out, since every estimate there is 0/0.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Paired is null unless r = 1, and null when every entrapment is foreign. With no partner anywhere it collapses
        /// to 2 N_E / (N_T + N_E), which ignores r and exceeds combined for any r above 1.
        /// </para>
        /// <para>
        /// A partner counts as discovered by its own q-value, not by its score. "Scored lower" is strict, so a tied
        /// partner adds nothing beyond the entrapment itself.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Any argument is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is not a positive finite number.</exception>
        public static IReadOnlyList<EntrapmentFdpPoint> Sweep(
            IReadOnlyCollection<ScoredIdentification> targets,
            IReadOnlyCollection<ScoredEntrapment> entrapments,
            double r,
            IReadOnlyList<double> qValueThresholds)
        {
            ArgumentNullException.ThrowIfNull(targets);
            ArgumentNullException.ThrowIfNull(entrapments);
            ArgumentNullException.ThrowIfNull(qValueThresholds);
            ValidateRatio(r);

            bool pairedApplies = r == 1.0 && (entrapments.Count == 0 || entrapments.Any(entrapment => !entrapment.IsForeign));

            var points = new List<EntrapmentFdpPoint>(qValueThresholds.Count);
            foreach (double threshold in qValueThresholds)
            {
                int targetCount = targets.Count(target => target.QValue <= threshold);
                int entrapmentCount = 0, partnerBelowCutoff = 0, partnerDiscoveredButLower = 0;
                foreach (var entrapment in entrapments)
                {
                    if (entrapment.QValue > threshold)
                        continue;
                    entrapmentCount++;

                    if (entrapment.Partner is not { } partner || partner.QValue > threshold)
                        partnerBelowCutoff++;
                    else if (partner.Score < entrapment.Score)
                        partnerDiscoveredButLower++;
                }

                if (targetCount + entrapmentCount == 0)
                    continue;

                points.Add(new EntrapmentFdpPoint(
                    threshold,
                    targetCount,
                    entrapmentCount,
                    LowerBound(targetCount, entrapmentCount),
                    Combined(targetCount, entrapmentCount, r),
                    pairedApplies ? Paired(targetCount, entrapmentCount, partnerBelowCutoff, partnerDiscoveredButLower) : null));
            }
            return points;
        }

        /// <summary>
        /// Maps an entrapment generator's exclusion reason to what it means for the estimate.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// The reason is not recognised. The two kinds correct the estimate in opposite directions, so a new reason
        /// must be classified deliberately rather than defaulted.
        /// </exception>
        public static EntrapmentExclusionKind ClassifyExclusion(string reason) => reason switch
        {
            "unrepairableRunCollision" or "sharedWithTarget" or "initiatorMethionineCollision" => EntrapmentExclusionKind.NotReallyEntrapment,
            "ambiguous" => EntrapmentExclusionKind.Unpairable,
            _ => throw new ArgumentException(
                $"Unrecognised entrapment exclusion reason '{reason}'. It must be classified as removing a real target from the " +
                "entrapment count or as an unpairable target; the two bias the estimate in opposite directions.", nameof(reason)),
        };

        private static void ValidateRatio(double r)
        {
            if (!double.IsFinite(r) || r <= 0)
                throw new ArgumentOutOfRangeException(nameof(r), r, "The entrapment-to-target ratio must be a positive finite number.");
        }
    }
}
