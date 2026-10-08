using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
{
    /// <summary>
    /// Target-decoy q-values: (D + 1) / T down the ranking, capped at 1 and made monotone, the estimate for which
    /// target-decoy competition controls the FDR.
    /// </summary>
    public static class TargetDecoyQValues
    {
        /// <summary>
        /// One q-value per entry, in input order; higher scores rank first. Decoys receive the value at their rank. Among
        /// equal scores, decoys rank first (the conservative choice), then input order.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree or a score is not finite.</exception>
        public static double[] Compute(IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy)
        {
            ArgumentNullException.ThrowIfNull(scores);
            ArgumentNullException.ThrowIfNull(isDecoy);
            if (scores.Count != isDecoy.Count)
                throw new ArgumentException($"There are {scores.Count} scores but {isDecoy.Count} decoy labels.", nameof(isDecoy));
            if (!scores.All(double.IsFinite))
                throw new ArgumentException("Every score must be finite.", nameof(scores));

            // Score descending, decoys first at a tie, then index: the same total order as a stable LINQ sort, without its
            // per-comparison key machinery (the rescorer's seed step computes this for every feature and sign)
            double[] s = scores.ToArray();
            bool[] d = isDecoy.ToArray();
            int[] order = new int[s.Length];
            for (int i = 0; i < order.Length; i++)
                order[i] = i;
            Array.Sort(order, (a, b) =>
            {
                int c = s[b].CompareTo(s[a]);
                if (c != 0) return c;
                c = d[b].CompareTo(d[a]);
                return c != 0 ? c : a.CompareTo(b);
            });
            var ranked = new double[order.Length];
            int decoys = 0, targets = 0;
            for (int k = 0; k < order.Length; k++)
            {
                if (isDecoy[order[k]]) decoys++; else targets++;
                ranked[k] = targets == 0 ? 1 : Math.Min(1, (decoys + 1.0) / targets);
            }
            for (int k = order.Length - 2; k >= 0; k--)
                ranked[k] = Math.Min(ranked[k], ranked[k + 1]);

            var q = new double[order.Length];
            for (int k = 0; k < order.Length; k++)
                q[order[k]] = ranked[k];
            return q;
        }
    }

    /// <summary>Which entries a picked competition kept, and their q-values.</summary>
    public sealed class PickedResult
    {
        internal PickedResult(bool[] kept, double[] qValues)
        {
            Kept = kept;
            QValues = qValues;
        }

        /// <summary>True for the one entry that competed under each pair key, in input order.</summary>
        public bool[] Kept { get; }

        /// <summary>(D+1)/T over the kept entries, in input order; NaN for an entry that was not kept.</summary>
        public double[] QValues { get; }
    }

    /// <summary>
    /// The picked target-decoy competition (Savitski et al. 2015, Mol. Cell. Proteomics 14:2394). A target and its own decoy
    /// share a pair key, for example a protein accession with its decoy prefix removed. Only the better-scoring of each
    /// pair is kept, so a target and its decoy never both count, and q-values are computed over the kept entries alone.
    /// </summary>
    /// <remarks>
    /// The key is the caller's: this class knows nothing of accessions or decoy prefixes. An entry with a key of its own
    /// competes alone. A target tied with its decoy loses, the conservative choice. Entries that were not kept get no
    /// q-value; they are not rescued with a classic one.
    /// </remarks>
    public static class PickedTargetDecoy
    {
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, a key is null, or a score is not finite.</exception>
        public static PickedResult Compete(IReadOnlyList<string> pairKeys, IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy)
        {
            ArgumentNullException.ThrowIfNull(pairKeys);
            ArgumentNullException.ThrowIfNull(scores);
            ArgumentNullException.ThrowIfNull(isDecoy);
            if (scores.Count != pairKeys.Count)
                throw new ArgumentException($"There are {pairKeys.Count} pair keys but {scores.Count} scores.", nameof(scores));
            if (isDecoy.Count != pairKeys.Count)
                throw new ArgumentException($"There are {pairKeys.Count} pair keys but {isDecoy.Count} decoy labels.", nameof(isDecoy));
            if (pairKeys.Any(key => key is null))
                throw new ArgumentException("Every entry needs a pair key.", nameof(pairKeys));
            if (!scores.All(double.IsFinite))
                throw new ArgumentException("Every score must be finite.", nameof(scores));

            // The winner under each key: highest score, then decoy over target, then earliest input
            var winners = Enumerable.Range(0, pairKeys.Count)
                .GroupBy(i => pairKeys[i], StringComparer.Ordinal)
                .Select(group => group.OrderByDescending(i => scores[i]).ThenByDescending(i => isDecoy[i]).ThenBy(i => i).First())
                .Order()
                .ToArray();

            double[] winnerQ = TargetDecoyQValues.Compute(winners.Select(i => scores[i]).ToArray(), winners.Select(i => isDecoy[i]).ToArray());
            var kept = new bool[pairKeys.Count];
            var qValues = Enumerable.Repeat(double.NaN, pairKeys.Count).ToArray();
            for (int w = 0; w < winners.Length; w++)
            {
                kept[winners[w]] = true;
                qValues[winners[w]] = winnerQ[w];
            }
            return new PickedResult(kept, qValues);
        }
    }
}
