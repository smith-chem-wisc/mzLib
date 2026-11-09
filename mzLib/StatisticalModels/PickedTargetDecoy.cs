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
        /// One q-value per entry, in input order; higher scores rank first. Decoys receive the value at their rank. Equal
        /// scores are one rank: the whole tie block is counted before any of it gets its value, so a tied target never ranks
        /// ahead of a tied decoy (the conservative choice) and neither labels nor row order change any value.
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

            // Only read, so arrays are used as given: the seed step calls this for every feature and sign of every fold
            double[] s = scores as double[] ?? scores.ToArray();
            bool[] d = isDecoy as bool[] ?? isDecoy.ToArray();
            var q = new double[s.Length];
            Into(s, d, s.Length, new int[s.Length], new double[s.Length], q);
            return q;
        }

        /// <summary>
        /// <see cref="Compute"/> over the first <paramref name="count"/> entries, written into <paramref name="q"/>, with the
        /// caller's work buffers: for a caller that ranks many score vectors of one length. Every buffer holds at least
        /// <paramref name="count"/> entries; the scores must be finite.
        /// </summary>
        internal static void Into(double[] s, bool[] d, int count, int[] order, double[] ranked, double[] q)
        {
            // Score descending, decoys first at a tie, then index: a total order, so any sort gives this one permutation,
            // the same as a stable LINQ sort, without its per-comparison key machinery
            for (int i = 0; i < count; i++)
                order[i] = i;
            order.AsSpan(0, count).Sort((a, b) =>
            {
                int c = s[b].CompareTo(s[a]);
                if (c != 0) return c;
                c = d[b].CompareTo(d[a]);
                return c != 0 ? c : a.CompareTo(b);
            });
            int decoys = 0, targets = 0;
            for (int k = 0; k < count; )
            {
                int end = k;
                do
                {
                    if (d[order[end]]) decoys++; else targets++;
                    end++;
                }
                while (end < count && s[order[end]] == s[order[k]]);
                double value = targets == 0 ? 1 : Math.Min(1, (decoys + 1.0) / targets);
                for (; k < end; k++)
                    ranked[k] = value;
            }
            for (int k = count - 2; k >= 0; k--)
                ranked[k] = Math.Min(ranked[k], ranked[k + 1]);

            for (int k = 0; k < count; k++)
                q[order[k]] = ranked[k];
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
