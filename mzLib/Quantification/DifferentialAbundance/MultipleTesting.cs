using System;
using System.Collections.Generic;
using System.Linq;

namespace Quantification.DifferentialAbundance
{
    /// <summary>Multiple-testing adjustment of p-values.</summary>
    /// <remarks>
    /// Not to be confused with the target-decoy q-values elsewhere in mzLib, which estimate the false
    /// discovery rate among identifications from decoy counts. These adjust a family of p-values.
    /// </remarks>
    public static class MultipleTesting
    {
        /// <summary>
        /// Benjamini–Hochberg (1995) step-up adjustment: for the i-th smallest of m p-values,
        /// min over j ≥ i of p_(j)·m/j, capped at 1. Non-finite entries are left NaN and are not counted in m,
        /// so the family is exactly the set of features that were tested.
        /// </summary>
        /// <returns>Adjusted values in the input order.</returns>
        public static double[] BenjaminiHochberg(IReadOnlyList<double> pValues)
        {
            ArgumentNullException.ThrowIfNull(pValues);
            var result = Enumerable.Repeat(double.NaN, pValues.Count).ToArray();
            var order = Enumerable.Range(0, pValues.Count)
                .Where(i => double.IsFinite(pValues[i]))
                .OrderBy(i => pValues[i]).ThenBy(i => i)
                .ToArray();
            int m = order.Length;
            double running = 1.0;
            for (int rank = m; rank >= 1; rank--)
            {
                int i = order[rank - 1];
                if (pValues[i] < 0 || pValues[i] > 1)
                    throw new ArgumentOutOfRangeException(nameof(pValues), pValues[i], $"p-value at index {i} is outside [0, 1].");
                running = Math.Min(running, pValues[i] * m / rank);
                result[i] = running;
            }
            return result;
        }
    }
}
