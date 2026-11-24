using System;
using System.Collections.Generic;
using System.Linq;

namespace StatisticalModels
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

        /// <summary>
        /// Benjamini–Bogomolov (2014) selective inference over families of hypotheses: when R of m families have been
        /// selected (for example by <see cref="BenjaminiHochberg"/> on one p per family), each selected family's p-values
        /// are adjusted by <see cref="BenjaminiHochberg"/> within the family and compared at q·R/m, not at q. That controls
        /// the expected average false discovery rate over the selected families, under the conditions Benjamini and
        /// Bogomolov give. Families that were not selected are not tested.
        /// </summary>
        /// <param name="withinFamilyPValues">Each family's p-values. Every family considered for selection is listed, so m is
        /// their number. Non-finite entries are left NaN and not counted, as in <see cref="BenjaminiHochberg"/>.</param>
        /// <param name="selected">Whether each family was selected, in the same order.</param>
        /// <param name="q">The level, in (0, 1).</param>
        public static SelectiveAdjustment BenjaminiBogomolov(IReadOnlyList<IReadOnlyList<double>> withinFamilyPValues,
            IReadOnlyList<bool> selected, double q = 0.05)
        {
            ArgumentNullException.ThrowIfNull(withinFamilyPValues);
            ArgumentNullException.ThrowIfNull(selected);
            if (selected.Count != withinFamilyPValues.Count)
                throw new ArgumentException($"{selected.Count} selection flags for {withinFamilyPValues.Count} families.", nameof(selected));
            if (!(q > 0 && q < 1)) throw new ArgumentOutOfRangeException(nameof(q), q, "q must lie strictly between 0 and 1.");
            int m = withinFamilyPValues.Count, r = selected.Count(s => s);
            var adjusted = new double[m][];
            for (int f = 0; f < m; f++)
            {
                var family = withinFamilyPValues[f] ?? throw new ArgumentNullException(nameof(withinFamilyPValues), $"Family {f} is null.");
                adjusted[f] = selected[f] ? BenjaminiHochberg(family) : Enumerable.Repeat(double.NaN, family.Count).ToArray();
            }
            return new SelectiveAdjustment { Families = m, Selected = r, Level = r > 0 ? q * r / m : double.NaN, Adjusted = adjusted };
        }
    }

    /// <summary>The result of <see cref="MultipleTesting.BenjaminiBogomolov"/>.</summary>
    public sealed class SelectiveAdjustment
    {
        internal SelectiveAdjustment() { }

        /// <summary>m: the families considered for selection.</summary>
        public int Families { get; internal set; }
        /// <summary>R: the families selected.</summary>
        public int Selected { get; internal set; }
        /// <summary>q·R/m, the level a selected family's adjusted p-values are compared at; NaN when none was selected.</summary>
        public double Level { get; internal set; }
        /// <summary>Per family, in input order: its within-family Benjamini–Hochberg adjusted p-values in its input order
        /// when it was selected; NaN in every entry when it was not.</summary>
        public IReadOnlyList<double[]> Adjusted { get; internal set; } = Array.Empty<double[]>();
    }
}
