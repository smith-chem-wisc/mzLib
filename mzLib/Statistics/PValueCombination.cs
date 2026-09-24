using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;

namespace Statistics
{
    /// <summary>Which rule combined the p-values.</summary>
    public enum PValueCombinationMethod
    {
        /// <summary>Fisher (1925): X = −2 Σ ln pᵢ, referred to χ² with 2k degrees of freedom.</summary>
        Fisher,
        /// <summary>Stouffer (1949), optionally weighted: Z = Σ wᵢ zᵢ / √(Σ wᵢ²), zᵢ = Φ⁻¹(1 − pᵢ), referred to N(0, 1).</summary>
        Stouffer,
    }

    /// <summary>One p-value combined from several independent tests of the same hypothesis.</summary>
    public sealed class CombinedPValue
    {
        internal CombinedPValue() { }

        /// <summary>The rule used. Named on the result because the two rules answer different questions.</summary>
        public PValueCombinationMethod Method { get; internal set; }
        /// <summary>Number of finite p-values combined. Non-finite inputs are omitted and not counted.</summary>
        public int Studies { get; internal set; }
        /// <summary>Fisher's X (χ² with 2k df) or Stouffer's Z (standard normal). NaN when nothing was combined.</summary>
        public double Statistic { get; internal set; }
        /// <summary>Combined p-value, upper tail of the statistic's reference distribution. NaN when nothing was combined.</summary>
        public double PValue { get; internal set; }
    }

    /// <summary>
    /// Combination of p-values from independent tests (e.g. one per dataset) into one p-value.
    /// </summary>
    /// <remarks>
    /// Both rules assume the tests are INDEPENDENT and test the same one-sided hypothesis. Fisher is
    /// sensitive to one very small p-value; Stouffer weighs every study, and with weights (e.g. √n) it
    /// lets larger studies count more. For two-sided tests whose effects may point in different
    /// directions, convert each to a one-sided p-value with <see cref="OneSided"/> first, so that
    /// opposite effects cancel instead of reinforcing each other. Neither rule estimates an effect size:
    /// for that, pool estimates with <see cref="RandomEffectsMeta"/>.
    /// </remarks>
    public static class PValueCombination
    {
        /// <summary>Fisher's method: X = −2 Σ ln pᵢ ~ χ²(2k) under the joint null.</summary>
        /// <param name="pValues">One p-value per test, in [0, 1]. Non-finite values (tests that could not be run) are omitted.</param>
        public static CombinedPValue Fisher(IReadOnlyList<double> pValues)
        {
            var p = Finite(pValues);
            var result = new CombinedPValue { Method = PValueCombinationMethod.Fisher, Studies = p.Length,
                Statistic = double.NaN, PValue = double.NaN };
            if (p.Length == 0) return result;
            double x = -2 * p.Sum(Math.Log);
            result.Statistic = x;
            // Upper tail of χ²(2k) is Q(k, x/2), computed directly so a tiny p-value keeps its digits.
            result.PValue = double.IsPositiveInfinity(x) ? 0 : SpecialFunctions.GammaUpperRegularized(p.Length, x / 2);
            return result;
        }

        /// <summary>Stouffer's method: Z = Σ wᵢ zᵢ / √(Σ wᵢ²), zᵢ = Φ⁻¹(1 − pᵢ), Z ~ N(0, 1) under the joint null.</summary>
        /// <param name="pValues">One ONE-SIDED p-value per test, in [0, 1]. Non-finite values are omitted with their weight.</param>
        /// <param name="weights">Optional weight per test (finite, positive), e.g. √n. Unweighted when null.</param>
        public static CombinedPValue Stouffer(IReadOnlyList<double> pValues, IReadOnlyList<double>? weights = null)
        {
            ArgumentNullException.ThrowIfNull(pValues);
            if (weights != null && weights.Count != pValues.Count)
                throw new ArgumentException($"{weights.Count} weights for {pValues.Count} p-values.", nameof(weights));
            Finite(pValues);
            double num = 0, den = 0;
            int k = 0;
            for (int i = 0; i < pValues.Count; i++)
            {
                if (!double.IsFinite(pValues[i])) continue;
                double w = weights?[i] ?? 1.0;
                if (!(w > 0) || !double.IsFinite(w))
                    throw new ArgumentException($"Weight {i} must be finite and positive.", nameof(weights));
                num += w * Normal.InvCDF(0, 1, 1 - pValues[i]);
                den += w * w;
                k++;
            }
            var result = new CombinedPValue { Method = PValueCombinationMethod.Stouffer, Studies = k,
                Statistic = double.NaN, PValue = double.NaN };
            if (k == 0) return result;
            double z = num / Math.Sqrt(den);
            result.Statistic = z;
            result.PValue = double.IsNaN(z) ? double.NaN : Normal.CDF(0, 1, -z);
            return result;
        }

        /// <summary>
        /// Converts a two-sided p-value into the one-sided p-value for "effect &gt; 0", using the sign of the
        /// estimated effect: p/2 when the effect is positive, 1 − p/2 when negative, and NaN when the effect is
        /// zero or not finite (no direction to test).
        /// </summary>
        public static double OneSided(double twoSidedP, double effect)
        {
            if (!double.IsFinite(twoSidedP) || !double.IsFinite(effect) || effect == 0) return double.NaN;
            if (twoSidedP < 0 || twoSidedP > 1)
                throw new ArgumentOutOfRangeException(nameof(twoSidedP), twoSidedP, "p-value is outside [0, 1].");
            return effect > 0 ? twoSidedP / 2 : 1 - twoSidedP / 2;
        }

        private static double[] Finite(IReadOnlyList<double> pValues)
        {
            ArgumentNullException.ThrowIfNull(pValues);
            var p = pValues.Where(double.IsFinite).ToArray();
            for (int i = 0; i < p.Length; i++)
                if (p[i] < 0 || p[i] > 1)
                    throw new ArgumentOutOfRangeException(nameof(pValues), p[i], "A p-value is outside [0, 1].");
            return p;
        }
    }
}
