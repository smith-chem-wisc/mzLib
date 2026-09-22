using System;
using MathNet.Numerics;

namespace Quantification.DifferentialAbundance
{
    /// <summary>
    /// The trigamma function and its inverse, which empirical-Bayes variance moderation needs and
    /// MathNet.Numerics 5.0.0 does not provide (it has <see cref="SpecialFunctions.DiGamma"/> only).
    /// </summary>
    /// <remarks>
    /// Both polygammas use the upward recurrence to move the argument above 10 and then the standard
    /// asymptotic series, which is accurate to about 1e-15 there. They are defined for x &gt; 0 only,
    /// which is the only domain moderation calls them on (half a degrees-of-freedom value).
    /// </remarks>
    internal static class Polygamma
    {
        private const double AsymptoticThreshold = 10.0;
        private const int MaxNewtonIterations = 100;
        private const double NewtonRelativeTolerance = 1e-12;
        /// <summary>Below this, Trigamma(y) = x has y = 1/x to double precision.</summary>
        private const double SmallArgumentLimit = 1e-6;

        /// <summary>The trigamma function, the second derivative of ln Γ(x), for x &gt; 0.</summary>
        public static double Trigamma(double x)
        {
            if (!(x > 0)) throw new ArgumentOutOfRangeException(nameof(x), x, "Trigamma is implemented for x > 0 only.");
            double shift = 0;
            while (x < AsymptoticThreshold)
            {
                shift += 1.0 / (x * x);
                x += 1;
            }
            double r = 1.0 / x, r2 = r * r;
            // psi1(x) ~ 1/x + 1/(2x^2) + sum_k B_2k / x^(2k+1)
            double series = r + r2 / 2 + r * r2 * (1.0 / 6 + r2 * (-1.0 / 30 + r2 * (1.0 / 42 + r2 * (-1.0 / 30 + r2 * (5.0 / 66)))));
            return shift + series;
        }

        /// <summary>The tetragamma function, the third derivative of ln Γ(x), for x &gt; 0.</summary>
        public static double Tetragamma(double x)
        {
            if (!(x > 0)) throw new ArgumentOutOfRangeException(nameof(x), x, "Tetragamma is implemented for x > 0 only.");
            double shift = 0;
            while (x < AsymptoticThreshold)
            {
                shift -= 2.0 / (x * x * x);
                x += 1;
            }
            double r = 1.0 / x, r2 = r * r;
            // psi2(x) ~ -1/x^2 - 1/x^3 - sum_k (2k+1) B_2k / x^(2k+2)
            double series = -r2 - r * r2 - r2 * r2 * (0.5 + r2 * (-1.0 / 6 + r2 * (1.0 / 6 + r2 * (-3.0 / 10 + r2 * (5.0 / 6)))));
            return shift + series;
        }

        /// <summary>
        /// Solves Trigamma(y) = x for y, for x &gt; 0. Trigamma is strictly decreasing on (0, ∞), so the
        /// solution is unique.
        /// </summary>
        /// <remarks>
        /// Newton's method on 1/Trigamma(y), which is close to linear in y (≈ y - 1/2 for large y) and so
        /// converges without overshooting into y &lt;= 0. The two limits set the starting value: Trigamma(y) ≈ 1/y² as y → 0 and ≈ 1/y as y → ∞; only x &lt; 1e-6 is returned
        /// from the limit directly, where it is exact to double precision.
        /// </remarks>
        public static double TrigammaInverse(double x)
        {
            if (!(x > 0) || !double.IsFinite(x)) throw new ArgumentOutOfRangeException(nameof(x), x, "TrigammaInverse is defined for finite x > 0 only.");
            if (x < SmallArgumentLimit) return 1.0 / x;

            // Start near the solution from the matching limit, so Newton needs only a few steps at either end.
            double y = x > 1 ? 1.0 / Math.Sqrt(x) : 0.5 + 1.0 / x;
            for (int i = 0; i < MaxNewtonIterations; i++)
            {
                double tri = Trigamma(y);
                // f(y) = 1/tri - 1/x ; f'(y) = -tetra / tri^2 ; step = f / f'
                double step = tri * (1 - tri / x) / Tetragamma(y);
                y += step;
                if (Math.Abs(step) / y < NewtonRelativeTolerance) break;
            }
            return y;
        }
    }
}
