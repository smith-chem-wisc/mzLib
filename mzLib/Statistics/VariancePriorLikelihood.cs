using System;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace Statistics
{
    /// <summary>
    /// Maximum marginal likelihood of the variance prior, with each feature's own residual degrees of freedom.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Under the model of <see cref="EmpiricalBayes"/>, s²_g | σ²_g ~ σ²_g χ²(d_g)/d_g and 1/σ²_g ~ χ²(d0)/(d0 s0²_g),
    /// so marginally s²_g / s0²_g ~ F(d_g, d0). The log-likelihood of (d0, s0²) is the sum of the log F
    /// densities, each at its OWN d_g. Unequal residual df therefore enter exactly, through the likelihood:
    /// a feature with fewer df has a flatter F density and so constrains the prior less, and no weights have
    /// to be chosen. When d0 = ∞ the prior is a point mass and s²_g / s0²_g ~ χ²(d_g)/d_g.
    /// </para>
    /// <para>
    /// log s0²_g = B_g γ, with B the intercept (no trend) or a natural cubic spline basis in the covariate.
    /// For fixed d0 the log-likelihood is concave in log s0², hence in γ, so γ is found by Newton's method
    /// with step halving. d0 is found by a grid over log d0 followed by Brent's method, and d0 = ∞ is kept
    /// when its likelihood is at least as high.
    /// </para>
    /// <para>
    /// This is an independent estimator, derived from the model, not a port of limma's
    /// <c>fitFDistUnequalDF1</c>. The two agree in what they estimate but not in procedure (limma adds a
    /// df-based down-weighting and a lowess trend), so their numbers differ.
    /// </para>
    /// </remarks>
    internal static class VariancePriorLikelihood
    {
        /// <summary>Range searched for log d0 before refinement: d0 from about 0.05 to 1.2 million.</summary>
        private const double LogDfLow = -3, LogDfHigh = 14, LogDfStep = 0.25;

        /// <summary>Fits (γ, d0). Returns the prior df and the fitted log s0² per feature.</summary>
        internal static (double d0, double[] logScale) Fit(double[] s2, double[] df, Matrix<double> basis)
        {
            int n = s2.Length;
            var logS2 = s2.Select(v => Math.Log(v)).ToArray();
            // Start γ at the least-squares fit of log s² on the basis; every inner fit then warm-starts from it.
            var start = basis.QR(MathNet.Numerics.LinearAlgebra.Factorization.QRMethod.Thin)
                             .Solve(Vector<double>.Build.DenseOfArray(logS2));

            double Neg(double u) => -Profile(s2, df, basis, Math.Exp(u), start).logLik;

            double bestU = LogDfLow, best = double.MaxValue;
            for (double u = LogDfLow; u <= LogDfHigh + 1e-12; u += LogDfStep)
            {
                double v = Neg(u);
                if (v < best) { best = v; bestU = u; }
            }
            double uStar = Brent(Neg, Math.Max(LogDfLow, bestU - LogDfStep), Math.Min(LogDfHigh, bestU + LogDfStep), 1e-10);
            var finite = Profile(s2, df, basis, Math.Exp(uStar), start);
            var infinite = Profile(s2, df, basis, double.PositiveInfinity, start);
            var chosen = infinite.logLik >= finite.logLik ? infinite : finite;
            double d0 = infinite.logLik >= finite.logLik ? double.PositiveInfinity : Math.Exp(uStar);
            return (d0, (basis * chosen.gamma).ToArray());
        }

        /// <summary>γ maximizing the log-likelihood at fixed d0, and that maximum.</summary>
        internal static (double logLik, Vector<double> gamma) Profile(double[] s2, double[] df, Matrix<double> basis,
            double d0, Vector<double> start)
        {
            var gamma = start.Clone();
            double ll = LogLik(s2, df, basis, d0, gamma);
            for (int iter = 0; iter < 100; iter++)
            {
                var (grad, hess) = Derivatives(s2, df, basis, d0, gamma);
                var step = (-hess).Cholesky().Solve(grad);
                // Concave objective: the full Newton step is almost always accepted; halve it if not.
                double t = 1, llNext;
                Vector<double> next;
                do { next = gamma + t * step; llNext = LogLik(s2, df, basis, d0, next); t /= 2; }
                while (llNext < ll && t > 1e-10);
                if (llNext < ll) break;
                double gain = llNext - ll;
                gamma = next;
                ll = llNext;
                if (gain <= 1e-13 * (Math.Abs(ll) + 1)) break;
            }
            return (ll, gamma);
        }

        /// <summary>Full log-likelihood Σ log f(s²_g), including every constant, so that values at different d0 compare.</summary>
        internal static double LogLik(double[] s2, double[] df, Matrix<double> basis, double d0, Vector<double> gamma)
        {
            var tau = basis * gamma;
            double sum = 0;
            bool inf = double.IsPositiveInfinity(d0);
            double lgD0 = inf ? 0 : SpecialFunctions.GammaLn(d0 / 2);
            for (int g = 0; g < s2.Length; g++)
            {
                double d = df[g], x = s2[g] * Math.Exp(-tau[g]);
                double common = -SpecialFunctions.GammaLn(d / 2) + (d / 2 - 1) * Math.Log(x) - tau[g];
                if (inf)
                    sum += common + d / 2 * Math.Log(d / 2) - d / 2 * x;
                else
                    sum += common + SpecialFunctions.GammaLn((d + d0) / 2) - lgD0 + d / 2 * Math.Log(d / d0)
                           - (d + d0) / 2 * Log1P(d * x / d0);
            }
            return sum;
        }

        /// <summary>Gradient and Hessian of the log-likelihood in γ at fixed d0.</summary>
        private static (Vector<double> grad, Matrix<double> hess) Derivatives(double[] s2, double[] df, Matrix<double> basis,
            double d0, Vector<double> gamma)
        {
            int k = basis.ColumnCount;
            var tau = basis * gamma;
            var grad = Vector<double>.Build.Dense(k);
            var hess = Matrix<double>.Build.Dense(k, k);
            bool inf = double.IsPositiveInfinity(d0);
            for (int g = 0; g < s2.Length; g++)
            {
                double d = df[g], x = s2[g] * Math.Exp(-tau[g]);
                double first, second;
                if (inf)
                {
                    // ℓ(τ) = −(d/2)τ − (d/2) s² e^{−τ} + const
                    first = -d / 2 + d / 2 * x;
                    second = -d / 2 * x;
                }
                else
                {
                    // ℓ(τ) = −(d/2)τ − ((d + d0)/2) log(1 + (d/d0) s² e^{−τ}) + const
                    double a = d * x / d0, q = a / (1 + a);
                    first = -d / 2 + (d + d0) / 2 * q;
                    second = -(d + d0) / 2 * q * (1 - q);
                }
                for (int i = 0; i < k; i++)
                {
                    grad[i] += first * basis[g, i];
                    for (int j = 0; j < k; j++) hess[i, j] += second * basis[g, i] * basis[g, j];
                }
            }
            return (grad, hess);
        }
        /// <summary>log(1 + y), accurate for small y (Goldberg's correction).</summary>
        private static double Log1P(double y)
        {
            double u = 1 + y;
            return u == 1 ? y : Math.Log(u) * y / (u - 1);
        }

        /// <summary>Brent's (1973) derivative-free minimizer on [a, b].</summary>
        private static double Brent(Func<double, double> f, double a, double b, double tol)
        {
            const double golden = 0.3819660112501051;
            double x = a + golden * (b - a), w = x, v = x, fx = f(x), fw = fx, fv = fx, d = 0, e = 0;
            for (int iter = 0; iter < 200; iter++)
            {
                double mid = 0.5 * (a + b), tol1 = tol * Math.Abs(x) + 1e-15, tol2 = 2 * tol1;
                if (Math.Abs(x - mid) <= tol2 - 0.5 * (b - a)) break;
                bool goldenStep = true;
                if (Math.Abs(e) > tol1)
                {
                    double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
                    q = 2 * (q - r);
                    if (q > 0) p = -p; else q = -q;
                    double eTemp = e;
                    e = d;
                    if (Math.Abs(p) < Math.Abs(0.5 * q * eTemp) && p > q * (a - x) && p < q * (b - x))
                    {
                        d = p / q;
                        double u0 = x + d;
                        if (u0 - a < tol2 || b - u0 < tol2) d = mid >= x ? tol1 : -tol1;
                        goldenStep = false;
                    }
                }
                if (goldenStep) { e = (x >= mid ? a : b) - x; d = golden * e; }
                double u = Math.Abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1);
                double fu = f(u);
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
                    else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
                }
            }
            return x;
        }
    }
}
