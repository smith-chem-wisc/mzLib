using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace StatisticalModels
{
    /// <summary>
    /// A fitted linear discriminant: <c>score = Bias + Σ Weights[i] · x[i]</c>, in the original feature units.
    /// Higher scores look more like the positive class.
    /// </summary>
    public sealed class LinearDiscriminantFit
    {
        internal LinearDiscriminantFit(double[] weights, double bias)
        {
            _weights = weights;
            Bias = bias;
        }

        private readonly double[] _weights;

        /// <summary>One weight per feature; exactly 0 for a feature that was constant in training.</summary>
        public IReadOnlyList<double> Weights => _weights;

        public double Bias { get; }

        public double Score(ReadOnlySpan<double> features)
        {
            if (features.Length != _weights.Length)
                throw new ArgumentException($"Expected {_weights.Length} features; got {features.Length}.", nameof(features));
            double score = Bias;
            for (int i = 0; i < _weights.Length; i++)
                score += _weights[i] * features[i];
            return score;
        }
    }

    /// <summary>
    /// Fisher's linear discriminant with a ridge on the pooled within-class covariance, fitted on standardized features.
    /// It is closed-form and deterministic, and it stays finite when the classes are separable, which semi-supervised
    /// target-decoy training makes them almost by construction. An unpenalized logistic regression would report
    /// separation there instead of a direction.
    /// </summary>
    public static class LinearDiscriminant
    {
        /// <summary>Ridge added to the pooled covariance of the standardized features.</summary>
        public const double DefaultRidge = 1e-3;

        /// <summary>A feature whose standard deviation is at or below this is treated as constant.</summary>
        internal const double ConstantTolerance = 1e-12;

        /// <param name="features">One row per observation; every row the same length.</param>
        /// <param name="isPositive">One label per row; both classes must be present.</param>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">
        /// Lengths disagree, rows are ragged, a value is not finite, or a class is missing.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="ridge"/> is negative or not finite.</exception>
        public static LinearDiscriminantFit Fit(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, double ridge = DefaultRidge)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isPositive);
            if (!(ridge >= 0) || !double.IsFinite(ridge))
                throw new ArgumentOutOfRangeException(nameof(ridge), ridge, "The ridge must be a finite, non-negative number.");
            if (features.Count != isPositive.Count)
                throw new ArgumentException($"There are {features.Count} feature rows but {isPositive.Count} labels.", nameof(isPositive));
            if (features.Count == 0)
                throw new ArgumentException("There are no observations.", nameof(features));
            int p = features[0].Length;
            foreach (var row in features)
            {
                if (row is null || row.Length != p)
                    throw new ArgumentException("Every feature row must have the same length.", nameof(features));
                if (!row.All(double.IsFinite))
                    throw new ArgumentException("Every feature value must be finite.", nameof(features));
            }
            int positives = isPositive.Count(b => b);
            if (positives == 0 || positives == isPositive.Count)
                throw new ArgumentException("Both classes must be present.", nameof(isPositive));

            // Standardize; constant features take no part. Row by row (each row read once, not once per feature); every
            // feature's sums still add the rows in the same order, so the result is unchanged
            int n = features.Count;
            var mean = new double[p];
            var sd = new double[p];
            for (int i = 0; i < n; i++)
            {
                var row = features[i];
                for (int j = 0; j < p; j++) mean[j] += row[j];
            }
            for (int j = 0; j < p; j++) mean[j] /= n;
            var squares = new double[p];
            for (int i = 0; i < n; i++)
            {
                var row = features[i];
                for (int j = 0; j < p; j++) squares[j] += (row[j] - mean[j]) * (row[j] - mean[j]);
            }
            for (int j = 0; j < p; j++)
                sd[j] = Math.Sqrt(squares[j] / Math.Max(1, n - 1));
            int[] active = Enumerable.Range(0, p).Where(j => sd[j] > ConstantTolerance).ToArray();
            var weights = new double[p];
            if (active.Length == 0)
                return new LinearDiscriminantFit(weights, 0);

            int k = active.Length;
            var classMean = new double[2, k];
            var count = new int[2];
            for (int i = 0; i < n; i++)
            {
                int c = isPositive[i] ? 1 : 0;
                count[c]++;
                for (int a = 0; a < k; a++)
                    classMean[c, a] += (features[i][active[a]] - mean[active[a]]) / sd[active[a]];
            }
            for (int c = 0; c < 2; c++)
                for (int a = 0; a < k; a++)
                    classMean[c, a] /= count[c];

            // The pooled scatter in a flat array, upper triangle only, then mirrored: each cell adds the same products in the same
            // row order (centered[a] * centered[b] is exactly centered[b] * centered[a]), so the matrix is unchanged; the matrix
            // indexer over every cell was most of the fit's time
            var scatter = new double[k * k];
            var centered = new double[k];
            for (int i = 0; i < n; i++)
            {
                int c = isPositive[i] ? 1 : 0;
                var row = features[i];
                for (int a = 0; a < k; a++)
                    centered[a] = (row[active[a]] - mean[active[a]]) / sd[active[a]] - classMean[c, a];
                for (int a = 0; a < k; a++)
                {
                    double x = centered[a];
                    int at = a * k;
                    for (int b = a; b < k; b++)
                        scatter[at + b] += x * centered[b];
                }
            }
            var within = Matrix<double>.Build.Dense(k, k, (a, b) => a <= b ? scatter[a * k + b] : scatter[b * k + a]);
            within = within.Divide(Math.Max(1, n - 2));
            for (int a = 0; a < k; a++)
                within[a, a] += ridge;

            var difference = Vector<double>.Build.Dense(k, a => classMean[1, a] - classMean[0, a]);
            Vector<double> standardizedWeights = within.Cholesky().Solve(difference);

            // Back to original units; the bias puts 0 midway between the class means
            double bias = 0;
            for (int a = 0; a < k; a++)
            {
                int j = active[a];
                weights[j] = standardizedWeights[a] / sd[j];
                bias -= weights[j] * mean[j] + standardizedWeights[a] * (classMean[0, a] + classMean[1, a]) / 2;
            }
            return new LinearDiscriminantFit(weights, bias);
        }
    }
}
