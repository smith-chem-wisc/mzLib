using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace StatisticalModels
{
    /// <summary>
    /// A small feed-forward binary classifier: standardised inputs, tanh hidden layers and a sigmoid output, trained on
    /// cross-entropy with Adam in mini-batches. Training is deterministic for a seed. It gives the target-decoy rescorer a
    /// non-linear model, as DIA-NN's classifier ensemble does (Demichev et al. 2020, Nat. Methods 17:41).
    /// </summary>
    public sealed class MultilayerPerceptron
    {
        private readonly double[] _mean;
        private readonly double[] _sd;
        // Layer l: [inputs, outputs] row-major in one flat array (element (i, o) at i * outputs + o); a 2-D array's
        // two-index bounds checks were most of training's time
        private readonly double[][] _weights;
        private readonly int[] _inputs;
        private readonly int[] _outputs;
        private readonly double[][] _biases;

        private MultilayerPerceptron(double[] mean, double[] sd, double[][] weights, int[] inputs, int[] outputs, double[][] biases)
        {
            _mean = mean;
            _sd = sd;
            _weights = weights;
            _inputs = inputs;
            _outputs = outputs;
            _biases = biases;
        }

        /// <summary>The probability that <paramref name="features"/> belong to the positive class.</summary>
        /// <exception cref="ArgumentException">The feature count differs from training.</exception>
        public double Predict(double[] features) => Sigmoid(PredictLogit(features));

        /// <summary>The log-odds that <paramref name="features"/> belong to the positive class: the output before the sigmoid,
        /// which keeps its resolution where the probability has rounded to 0 or 1.</summary>
        /// <exception cref="ArgumentException">The feature count differs from training.</exception>
        public double PredictLogit(double[] features)
        {
            ArgumentNullException.ThrowIfNull(features);
            if (features.Length != _mean.Length)
                throw new ArgumentException($"Expected {_mean.Length} features; got {features.Length}.", nameof(features));
            Span<double> first = stackalloc double[Width];
            Span<double> second = stackalloc double[Width];
            return Forward(features, first, second);
        }

        /// <summary>Features the network was trained on.</summary>
        internal int InputCount => _mean.Length;

        /// <summary>The widest layer, inputs included: the size of each buffer <see cref="Forward"/> needs.</summary>
        internal int Width => Math.Max(_mean.Length, _outputs.Max());

        /// <summary>
        /// <see cref="PredictLogit(double[])"/> on two caller buffers of at least <see cref="Width"/>, allocating nothing:
        /// prediction is the hot loop of a DIA rescoring. Same arithmetic, in the same order.
        /// </summary>
        internal double Forward(double[] features, Span<double> first, Span<double> second)
        {
            int width = _mean.Length;
            for (int j = 0; j < width; j++)
                first[j] = (features[j] - _mean[j]) / _sd[j];
            bool inFirst = true;
            for (int l = 0; l < _weights.Length; l++)
            {
                var w = _weights[l];
                bool last = l == _weights.Length - 1;
                int outputs = _outputs[l];
                Span<double> input = (inFirst ? first : second)[..width];
                Span<double> output = (inFirst ? second : first)[..outputs];
                // Each output from its bias over the inputs in increasing order, as a dot product per output would sum it,
                // but vectorised across the outputs (AddScaled)
                _biases[l].AsSpan().CopyTo(output);
                for (int i = 0; i < width; i++)
                    AddScaled(output, w.AsSpan(i * outputs, outputs), input[i]);
                if (!last)
                    for (int o = 0; o < outputs; o++)
                        output[o] = Math.Tanh(output[o]);
                width = outputs;
                inFirst = !inFirst;
            }
            return (inFirst ? first : second)[0];
        }

        /// <summary>
        /// <see cref="Forward"/> for <see cref="Vector{T}.Count"/> rows at once, one per lane: each lane gets exactly its row's
        /// operations in the same order (standardise; per output, the bias then one product and one addition per input in
        /// increasing order; tanh per lane), so each lane's logit is bit-identical to <see cref="Forward"/>. Narrow layers
        /// (down to the single output) use whole vectors this way.
        /// </summary>
        internal Vector<double> ForwardLanes(IReadOnlyList<double[]> features, ReadOnlySpan<int> rows, Vector<double>[] first, Vector<double>[] second)
        {
            Span<double> lane = stackalloc double[Vector<double>.Count];
            int width = _mean.Length;
            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < lane.Length; k++)
                    lane[k] = (features[rows[k]][j] - _mean[j]) / _sd[j];
                first[j] = new Vector<double>(lane);
            }
            bool inFirst = true;
            for (int l = 0; l < _weights.Length; l++)
            {
                var w = _weights[l];
                var bias = _biases[l];
                bool last = l == _weights.Length - 1;
                int outputs = _outputs[l];
                var input = inFirst ? first : second;
                var output = inFirst ? second : first;
                for (int o = 0; o < outputs; o++)
                    output[o] = new Vector<double>(bias[o]);
                for (int i = 0; i < width; i++)
                {
                    var x = input[i];
                    int at = i * outputs;
                    for (int o = 0; o < outputs; o++)
                        output[o] += x * new Vector<double>(w[at + o]);
                }
                if (!last)
                    for (int o = 0; o < outputs; o++)
                    {
                        output[o].CopyTo(lane);
                        for (int k = 0; k < lane.Length; k++)
                            lane[k] = Math.Tanh(lane[k]);
                        output[o] = new Vector<double>(lane);
                    }
                width = outputs;
                inFirst = !inFirst;
            }
            return (inFirst ? first : second)[0];
        }

        /// <param name="hiddenLayers">Units in each hidden layer, input side first.</param>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">Lengths disagree, rows are ragged or not finite, or a class is missing.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A layer has no units, or a count or rate is not positive.</exception>
        public static MultilayerPerceptron Train(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, IReadOnlyList<int> hiddenLayers,
            int epochs, int seed, int batchSize = 50, double learningRate = 0.003)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(isPositive);
            ArgumentNullException.ThrowIfNull(hiddenLayers);
            if (features.Count != isPositive.Count)
                throw new ArgumentException($"There are {features.Count} rows but {isPositive.Count} labels.", nameof(isPositive));
            if (features.Count == 0 || isPositive.All(b => b) || !isPositive.Any(b => b))
                throw new ArgumentException("Both classes must be present.", nameof(isPositive));
            if (hiddenLayers.Any(units => units < 1))
                throw new ArgumentOutOfRangeException(nameof(hiddenLayers), "Every hidden layer needs at least one unit.");
            ArgumentOutOfRangeException.ThrowIfLessThan(epochs, 1);
            ArgumentOutOfRangeException.ThrowIfLessThan(batchSize, 1);
            if (!(learningRate > 0))
                throw new ArgumentOutOfRangeException(nameof(learningRate), learningRate, "The learning rate must be positive.");
            int p = features[0].Length;
            if (features.Any(row => row is null || row.Length != p || !row.All(double.IsFinite)))
                throw new ArgumentException("Every row needs the same number of finite features.", nameof(features));

            // Standardisation from the training rows; a constant feature is centred only
            var mean = new double[p];
            var sd = new double[p];
            for (int j = 0; j < p; j++)
            {
                mean[j] = features.Average(row => row[j]);
                double variance = features.Sum(row => (row[j] - mean[j]) * (row[j] - mean[j])) / Math.Max(1, features.Count - 1);
                sd[j] = variance > 1e-24 ? Math.Sqrt(variance) : 1;
            }

            var random = new Random(seed);
            int[] sizes = [p, .. hiddenLayers, 1];
            var weights = new double[sizes.Length - 1][];
            var biases = new double[sizes.Length - 1][];
            for (int l = 0; l < weights.Length; l++)
            {
                double limit = Math.Sqrt(6.0 / (sizes[l] + sizes[l + 1])); // Xavier/Glorot uniform
                weights[l] = new double[sizes[l] * sizes[l + 1]];
                for (int i = 0; i < sizes[l]; i++)
                    for (int o = 0; o < sizes[l + 1]; o++)
                        weights[l][i * sizes[l + 1] + o] = (2 * random.NextDouble() - 1) * limit;
                biases[l] = new double[sizes[l + 1]];
            }
            var net = new MultilayerPerceptron(mean, sd, weights, sizes[..^1], sizes[1..], biases);
            net.Fit(features, isPositive, epochs, batchSize, learningRate, random);            return net;
        }

        /// <summary>An ensemble of <paramref name="members"/> networks, each seeded differently; it predicts their mean.</summary>
        public static MultilayerPerceptronEnsemble TrainEnsemble(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, int members,
            IReadOnlyList<int> hiddenLayers, int epochs, int seed, int batchSize = 50, double learningRate = 0.003)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(members, 1);
            // Members train in parallel; each has its own seed, so the result is the same as training them in turn
            var trained = new MultilayerPerceptron[members];
            System.Threading.Tasks.Parallel.For(0, members, m => trained[m] = Train(features, isPositive, hiddenLayers, epochs, seed + 7919 * m, batchSize, learningRate));
            return new MultilayerPerceptronEnsemble(trained);
        }

        private void Fit(IReadOnlyList<double[]> features, IReadOnlyList<bool> isPositive, int epochs, int batchSize, double learningRate, Random random)
        {
            int layers = _weights.Length;
            var mW = _weights.Select(w => new double[w.Length]).ToArray();
            var vW = _weights.Select(w => new double[w.Length]).ToArray();
            var mB = _biases.Select(b => new double[b.Length]).ToArray();
            var vB = _biases.Select(b => new double[b.Length]).ToArray();
            var gW = _weights.Select(w => new double[w.Length]).ToArray();
            var gB = _biases.Select(b => new double[b.Length]).ToArray();
            const double beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8;
            int step = 0;
            int[] order = Enumerable.Range(0, features.Count).ToArray();
            // Buffers reused for every row: activations per layer, and the backpropagated error at each layer's output
            var activations = new double[layers + 1][];
            activations[0] = new double[_mean.Length];
            for (int l = 0; l < layers; l++)
                activations[l + 1] = new double[_outputs[l]];
            var deltas = activations.Select(a => new double[a.Length]).ToArray();

            for (int epoch = 0; epoch < epochs; epoch++)
            {
                for (int i = order.Length - 1; i > 0; i--) // Fisher–Yates with the seeded generator
                {
                    int k = random.Next(i + 1);
                    (order[i], order[k]) = (order[k], order[i]);
                }
                for (int start = 0; start < order.Length; start += batchSize)
                {
                    int end = Math.Min(order.Length, start + batchSize);
                    foreach (var g in gW) Array.Clear(g);
                    foreach (var g in gB) Array.Clear(g);
                    for (int r = start; r < end; r++)
                    {
                        int row = order[r];
                        StandardiseInto(features[row], activations[0]);
                        for (int l = 0; l < layers; l++)
                            LayerInto(activations[l], l, last: l == layers - 1, activations[l + 1]);
                        // Cross-entropy on the sigmoid output: dLoss/dz = p - y
                        double[] delta = deltas[layers];
                        delta[0] = Sigmoid(activations[layers][0]) - (isPositive[row] ? 1 : 0);
                        for (int l = layers - 1; l >= 0; l--)
                        {
                            double[] input = activations[l];
                            double[] g = gW[l], w = _weights[l];
                            int outputs = delta.Length;
                            for (int o = 0; o < delta.Length; o++)
                                gB[l][o] += delta[o];
                            // Row by row (contiguous); every gradient cell gets the same single addition as before
                            for (int i = 0; i < input.Length; i++)
                            {
                                AddScaled(g.AsSpan(i * outputs, outputs), delta, input[i]);
                            }
                            if (l == 0)
                                break;
                            double[] previous = deltas[l];
                            for (int i = 0; i < input.Length; i++)
                            {
                                double sum = 0;
                                int at = i * outputs;
                                for (int o = 0; o < outputs; o++)
                                    sum += w[at + o] * delta[o];
                                previous[i] = sum * (1 - input[i] * input[i]); // tanh'
                            }
                            delta = previous;
                        }
                    }

                    step++;
                    double n = end - start;
                    double correction1 = 1 - Math.Pow(beta1, step), correction2 = 1 - Math.Pow(beta2, step);
                    for (int l = 0; l < layers; l++)
                    {
                        // Cell by cell in the same (i, o) order as before; the flat index is i * outputs + o
                        double[] wl = _weights[l], gl = gW[l], ml = mW[l], vl = vW[l];
                        int k = 0;
                        if (Vector.IsHardwareAccelerated)
                        {
                            // The same operations per cell, in the same order, Vector<double>.Count cells at a time
                            var vn = new Vector<double>(n);
                            var vBeta1 = new Vector<double>(beta1);
                            var vOneMinusBeta1 = new Vector<double>(1 - beta1);
                            var vBeta2 = new Vector<double>(beta2);
                            var vOneMinusBeta2 = new Vector<double>(1 - beta2);
                            var vRate = new Vector<double>(learningRate);
                            var vCorrection1 = new Vector<double>(correction1);
                            var vCorrection2 = new Vector<double>(correction2);
                            var vEpsilon = new Vector<double>(epsilon);
                            for (; k <= wl.Length - Vector<double>.Count; k += Vector<double>.Count)
                            {
                                var g = new Vector<double>(gl, k) / vn;
                                var m = vBeta1 * new Vector<double>(ml, k) + vOneMinusBeta1 * g;
                                var v = vBeta2 * new Vector<double>(vl, k) + vOneMinusBeta2 * g * g;
                                m.CopyTo(ml, k);
                                v.CopyTo(vl, k);
                                (new Vector<double>(wl, k) - vRate * (m / vCorrection1) / (Vector.SquareRoot(v / vCorrection2) + vEpsilon)).CopyTo(wl, k);
                            }
                        }
                        for (; k < wl.Length; k++)
                        {
                            double g = gl[k] / n;
                            ml[k] = beta1 * ml[k] + (1 - beta1) * g;
                            vl[k] = beta2 * vl[k] + (1 - beta2) * g * g;
                            wl[k] -= learningRate * (ml[k] / correction1) / (Math.Sqrt(vl[k] / correction2) + epsilon);
                        }
                        for (int o = 0; o < _biases[l].Length; o++)
                        {
                            double g = gB[l][o] / n;
                            mB[l][o] = beta1 * mB[l][o] + (1 - beta1) * g;
                            vB[l][o] = beta2 * vB[l][o] + (1 - beta2) * g * g;
                            _biases[l][o] -= learningRate * (mB[l][o] / correction1) / (Math.Sqrt(vB[l][o] / correction2) + epsilon);
                        }
                    }
                }
            }
        }

        private void StandardiseInto(double[] features, double[] x)
        {
            for (int j = 0; j < x.Length; j++)
                x[j] = (features[j] - _mean[j]) / _sd[j];
        }

        private void LayerInto(double[] input, int l, bool last, double[] output)
        {
            // Row by row through the weights (contiguous), accumulating each output from its bias over the inputs in
            // increasing order: the same sums, in the same order, as one output at a time
            var w = _weights[l];
            var bias = _biases[l];
            for (int o = 0; o < output.Length; o++)
                output[o] = bias[o];
            int outputs = output.Length;
            for (int i = 0; i < input.Length; i++)
                AddScaled(output, w.AsSpan(i * outputs, outputs), input[i]);
            if (!last)
                for (int o = 0; o < output.Length; o++)
                    output[o] = Math.Tanh(output[o]);
        }

        /// <summary>
        /// destination[o] += x * source[o], vectorised: each cell gets the same one product and one addition as the scalar
        /// loop (no fused multiply-add), so the result is bit-identical.
        /// </summary>
        private static void AddScaled(Span<double> destination, ReadOnlySpan<double> source, double x)
        {
            int o = 0;
            if (Vector.IsHardwareAccelerated)
            {
                var vx = new Vector<double>(x);
                for (; o <= destination.Length - Vector<double>.Count; o += Vector<double>.Count)
                    (new Vector<double>(destination[o..]) + new Vector<double>(source[o..]) * vx).CopyTo(destination[o..]);
            }
            for (; o < destination.Length; o++)
                destination[o] += x * source[o];
        }

        private static double Sigmoid(double z) => 1 / (1 + Math.Exp(-z));
    }

    /// <summary>An averaged ensemble of <see cref="MultilayerPerceptron"/>s.</summary>
    public sealed class MultilayerPerceptronEnsemble
    {
        internal MultilayerPerceptronEnsemble(IReadOnlyList<MultilayerPerceptron> members) => Members = members;

        public IReadOnlyList<MultilayerPerceptron> Members { get; }

        /// <summary>The members' mean probability.</summary>
        public double Predict(double[] features) => Members.Average(m => m.Predict(features));

        /// <summary>
        /// The log-odds of the members' mean probability, computed from their logits so that it never saturates: a confident
        /// member's probability rounds to exactly 1 in double precision, and every confident row would tie.
        /// </summary>
        public double PredictLogit(double[] features)
        {
            ArgumentNullException.ThrowIfNull(features);
            if (Members.Count == 0 || features.Length != Members[0].InputCount)
                return LogitOfMeanProbability(Members.Select(m => m.PredictLogit(features)).ToArray()); // the members report the error
            int width = Members.Max(m => m.Width);
            Span<double> first = stackalloc double[width];
            Span<double> second = stackalloc double[width];
            Span<double> logits = stackalloc double[Members.Count];
            for (int m = 0; m < Members.Count; m++)
                logits[m] = Members[m].Forward(features, first, second);
            return LogitOfMeanProbability(logits);
        }

        /// <summary>
        /// <see cref="PredictLogit"/> for each of <paramref name="rows"/> (indices into <paramref name="features"/>), written to
        /// <paramref name="into"/> at the same index. Bit-identical to calling it row by row, and several times cheaper: rows
        /// are scored a vector's width at a time (<see cref="MultilayerPerceptron.ForwardLanes"/>), the rest one at a time.
        /// </summary>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">A row's feature count differs from training.</exception>
        public void PredictLogits(IReadOnlyList<double[]> features, IReadOnlyList<int> rows, double[] into)
        {
            ArgumentNullException.ThrowIfNull(features);
            ArgumentNullException.ThrowIfNull(rows);
            ArgumentNullException.ThrowIfNull(into);
            int lanes = Vector<double>.Count, members = Members.Count, k = 0;
            if (Vector.IsHardwareAccelerated && members > 0 && rows.Count >= lanes)
            {
                int width = Members.Max(m => m.Width), inputs = Members[0].InputCount;
                var first = new Vector<double>[width];
                var second = new Vector<double>[width];
                Span<int> block = stackalloc int[lanes];
                Span<double> logits = stackalloc double[lanes * members];
                for (; k <= rows.Count - lanes; k += lanes)
                {
                    bool valid = true;
                    for (int lane = 0; lane < lanes; lane++)
                    {
                        block[lane] = rows[k + lane];
                        valid &= features[block[lane]] is { } row && row.Length == inputs;
                    }
                    if (!valid)
                        break; // the row-by-row path below reports the bad row
                    for (int m = 0; m < members; m++)
                    {
                        var z = Members[m].ForwardLanes(features, block, first, second);
                        for (int lane = 0; lane < lanes; lane++)
                            logits[lane * members + m] = z[lane];
                    }
                    for (int lane = 0; lane < lanes; lane++)
                        into[block[lane]] = LogitOfMeanProbability(logits.Slice(lane * members, members));
                }
            }
            for (; k < rows.Count; k++)
                into[rows[k]] = PredictLogit(features[rows[k]]);
        }

        /// <summary>log(mean sigmoid(z)) - log(mean sigmoid(-z)), each a log-sum-exp of log-sigmoids.</summary>
        internal static double LogitOfMeanProbability(IReadOnlyList<double> logits) =>
            LogSumExp(logits.Select(z => -Softplus(-z))) - LogSumExp(logits.Select(z => -Softplus(z)));

        /// <summary><see cref="LogitOfMeanProbability(IReadOnlyList{double})"/> without allocating, summed in the same order.</summary>
        internal static double LogitOfMeanProbability(ReadOnlySpan<double> logits) => LogSumExpOfLogSigmoid(logits, -1) - LogSumExpOfLogSigmoid(logits, 1);

        /// <summary>log sum exp(-softplus(sign * z)) over the logits.</summary>
        private static double LogSumExpOfLogSigmoid(ReadOnlySpan<double> logits, double sign)
        {
            double max = double.NegativeInfinity;
            foreach (double z in logits)
                max = Math.Max(max, -Softplus(sign * z));
            double sum = 0;
            foreach (double z in logits)
                sum += Math.Exp(-Softplus(sign * z) - max);
            return max + Math.Log(sum);
        }

        private static double Softplus(double t) => Math.Max(t, 0) + Math.Log(1 + Math.Exp(-Math.Abs(t)));

        private static double LogSumExp(IEnumerable<double> values)
        {
            double[] v = values.ToArray();
            double max = v.Max();
            return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
        }
    }
}
