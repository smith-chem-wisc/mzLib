using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// A small feed-forward network (tanh hidden layers, sigmoid output, cross-entropy, Adam) and an averaged ensemble of
/// them. It gives the target-decoy rescorer a non-linear model, as DIA-NN's classifier does.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MultilayerPerceptronTests
{
    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    /// <summary>Positives inside a ring, negatives outside: no straight line separates them.</summary>
    private static (List<double[]> X, List<bool> Y) Ring(int n, int seed)
    {
        var random = new Random(seed);
        var x = new List<double[]>();
        var y = new List<bool>();
        for (int i = 0; i < n; i++)
        {
            double a = 3 * Gaussian(random), b = 3 * Gaussian(random);
            x.Add([a, b]);
            y.Add(Math.Sqrt(a * a + b * b) < 3);
        }
        return (x, y);
    }

    private static double Accuracy(Func<double[], double> predict, List<double[]> x, List<bool> y) =>
        x.Select((v, i) => (predict(v) > 0.5) == y[i]).Count(ok => ok) / (double)x.Count;

    [Test]
    public void TheNetworkLearnsWhatNoLineCan()
    {
        var (x, y) = Ring(3000, 1);
        var (tx, ty) = Ring(1000, 2);

        var net = MultilayerPerceptron.Train(x, y, hiddenLayers: [16, 8], epochs: 60, seed: 7);
        var line = LinearDiscriminant.Fit(x, y);

        Assert.That(Accuracy(net.Predict, tx, ty), Is.GreaterThan(0.93), "network on held-out data");
        Assert.That(Accuracy(v => line.Score(v) > 0 ? 1 : 0, tx, ty), Is.LessThan(0.75), "a line cannot do it");
    }

    /// <summary>
    /// Training and prediction to the bit, pinned on the DIA rescorer's architecture: the values were produced by the
    /// network before its weights moved from 2-D arrays to flat ones, which must change nothing but speed.
    /// </summary>
    [Test]
    public void TrainingIsPinnedToTheBit()
    {
        var (x, y) = Ring(400, 5);
        var random = new Random(9);
        var rows = x.Select(v => new[] { v[0], v[1], Gaussian(random), v[0] * v[1], Gaussian(random), 1.0 }).ToList();

        var net = MultilayerPerceptron.Train(rows, y, [25, 20, 15, 10, 5], epochs: 3, seed: 13);

        long[] bits = rows.Take(4).Select(r => BitConverter.DoubleToInt64Bits(net.PredictLogit(r))).ToArray();
        TestContext.Out.WriteLine("pinned: " + string.Join(", ", bits));
        Assert.That(bits, Is.EqualTo(PinnedLogitBits));
    }

    /// <summary>
    /// Scoring many rows at once (a row per vector lane) gives each row exactly its one-at-a-time logit, for any count of rows,
    /// including those left over when the count is not a multiple of the vector width.
    /// </summary>
    [Test]
    public void ManyRowsAtOnceScoreAsOneAtATime()
    {
        var (x, y) = Ring(300, 6);
        var random = new Random(4);
        var rows = x.Select(v => new[] { v[0], v[1], Gaussian(random), v[0] * v[1], Gaussian(random), 1.0 }).ToList();
        var ensemble = MultilayerPerceptron.TrainEnsemble(rows, y, members: 3, [25, 20, 15, 10, 5], epochs: 2, seed: 3);

        foreach (int count in new[] { 0, 1, 3, 4, 7, 8, 9, 64, 299 })
        {
            int[] which = Enumerable.Range(0, count).Select(k => (k * 37) % rows.Count).ToArray();
            var batch = new double[rows.Count];
            ensemble.PredictLogits(rows, which, batch);
            foreach (int i in which)
                Assert.That(BitConverter.DoubleToInt64Bits(batch[i]), Is.EqualTo(BitConverter.DoubleToInt64Bits(ensemble.PredictLogit(rows[i]))), $"{count} rows, row {i}");
        }
    }

    private static readonly long[] PinnedLogitBits = [-4618388307932396066, -4622205833387768178, -4618959286307770040, -4633051758969041453];

    [Test]
    public void TrainingIsDeterministicForASeed()
    {
        var (x, y) = Ring(500, 3);

        var a = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 11);
        var b = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 11);
        var c = MultilayerPerceptron.Train(x, y, [8, 4], epochs: 5, seed: 12);

        Assert.That(x.Select(a.Predict), Is.EqualTo(x.Select(b.Predict)));
        Assert.That(x.Select(a.Predict), Is.Not.EqualTo(x.Select(c.Predict)));
    }

    /// <summary>Outputs are probabilities, and the ensemble is the mean of its members.</summary>
    [Test]
    public void TheEnsembleAveragesItsMembers()
    {
        var (x, y) = Ring(500, 4);

        var ensemble = MultilayerPerceptron.TrainEnsemble(x, y, members: 3, [8, 4], epochs: 5, seed: 5);

        var v = x[17];
        Assert.That(ensemble.Members, Has.Count.EqualTo(3));
        Assert.That(ensemble.Predict(v), Is.EqualTo(ensemble.Members.Average(m => m.Predict(v))).Within(1e-12));
        Assert.That(x.Select(ensemble.Predict), Is.All.InRange(0.0, 1.0));
    }

    /// <summary>Members train in parallel, and the ensemble is still the same for the same seed.</summary>
    [Test]
    public void TheEnsembleIsDeterministicThoughItsMembersTrainInParallel()
    {
        var (x, y) = Ring(400, 6);

        var a = MultilayerPerceptron.TrainEnsemble(x, y, members: 4, [8, 4], epochs: 3, seed: 9);
        var b = MultilayerPerceptron.TrainEnsemble(x, y, members: 4, [8, 4], epochs: 3, seed: 9);

        Assert.That(x.Select(a.Predict), Is.EqualTo(x.Select(b.Predict)));
    }

    /// <summary>
    /// The ensemble's log-odds, computed from the members' raw logits. Where probabilities are moderate it is the logit of
    /// the mean probability.
    /// </summary>
    [Test]
    public void TheEnsembleLogitIsTheLogOddsOfItsMeanProbability()
    {
        var (x, y) = Ring(500, 8);
        var ensemble = MultilayerPerceptron.TrainEnsemble(x, y, members: 3, [8, 4], epochs: 5, seed: 3);

        foreach (var v in x.Take(50))
        {
            double p = ensemble.Predict(v);
            Assert.That(ensemble.PredictLogit(v), Is.EqualTo(Math.Log(p / (1 - p))).Within(1e-9));
        }
    }

    /// <summary>
    /// A confident member's probability rounds to exactly 1 in double precision (from a logit of about 37), so every
    /// confident row would get the same score and its order would be lost. On the whole-proteome DIA search, saturated
    /// targets and decoys tied and the folds' normalisation decided between them. Computed from the logits, confident
    /// rows stay distinct and ordered.
    /// </summary>
    [Test]
    public void ConfidentPredictionsStayDistinctAndOrdered()
    {
        double a = MultilayerPerceptronEnsemble.LogitOfMeanProbability([40, 45]);
        double b = MultilayerPerceptronEnsemble.LogitOfMeanProbability([41, 45]);
        double c = MultilayerPerceptronEnsemble.LogitOfMeanProbability([60, 70]);

        Assert.That(1 / (1 + Math.Exp(-40.0)), Is.EqualTo(1.0), "the naive probability has saturated");
        Assert.That(new[] { a, b, c }, Is.All.Matches<double>(double.IsFinite));
        Assert.That(a, Is.LessThan(b));
        Assert.That(b, Is.LessThan(c));
        // log-odds of the mean probability: 1 - p averages exp(-40) and exp(-45)
        Assert.That(a, Is.EqualTo(-Math.Log((Math.Exp(-40) + Math.Exp(-45)) / 2)).Within(1e-9));
        Assert.That(MultilayerPerceptronEnsemble.LogitOfMeanProbability([-60, -70]), Is.EqualTo(-c).Within(1e-9), "symmetric for negatives");
    }

    [Test]
    public void ArgumentsAreChecked()
    {
        var (x, y) = Ring(50, 5);
        Assert.Throws<ArgumentNullException>(() => MultilayerPerceptron.Train(null!, y, [4], 1, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, y.Take(10).ToList(), [4], 1, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, x.Select(_ => true).ToList(), [4], 1, 1), "both classes needed");
        Assert.Throws<ArgumentOutOfRangeException>(() => MultilayerPerceptron.Train(x, y, [0], 1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultilayerPerceptron.Train(x, y, [4], 0, 1));
        Assert.Throws<ArgumentException>(() => MultilayerPerceptron.Train(x, y, [4], 1, 1).Predict([1.0]), "feature count");
    }

    /// <summary>
    /// Prediction is the hot loop of a DIA rescoring (12 networks x millions of rows x 3 folds), so it must not allocate per
    /// row; these values, from the allocating version it replaced, pin that it computes exactly the same thing.
    /// </summary>
    [Test]
    public void PredictionsArePinnedExactly()
    {
        var r = new Random(3);
        var x = Enumerable.Range(0, 400).Select(i => Enumerable.Range(0, 7).Select(_ => r.NextDouble() * 4 - 2 + (i % 2 == 0 ? 0.7 : 0)).ToArray()).ToList();
        var y = Enumerable.Range(0, 400).Select(i => i % 2 == 0).ToList();
        var ensemble = MultilayerPerceptron.TrainEnsemble(x, y, 3, new[] { 25, 20, 15, 10, 5 }, 2, seed: 11);
        double[] ensembleLogits = [0.3470434694645611, -0.7743676515512511, 0.9789297911956628, -0.9852972843394292, -0.6267818087740666];
        double[] memberLogits = [0.7590940525529443, -0.3878021032431866, 0.833344523382401, -0.9729353089861892, -0.5548285858996359];
        for (int i = 0; i < 5; i++)
        {
            Assert.That(ensemble.PredictLogit(x[i]), Is.EqualTo(ensembleLogits[i]), $"ensemble, row {i}");
            Assert.That(ensemble.Members[0].PredictLogit(x[i]), Is.EqualTo(memberLogits[i]), $"member, row {i}");
        }
    }
}
