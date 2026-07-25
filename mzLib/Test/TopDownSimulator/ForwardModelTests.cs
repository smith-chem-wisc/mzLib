using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class ForwardModelTests
{
    private static ProteoformModel MakeModel(double mass, double rtCenter, double muZ, double abundance)
    {
        return new ProteoformModel(
            MonoisotopicMass: mass,
            Abundance: abundance,
            RtProfile: new EmgProfile(Mu: rtCenter, Sigma: 0.2, Tau: 0.1),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: muZ, SigmaZ: 1.5));
    }

    /// <summary>
    /// The original scan-outer formulation of Rasterize, kept here as the reference the optimized
    /// (species, charge)-outer implementation must reproduce exactly.
    /// </summary>
    private static double[,] RasterizeScanOuter(
        ProteoformModel[] models, int minCharge, int maxCharge, double sigmaMz,
        double[] scanTimes, double[] mzGrid)
    {
        const double evaluationWindowInSigmas = 6.0;
        const double minimumContributionWeight = 1e-18;

        var result = new double[scanTimes.Length, mzGrid.Length];
        var kernels = models.Select(m => new IsotopeEnvelopeKernel(m.MonoisotopicMass)).ToArray();
        double mzPadding = evaluationWindowInSigmas * sigmaMz;

        for (int s = 0; s < scanTimes.Length; s++)
        {
            for (int p = 0; p < models.Length; p++)
            {
                var model = models[p];
                var kernel = kernels[p];

                double rt = model.RtProfile.Evaluate(scanTimes[s]);
                if (rt <= 0) continue;

                double rtScaledAbundance = model.Abundance * rt;
                if (rtScaledAbundance <= minimumContributionWeight) continue;

                for (int z = minCharge; z <= maxCharge; z++)
                {
                    double weight = rtScaledAbundance * model.ChargeDistribution.Evaluate(z);
                    if (weight <= minimumContributionWeight) continue;

                    var (minMz, maxMz) = kernel.GetMzBounds(z);
                    for (int b = 0; b < mzGrid.Length; b++)
                    {
                        if (mzGrid[b] < minMz - mzPadding || mzGrid[b] > maxMz + mzPadding) continue;
                        result[s, b] += weight * kernel.Evaluate(mzGrid[b], z, sigmaMz);
                    }
                }
            }
        }

        return result;
    }

    [Test]
    public void RasterizeMatchesTheScanOuterReferenceBitForBit()
    {
        var models = new[]
        {
            MakeModel(mass: 10000.0, rtCenter: 20.0, muZ: 8, abundance: 1.5e6),
            MakeModel(mass: 10500.0, rtCenter: 20.1, muZ: 9, abundance: 4.2e5),
            MakeModel(mass: 12000.0, rtCenter: 21.5, muZ: 10, abundance: 9.1e5),
        };

        const int minCharge = 6;
        const int maxCharge = 12;
        const double sigmaMz = 0.012;

        double[] scanTimes = Enumerable.Range(0, 25).Select(i => 19.5 + i * 0.1).ToArray();
        double[] mzGrid = new GridRasterizer()
            .BuildCentroidMzAxis(models, minCharge, maxCharge, sigmaMz);

        var actual = new ForwardModel(models, minCharge, maxCharge, sigmaMz).Rasterize(scanTimes, mzGrid);
        var expected = RasterizeScanOuter(models, minCharge, maxCharge, sigmaMz, scanTimes, mzGrid);

        Assert.That(actual.GetLength(0), Is.EqualTo(expected.GetLength(0)));
        Assert.That(actual.GetLength(1), Is.EqualTo(expected.GetLength(1)));

        for (int s = 0; s < expected.GetLength(0); s++)
        {
            for (int b = 0; b < expected.GetLength(1); b++)
            {
                Assert.That(actual[s, b], Is.EqualTo(expected[s, b]),
                    $"Cell [{s},{b}] diverged from the reference implementation.");
            }
        }
    }

    [Test]
    public void IntensityPeaksAtRtCenterAndCentroid()
    {
        var model = MakeModel(mass: 10000.0, rtCenter: 20.0, muZ: 8, abundance: 1e6);
        var fm = new ForwardModel(new[] { model }, minCharge: 5, maxCharge: 12, sigmaMz: 0.01);

        var kernel = new IsotopeEnvelopeKernel(10000.0);
        double apexMz = kernel.CentroidMzs(8)[0];

        double onPeak = fm.Evaluate(t: 20.0, mz: apexMz);
        double offRt = fm.Evaluate(t: 25.0, mz: apexMz);
        double offMz = fm.Evaluate(t: 20.0, mz: apexMz + 1.0);

        Assert.That(onPeak, Is.GreaterThan(offRt * 10));
        Assert.That(onPeak, Is.GreaterThan(offMz * 10));
        Assert.That(onPeak, Is.GreaterThan(0));
    }

    [Test]
    public void RasterizeMatchesEvaluate()
    {
        var model = MakeModel(mass: 12000.0, rtCenter: 15.0, muZ: 10, abundance: 1.0);
        var fm = new ForwardModel(new[] { model }, minCharge: 7, maxCharge: 13, sigmaMz: 0.02);

        double[] scanTimes = { 14.5, 15.0, 15.5 };
        double[] mzGrid = { 1000.0, 1200.0, 1500.0 };
        double[,] grid = fm.Rasterize(scanTimes, mzGrid);

        for (int s = 0; s < scanTimes.Length; s++)
        for (int b = 0; b < mzGrid.Length; b++)
            Assert.That(grid[s, b], Is.EqualTo(fm.Evaluate(scanTimes[s], mzGrid[b])).Within(1e-12));
    }

    [Test]
    public void TwoProteoformsSumIndependently()
    {
        var a = MakeModel(mass: 10000.0, rtCenter: 20.0, muZ: 8, abundance: 1.0);
        var b = MakeModel(mass: 12000.0, rtCenter: 20.0, muZ: 9, abundance: 1.0);
        var both = new ForwardModel(new[] { a, b }, minCharge: 5, maxCharge: 14, sigmaMz: 0.02);
        var onlyA = new ForwardModel(new[] { a }, minCharge: 5, maxCharge: 14, sigmaMz: 0.02);
        var onlyB = new ForwardModel(new[] { b }, minCharge: 5, maxCharge: 14, sigmaMz: 0.02);

        var kernelA = new IsotopeEnvelopeKernel(10000.0);
        double mz = kernelA.CentroidMzs(8)[0];

        Assert.That(
            both.Evaluate(20.0, mz),
            Is.EqualTo(onlyA.Evaluate(20.0, mz) + onlyB.Evaluate(20.0, mz)).Within(1e-9));
    }
}
