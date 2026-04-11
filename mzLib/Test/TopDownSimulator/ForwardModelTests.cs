using NUnit.Framework;
using TopDownSimulator.Model;

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
