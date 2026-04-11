using System;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class IsotopeEnvelopeKernelTests
{
    [Test]
    public void IntensitiesNormalizedToUnity()
    {
        var kernel = new IsotopeEnvelopeKernel(monoisotopicMass: 15000.0);
        double sum = 0;
        for (int i = 0; i < kernel.IsotopologueCount; i++)
            sum += kernel.Intensity(i);
        Assert.That(sum, Is.EqualTo(1.0).Within(1e-9));
    }

    [Test]
    public void CentroidMzsMatchChargeTransform()
    {
        const double mass = 10000.0;
        var kernel = new IsotopeEnvelopeKernel(mass);
        int z = 8;
        double[] centroids = kernel.CentroidMzs(z);
        for (int i = 0; i < centroids.Length; i++)
        {
            double expected = kernel.NeutralMass(i).ToMz(z);
            Assert.That(centroids[i], Is.EqualTo(expected).Within(1e-12));
        }
    }

    [Test]
    public void EvaluatePeaksAtCentroidPositions()
    {
        var kernel = new IsotopeEnvelopeKernel(monoisotopicMass: 12000.0);
        int z = 10;
        double[] centroids = kernel.CentroidMzs(z);
        double sigmaMz = 0.01;

        // "Far" must be well away from every centroid (at z=10 centroids are ~0.1 m/z apart
        // and the envelope spans ~5 m/z), so probe 5 m/z below the monoisotopic centroid.
        double onCentroid = kernel.Evaluate(centroids[0], z, sigmaMz);
        double farOff = kernel.Evaluate(centroids[0] - 5.0, z, sigmaMz);
        Assert.That(onCentroid, Is.GreaterThan(farOff * 1e6));
    }

    [Test]
    public void IntegralOverMzEqualsTotalIntensityDivZ()
    {
        // Each isotopologue's Gaussian contributes ∫ w_i/(σ√2π)·exp(-x²/(2σ²)) dx = w_i.
        // Total ∫ φ(mz; M, z, σ) d(mz) ≈ Σ w_i = 1 because weights sum to 1.
        var kernel = new IsotopeEnvelopeKernel(monoisotopicMass: 8000.0);
        int z = 6;
        double sigma = 0.02;
        double[] centroids = kernel.CentroidMzs(z);
        double lo = centroids.Min() - 0.5;
        double hi = centroids.Max() + 0.5;
        double integral = TrapezoidalIntegrate(mz => kernel.Evaluate(mz, z, sigma), lo, hi, 40000);
        Assert.That(integral, Is.EqualTo(1.0).Within(2e-3));
    }

    private static double TrapezoidalIntegrate(Func<double, double> f, double a, double b, int steps)
    {
        double h = (b - a) / steps;
        double sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < steps; i++)
            sum += f(a + i * h);
        return sum * h;
    }
}
