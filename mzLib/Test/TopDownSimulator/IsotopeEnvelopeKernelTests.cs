using System;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class IsotopeEnvelopeKernelTests
{
    [Test]
    public void IsotopologuesAreOrderedByAscendingMassAtEveryCharge()
    {
        // Averagine's own tables are ordered most-intense-first; the kernel must not inherit that,
        // because callers index isotopologues by position and take neighbours as neighbours in m/z.
        foreach (double mass in new[] { 2500.0, 8000.0, 14800.0, 30000.0 })
        {
            var kernel = new IsotopeEnvelopeKernel(mass);
            Assert.That(kernel.IsotopologueCount, Is.GreaterThan(2));

            for (int i = 1; i < kernel.IsotopologueCount; i++)
            {
                Assert.That(kernel.NeutralMass(i), Is.GreaterThan(kernel.NeutralMass(i - 1)),
                    $"Neutral masses are out of order at index {i} for mass {mass}.");
            }

            foreach (int z in new[] { 1, 8, 20, 29 })
            {
                double[] centroids = kernel.CentroidMzs(z);
                for (int i = 1; i < centroids.Length; i++)
                {
                    Assert.That(centroids[i], Is.GreaterThan(centroids[i - 1]),
                        $"Centroid m/z values are out of order at index {i}, charge {z}, mass {mass}.");
                }

                var (minMz, maxMz) = kernel.GetMzBounds(z);
                Assert.That(minMz, Is.EqualTo(centroids[0]));
                Assert.That(maxMz, Is.EqualTo(centroids[^1]));
            }
        }
    }

    [Test]
    public void ReorderingKeepsEachMassPairedWithItsOwnIntensity()
    {
        const double mass = 14800.0;
        var kernel = new IsotopeEnvelopeKernel(mass);

        // Rebuild the expected (mass, weight) pairs straight from the averagine table the kernel
        // drew them from, so the pairing is checked against the source rather than against itself.
        int idx = new Averagine().GetMostIntenseMassIndex(mass);
        double[] theoreticalMasses = Averagine.AllMasses[idx];
        double[] theoreticalIntensities = Averagine.AllIntensities[idx];
        double shift = mass - (theoreticalMasses[0] - Averagine.DiffToMonoisotopic[idx]);
        double intensitySum = theoreticalIntensities.Sum();

        var expected = theoreticalMasses
            .Select((m, i) => (Mass: m + shift, Weight: theoreticalIntensities[i] / intensitySum))
            .OrderBy(pair => pair.Mass)
            .ToArray();

        Assert.That(kernel.IsotopologueCount, Is.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.That(kernel.NeutralMass(i), Is.EqualTo(expected[i].Mass).Within(1e-12));
            Assert.That(kernel.Intensity(i), Is.EqualTo(expected[i].Weight).Within(1e-15),
                $"Isotopologue at {expected[i].Mass} lost its intensity in the reordering.");
        }

        // The brightest isotopologue is still the averagine table's most-intense mass — it has just
        // moved out of index 0 and into the middle of the envelope.
        int brightest = 0;
        for (int i = 1; i < kernel.IsotopologueCount; i++)
        {
            if (kernel.Intensity(i) > kernel.Intensity(brightest))
                brightest = i;
        }

        Assert.That(brightest, Is.GreaterThan(0),
            "At 14.8 kDa the monoisotopic peak is not the tallest, so index 0 must not be the brightest.");
        Assert.That(kernel.NeutralMass(brightest),
            Is.EqualTo(Averagine.MostIntenseMasses[idx] + shift).Within(1e-9));

        Assert.That(Enumerable.Range(0, kernel.IsotopologueCount).Sum(kernel.Intensity),
            Is.EqualTo(1.0).Within(1e-9));
    }

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
