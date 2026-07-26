using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class PeakWidthModelTests
{
    private const double K = 3.5e-7;

    [Test]
    public void OrbitrapWidthReproducesTheClosedFormAndIncreasesWithMz()
    {
        var model = new OrbitrapPeakWidth(K);

        foreach (double mz in new[] { 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0 })
            Assert.That(model.SigmaAt(mz), Is.EqualTo(K * Math.Pow(mz, 1.5)).Within(1e-15),
                $"σ at m/z {mz} does not match K·(m/z)^1.5.");

        double previous = 0;
        for (double mz = 200; mz <= 2000; mz += 50)
        {
            double sigma = model.SigmaAt(mz);
            Assert.That(sigma, Is.GreaterThan(previous), $"σ is not increasing at m/z {mz}.");
            previous = sigma;
        }
    }

    [Test]
    public void OrbitrapWidthBoundsItsSpanAtTheTopEnd()
    {
        var model = new OrbitrapPeakWidth(K);

        Assert.That(model.MaxSigmaOver(600, 1500), Is.EqualTo(model.SigmaAt(1500)).Within(1e-15));
        Assert.That(model.MaxSigmaOver(1500, 600), Is.EqualTo(model.SigmaAt(1500)).Within(1e-15),
            "The bound must not depend on argument order.");
    }

    [Test]
    public void ConstantWidthIsFlatAndBoundsItselfExactly()
    {
        var model = new ConstantPeakWidth(0.0072);

        Assert.That(model.SigmaAt(400), Is.EqualTo(0.0072));
        Assert.That(model.SigmaAt(2000), Is.EqualTo(0.0072));
        Assert.That(model.MaxSigmaOver(400, 2000), Is.EqualTo(0.0072));
    }

    [Test]
    public void NonPositiveWidthParametersAreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new ConstantPeakWidth(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ConstantPeakWidth(-1e-3));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ConstantPeakWidth(double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => new OrbitrapPeakWidth(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new OrbitrapPeakWidth(double.PositiveInfinity));
    }

    [Test]
    public void FromSigmaAtPassesThroughItsObservation()
    {
        var model = OrbitrapPeakWidth.FromSigmaAt(mz: 780.0, sigmaMz: 0.00721);
        Assert.That(model.SigmaAt(780.0), Is.EqualTo(0.00721).Within(1e-15));
    }

    [Test]
    public void PooledFitRecoversAKnownWidthLaw()
    {
        var truth = new OrbitrapPeakWidth(K);
        var measurements = SyntheticMeasurements(truth, jitter: 0.0);

        var fit = new PeakWidthModelFitter().Fit(measurements);

        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.K, Is.EqualTo(K).Within(1e-3).Percent);
        Assert.That(fit.FreeSlope, Is.EqualTo(OrbitrapPeakWidth.Exponent).Within(1e-6),
            "Noise-free data on the law itself must reproduce a slope of exactly 1.5.");
        Assert.That(fit.Clamped, Is.False);
    }

    [Test]
    public void PooledFitSurvivesScatterThatWouldDominateASingleRecord()
    {
        // Per-record sigma estimates on real data scatter by ~3x. Pooling is what makes k
        // recoverable anyway; a median of noisy scalars cannot see the m/z dependence at all.
        var truth = new OrbitrapPeakWidth(K);
        var measurements = SyntheticMeasurements(truth, jitter: 0.30);

        var fit = new PeakWidthModelFitter().Fit(measurements);

        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.K, Is.EqualTo(K).Within(15).Percent);
        Assert.That(fit.FreeSlope, Is.EqualTo(OrbitrapPeakWidth.Exponent).Within(0.35));
    }

    [Test]
    public void OneWildRecordDoesNotDistortTheFit()
    {
        var truth = new OrbitrapPeakWidth(K);
        var measurements = SyntheticMeasurements(truth, jitter: 0.05).ToList();

        // A record whose width estimate came out 50x too wide — the shape of a genuinely bad fit.
        for (int i = 0; i < 12; i++)
            measurements.Add(new PeakWidthMeasurement(900.0 + i, 50 * truth.SigmaAt(900.0), 1e6, 40));

        var withOutliers = new PeakWidthModelFitter().Fit(measurements);
        var withoutRejection = new PeakWidthModelFitter(
            new PeakWidthModelFitOptions(OutlierMadMultiple: 0)).Fit(measurements);

        Assert.That(withOutliers, Is.Not.Null);
        Assert.That(withoutRejection, Is.Not.Null);
        Assert.That(withOutliers!.MeasurementsUsed, Is.LessThan(withOutliers.MeasurementCount),
            "The outlier pass should have dropped something.");
        Assert.That(withOutliers.K, Is.EqualTo(K).Within(10).Percent);
        Assert.That(Math.Abs(withoutRejection!.K - K), Is.GreaterThan(Math.Abs(withOutliers.K - K)),
            "Rejecting outliers must land closer to the truth than keeping them.");
    }

    [Test]
    public void FitIsRefusedRatherThanExtrapolatedFromTooFewMeasurements()
    {
        var truth = new OrbitrapPeakWidth(K);
        var measurements = new[]
        {
            new PeakWidthMeasurement(800.0, truth.SigmaAt(800.0), 1e5, 40),
            new PeakWidthMeasurement(810.0, truth.SigmaAt(810.0), 1e5, 40),
        };

        Assert.That(new PeakWidthModelFitter().Fit(measurements), Is.Null);
    }

    [Test]
    public void FitIsRefusedWhenTheMzSpanCannotConstrainTheExponent()
    {
        // Plenty of measurements, all crammed into a narrow band. A power law fitted there is an
        // extrapolation dressed as a fit — the exponent is barely constrained, so k inherits
        // whatever the scatter happened to do. This is the shape of what a centroided run produces
        // once the coarse-sampling guard has thinned it: m/z 601-725, and a k ~4x an Orbitrap's.
        var truth = new OrbitrapPeakWidth(K);
        var narrow = Enumerable.Range(0, 60)
            .Select(i => 601.0 + i * 2.0)
            .Select(mz => new PeakWidthMeasurement(mz, truth.SigmaAt(mz), 1e5, 40))
            .ToArray();

        Assert.That(narrow, Has.Length.GreaterThanOrEqualTo(30), "Measurement count must not be the reason for refusal.");
        Assert.That(narrow.Max(m => m.Mz) / narrow.Min(m => m.Mz), Is.LessThan(1.5));
        Assert.That(new PeakWidthModelFitter().Fit(narrow), Is.Null);

        // The same measurements over a real span are fine, so it is the span that is disqualifying.
        var wide = Enumerable.Range(0, 60)
            .Select(i => 600.0 + i * 15.0)
            .Select(mz => new PeakWidthMeasurement(mz, truth.SigmaAt(mz), 1e5, 40))
            .ToArray();

        var fit = new PeakWidthModelFitter().Fit(wide);
        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.K, Is.EqualTo(K).Within(1e-3).Percent);
    }

    [Test]
    public void StandardErrorDoesNotDependOnHowFinelyEachWindowWasSampled()
    {
        // A window yields ONE sigma measurement however many profile points it was sampled at.
        // Treating its weight as a replication count shrinks the standard error by about
        // sqrt(samples), which would make every decision built on it a function of the detector's
        // sampling rate rather than of the data.
        var truth = new OrbitrapPeakWidth(K);
        var random = new Random(7);
        var offsets = Enumerable.Range(0, 120).Select(_ => 0.15 * (2 * random.NextDouble() - 1)).ToArray();

        PeakWidthModelFit FitWithPeakCount(int peakCount)
        {
            var measurements = offsets
                .Select((offset, i) =>
                {
                    double mz = 500.0 + i * 9.0;
                    return new PeakWidthMeasurement(mz, truth.SigmaAt(mz) * Math.Exp(offset), 1e5, peakCount);
                })
                .ToArray();

            var fit = new PeakWidthModelFitter().Fit(measurements);
            Assert.That(fit, Is.Not.Null);
            return fit!;
        }

        var sparse = FitWithPeakCount(4);
        var dense = FitWithPeakCount(60);

        Assert.That(dense.FreeSlope, Is.EqualTo(sparse.FreeSlope).Within(1e-9),
            "Uniform weights cancel in the point estimate, so only the standard error is at issue.");
        Assert.That(dense.FreeSlopeStandardError, Is.EqualTo(sparse.FreeSlopeStandardError).Within(1e-6).Percent,
            "The standard error must reflect the number of measurements, not the sampling rate.");
    }

    [Test]
    public void APerfectFitIsNotReportedAsContradictingTheWidthLaw()
    {
        // Noise-free data lies exactly on the law, so residuals are rounding error and the standard
        // error collapses toward zero. A ratio test alone would then read a tiny slope wobble as an
        // enormous number of standard errors and reject the best possible fit.
        var fit = new PeakWidthModelFitter().Fit(SyntheticMeasurements(new OrbitrapPeakWidth(K), jitter: 0.0));

        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.SlopeDeviation, Is.LessThan(1e-6));
        Assert.That(fit.ContradictsWidthLaw(), Is.False,
            $"A perfect fit was rejected: slope {fit.FreeSlope}, SE {fit.FreeSlopeStandardError}, " +
            $"{fit.SlopeDeviationInStandardErrors} standard errors.");
    }

    [Test]
    public void AWrongExponentIsReportedAsContradictingTheWidthLaw()
    {
        // sigma proportional to (m/z)^1.0 is the signature of width measured off centroided input:
        // the extraction window is half the isotopologue spacing, 0.5/z, and m/z ~ M/z, so centroid
        // scatter fills a window in direct proportion to m/z. This must be caught.
        var random = new Random(11);
        var linear = Enumerable.Range(0, 120)
            .Select(i =>
            {
                double mz = 500.0 + i * 9.0;
                double sigma = 1.2e-5 * mz * Math.Exp(0.10 * (2 * random.NextDouble() - 1));
                return new PeakWidthMeasurement(mz, sigma, 1e5, 40);
            })
            .ToArray();

        var fit = new PeakWidthModelFitter().Fit(linear);

        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.FreeSlope, Is.EqualTo(1.0).Within(0.15));
        Assert.That(fit.ContradictsWidthLaw(), Is.True,
            $"A slope of {fit.FreeSlope:F3} +/- {fit.FreeSlopeStandardError:F3} should have been rejected.");
    }

    [Test]
    public void OutlierRejectionRefusesRatherThanReadmittingTheMeasurementsItDropped()
    {
        // If too few survive the outlier pass, readmitting them would let precisely the
        // measurements identified as bad set k, while the log reported "n/n used" — indistinguishable
        // from having found no outliers at all.
        var truth = new OrbitrapPeakWidth(K);
        var random = new Random(3);
        var measurements = Enumerable.Range(0, 34)
            .Select(i =>
            {
                double mz = 500.0 + i * 33.0;
                double jitter = Math.Exp(0.08 * (2 * random.NextDouble() - 1));
                // 20 inliers, 14 wildly wrong — a clear majority survives the outlier pass, but
                // 20 is still short of MinimumMeasurements.
                double sigma = i < 20 ? truth.SigmaAt(mz) * jitter : truth.SigmaAt(mz) * 50 * jitter;
                return new PeakWidthMeasurement(mz, sigma, 1e5, 40);
            })
            .ToArray();

        Assert.That(new PeakWidthModelFitter().Fit(measurements), Is.Null);

        // Lowering the bar admits the same fit, confirming the count is what disqualified it rather
        // than some other guard.
        var permissive = new PeakWidthModelFitter(new PeakWidthModelFitOptions(MinimumMeasurements: 15));
        var fit = permissive.Fit(measurements);
        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.MeasurementsUsed, Is.EqualTo(20), "Exactly the inliers should have survived.");
        Assert.That(fit.MeasurementsUsed, Is.LessThan(fit.MeasurementCount),
            "Survivors must be reported as fewer than the total, not silently readmitted.");
    }

    [Test]
    public void SpanGuardIsAppliedToTheMeasurementsTheFitActuallyUsed()
    {
        // A handful of far-out measurements can satisfy the lever-arm guard and then be discarded
        // as outliers, leaving the narrow-band extrapolation the guard exists to refuse.
        var truth = new OrbitrapPeakWidth(K);
        var random = new Random(5);
        var narrowBand = Enumerable.Range(0, 100)
            .Select(i => 600.0 + i)
            // Real scatter, so the MAD is non-zero and the outlier pass has something to work with.
            .Select(mz => new PeakWidthMeasurement(
                mz, truth.SigmaAt(mz) * Math.Exp(0.06 * (2 * random.NextDouble() - 1)), 1e5, 40));

        // Far away in m/z and far off the law, so the outlier pass will drop exactly these.
        var farOutliers = Enumerable.Range(0, 10)
            .Select(i => new PeakWidthMeasurement(1500.0 + i, truth.SigmaAt(1500.0) * 40, 1e5, 40));

        var combined = narrowBand.Concat(farOutliers).ToArray();
        Assert.That(combined.Max(m => m.Mz) / combined.Min(m => m.Mz), Is.GreaterThan(1.5),
            "Before rejection the span looks wide enough, which is the trap being tested.");

        Assert.That(new PeakWidthModelFitter().Fit(combined), Is.Null);
    }

    [Test]
    public void ImplausibleFitIsClampedRatherThanShipped()
    {
        // Widths an order of magnitude past anything an Orbitrap produces.
        var absurd = Enumerable.Range(0, 40)
            .Select(i => new PeakWidthMeasurement(600.0 + i * 20, 5.0, 1e5, 40))
            .ToArray();

        var fit = new PeakWidthModelFitter().Fit(absurd);

        Assert.That(fit, Is.Not.Null);
        Assert.That(fit!.Clamped, Is.True);
        Assert.That(fit.Model.SigmaAt(800.0), Is.EqualTo(0.2).Within(1e-9));
    }

    /// <summary>
    /// Width measurements drawn from <paramref name="law"/> across a realistic acquisition range,
    /// multiplicatively perturbed by a deterministic pseudo-random factor.
    /// </summary>
    private static PeakWidthMeasurement[] SyntheticMeasurements(IPeakWidthModel law, double jitter)
    {
        var random = new Random(20260726);
        var measurements = new List<PeakWidthMeasurement>();

        for (double mz = 500; mz <= 1600; mz += 25)
        {
            for (int rep = 0; rep < 3; rep++)
            {
                double factor = jitter > 0 ? Math.Exp(jitter * (2 * random.NextDouble() - 1)) : 1.0;
                measurements.Add(new PeakWidthMeasurement(mz, law.SigmaAt(mz) * factor, 1e5, 40));
            }
        }

        return measurements.ToArray();
    }

    [Test]
    public void SimulatedPeaksAreWiderAtHighMzThanAtLow()
    {
        // The property a constant σ cannot have, measured on rendered profile data rather than
        // asserted against the model that produced it.
        var width = new OrbitrapPeakWidth(K);
        var kernel = new IsotopeEnvelopeKernel(10000.0);

        double lowMz = kernel.CentroidMzs(17)[0];
        double highMz = kernel.CentroidMzs(7)[0];
        Assert.That(lowMz, Is.LessThan(700), "Fixture expects the high-charge envelope near m/z 600.");
        Assert.That(highMz, Is.GreaterThan(1400), "Fixture expects the low-charge envelope near m/z 1500.");

        double lowFwhm = MeasureFwhm(kernel, charge: 17, width);
        double highFwhm = MeasureFwhm(kernel, charge: 7, width);

        Assert.That(highFwhm, Is.GreaterThan(lowFwhm * 2),
            $"FWHM at m/z {highMz:F0} ({highFwhm:F4}) should dwarf the one at m/z {lowMz:F0} ({lowFwhm:F4}).");
    }

    [Test]
    public void HighMzEnvelopesMergeUnderTheWidthLawButNotUnderAConstant()
    {
        // The whole point of the change: a real 25 kDa proteoform at high charge and high m/z has
        // isotopologues that partly merge, and handling that is what separates a good feature
        // finder from a naive one. A constant σ never produces the regime at all.
        const double mass = 25000.0;
        const int charge = 25;
        var kernel = new IsotopeEnvelopeKernel(mass);

        double apexMz = kernel.CentroidMzs(charge)[0];
        Assert.That(apexMz, Is.GreaterThan(950).And.LessThan(1100),
            "Fixture expects this proteoform's charge-25 envelope near m/z 1000.");

        // Chosen so both models give the same width at m/z 1000; only their m/z dependence differs.
        var law = OrbitrapPeakWidth.FromSigmaAt(1000.0, 0.0111);
        var constant = new ConstantPeakWidth(0.0111);

        int lowMzMaxima = CountLocalMaxima(kernel, charge, law);
        int highMzMaxima = CountLocalMaxima(new IsotopeEnvelopeKernel(mass * 1.6), 25, law);
        int highMzMaximaConstant = CountLocalMaxima(new IsotopeEnvelopeKernel(mass * 1.6), 25, constant);

        Assert.That(highMzMaxima, Is.LessThan(lowMzMaxima),
            "Under the width law, the same charge at higher m/z must resolve into fewer peaks.");
        Assert.That(highMzMaximaConstant, Is.GreaterThan(highMzMaxima),
            "A constant σ keeps the high-m/z envelope resolved — the defect this change fixes.");
    }

    /// <summary>Full width at half maximum of the monoisotopic peak, measured off a dense scan.</summary>
    private static double MeasureFwhm(IsotopeEnvelopeKernel kernel, int charge, IPeakWidthModel width)
    {
        double centre = kernel.CentroidMzs(charge)[0];
        double apex = kernel.Evaluate(centre, charge, width);
        double half = apex / 2;
        double step = width.SigmaAt(centre) / 200.0;

        double left = centre;
        while (kernel.Evaluate(left, charge, width) > half && centre - left < 10)
            left -= step;

        double right = centre;
        while (kernel.Evaluate(right, charge, width) > half && right - centre < 10)
            right += step;

        return right - left;
    }

    /// <summary>
    /// Distinct local maxima across the envelope, sampled finely enough that merging — not
    /// sampling — decides the count.
    /// </summary>
    private static int CountLocalMaxima(IsotopeEnvelopeKernel kernel, int charge, IPeakWidthModel width)
    {
        var (minMz, maxMz) = kernel.GetMzBounds(charge);
        double step = width.SigmaAt(minMz) / 50.0;
        var values = new List<double>();
        for (double mz = minMz - 3 * step; mz <= maxMz + 3 * step; mz += step)
            values.Add(kernel.Evaluate(mz, charge, width));

        double peak = values.Max();
        int maxima = 0;
        for (int i = 1; i < values.Count - 1; i++)
        {
            // Ignore the noise floor at the envelope's shoulders; only count structure that would
            // plausibly be centroided out of real data.
            if (values[i] < peak * 1e-3) continue;
            if (values[i] > values[i - 1] && values[i] >= values[i + 1])
                maxima++;
        }

        return maxima;
    }

    [Test]
    public void AbundanceFittedUnderAWidthModelRendersAtThatHeight()
    {
        // Pins the fitting/rendering mismatch closed: the peak heights a simulation produces must
        // be the ones the abundance was fitted to, not off by σ_fitted/σ_rendered.
        var width = new OrbitrapPeakWidth(K);
        var model = new ProteoformModel(
            MonoisotopicMass: 25000.0,
            Abundance: 1.5e6,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 24, SigmaZ: 2.0),
            Identifier: "P00000|ROUNDTRIP");

        double[] scanTimes = Enumerable.Range(0, 21).Select(i => 19.0 + i * 0.1).ToArray();
        var simulation = new Simulator().SimulateCentroided(new[] { model }, 20, 28, width, scanTimes);

        var extractor = new GroundTruthExtractor(simulation.Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var truth = extractor.Extract(model.MonoisotopicMass, rtCenter: 20.0, rtHalfWidth: 1.5,
            minCharge: 20, maxCharge: 28);

        var fitted = new AbundanceFitter().Fit(truth, width, model.RtProfile, model.ChargeDistribution);

        Assert.That(fitted.Abundance, Is.EqualTo(model.Abundance).Within(0.5).Percent);
    }
}
