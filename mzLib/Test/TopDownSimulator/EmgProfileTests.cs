using System;
using NUnit.Framework;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class EmgProfileTests
{
    [Test]
    public void PdfIntegratesToOne()
    {
        var emg = new EmgProfile(Mu: 20.0, Sigma: 0.3, Tau: 0.2);
        double integral = TrapezoidalIntegrate(emg.Evaluate, 15.0, 25.0, 20000);
        Assert.That(integral, Is.EqualTo(1.0).Within(1e-3));
    }

    [Test]
    public void TauZeroFallsBackToGaussian()
    {
        var emg = new EmgProfile(Mu: 10.0, Sigma: 0.5, Tau: 0.0);
        double valueAtMu = emg.Evaluate(10.0);
        double expected = 1.0 / (0.5 * Math.Sqrt(2 * Math.PI));
        Assert.That(valueAtMu, Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void HasRightwardTailWhenTauPositive()
    {
        var emg = new EmgProfile(Mu: 10.0, Sigma: 0.3, Tau: 0.5);
        // With tau > 0, the profile is skewed right: value at mu + delta
        // should exceed value at mu - delta for sufficiently large delta.
        Assert.That(emg.Evaluate(10.8), Is.GreaterThan(emg.Evaluate(9.2)));
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
