using Chemistry;
using NUnit.Framework;
using System;
using LogMassTransform = global::TopDownEngine.Envelopes.LogMassTransform;

namespace Test.TopDownEngine.Envelopes;

[TestFixture]
public class LogMassTransformTests
{
    [TestCase(12345.6789, 7)]
    [TestCase(18765.4321, 19)]
    [TestCase(50000.0, 42)]
    public void LogMz_KnownMassAndCharge_MatchesLogMassMinusLogCharge(double mass, int charge)
    {
        double mz = (mass / charge) + Constants.ProtonMass;

        double observed = LogMassTransform.LogMz(mz);
        double expected = Math.Log(mass) - Math.Log(charge);

        Assert.That(Math.Abs(observed - expected), Is.LessThanOrEqualTo(1e-12));
    }

    [Test]
    public void BuildTemplate_ReturnsNegativeLogChargeSeries()
    {
        double[] template = LogMassTransform.BuildTemplate(4);

        Assert.That(template, Has.Length.EqualTo(4));
        Assert.Multiple(() =>
        {
            Assert.That(template[0], Is.EqualTo(-Math.Log(1)).Within(1e-12));
            Assert.That(template[1], Is.EqualTo(-Math.Log(2)).Within(1e-12));
            Assert.That(template[2], Is.EqualTo(-Math.Log(3)).Within(1e-12));
            Assert.That(template[3], Is.EqualTo(-Math.Log(4)).Within(1e-12));
        });
    }
}
