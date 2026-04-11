using MzLibUtil;
using NUnit.Framework;
using TopDownEngine.Features;

namespace Test.TopDownEngine.Features;

[TestFixture]
public class FeatureBoxTests
{
    [Test]
    public void ConstructorRoundTrip_PreservesAllValues()
    {
        FeatureBox original = new(
            MzRange: new MzRange(500.1234, 500.5678),
            RtRange: new DoubleRange(12.34, 15.67),
            SeedIntensity: 12345.6,
            TotalIntensity: 67890.1,
            SourceFile: "sample1.raw",
            PeakCount: 42);

        FeatureBox roundTrip = new(
            MzRange: new MzRange(original.MzRange.Minimum, original.MzRange.Maximum),
            RtRange: new DoubleRange(original.RtRange.Minimum, original.RtRange.Maximum),
            SeedIntensity: original.SeedIntensity,
            TotalIntensity: original.TotalIntensity,
            SourceFile: original.SourceFile,
            PeakCount: original.PeakCount);

        Assert.That(roundTrip.MzRange.Minimum, Is.EqualTo(original.MzRange.Minimum));
        Assert.That(roundTrip.MzRange.Maximum, Is.EqualTo(original.MzRange.Maximum));
        Assert.That(roundTrip.RtRange.Minimum, Is.EqualTo(original.RtRange.Minimum));
        Assert.That(roundTrip.RtRange.Maximum, Is.EqualTo(original.RtRange.Maximum));
        Assert.That(roundTrip.SeedIntensity, Is.EqualTo(original.SeedIntensity));
        Assert.That(roundTrip.TotalIntensity, Is.EqualTo(original.TotalIntensity));
        Assert.That(roundTrip.SourceFile, Is.EqualTo(original.SourceFile));
        Assert.That(roundTrip.PeakCount, Is.EqualTo(original.PeakCount));
    }
}
