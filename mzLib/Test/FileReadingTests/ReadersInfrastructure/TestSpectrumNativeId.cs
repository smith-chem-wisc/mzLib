using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="SpectrumNativeId"/> called directly. MzIdentMLResultFile turns every spectrumID through it, so
    /// TestMzIdentMLResultFile covers it in context; these cases pin the shared entry point itself.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestSpectrumNativeId
    {
        [TestCase(null, -1)]
        [TestCase("", -1)]
        [TestCase("scan=42", 42)]
        [TestCase("controllerType=0 controllerNumber=1 scan=7", 7)]
        [TestCase("scan=81409 file=454", 81409)]
        [TestCase("index=0", 1)]
        [TestCase("index=16332", 16333)]
        [TestCase("mzMLid=controllerType=0", -1)]
        [TestCase("index=12x", -1)]
        [TestCase("scan=1,scan=2", 1)]
        public void OneBasedScanNumber_FromTheNativeId(string nativeId, int expected)
        {
            Assert.That(SpectrumNativeId.OneBasedScanNumberOf(nativeId), Is.EqualTo(expected));
        }
    }
}
