using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests;
[ExcludeFromCodeCoverage]
public class TestPixelStackRejection
{
    private SpectralAveragingOptions options;
    static double[] minMaxTest = { 10, 9, 8, 7, 6, 5 };
    static double[] minMaxExpected = { 9, 8, 7, 6 };

    static double[] percentileTest = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
    static double[] percentileExpected = new double[] { 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5 };

    static double[] sigmaTest = { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
    static double[] sigmaExpected = { 60, 50, 40, 30, 20, 10, 0 };

    static double[] winsorizedTest = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    static double[] winsorizedExpected = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

    static double[] testAveragedSigma = { 120, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 30, 15 };
    static double[] expectedAveragedSigma = testAveragedSigma[1..^1];

    private static double[] testThreshold = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    private static double[] expectedThreshold = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

    private static readonly object[] arguments =
    {
        new object[] { new TestCase(RejectionType.PercentileClipping, percentileTest, percentileExpected) },
        new object[] { new TestCase(RejectionType.WinsorizedSigmaClipping, winsorizedTest, winsorizedExpected) },
        new object[] { new TestCase(RejectionType.AveragedSigmaClipping, testAveragedSigma, expectedAveragedSigma) },
        new object[] { new TestCase(RejectionType.MinMaxClipping, minMaxTest, minMaxExpected) },
        new object[] { new TestCase(RejectionType.SigmaClipping, sigmaTest, sigmaExpected) },
        new object[] { new TestCase(RejectionType.NoRejection, minMaxTest, minMaxTest) }, 
        new object[] { new TestCase(RejectionType.BelowThresholdRejection, testThreshold, expectedThreshold)}
    };

    public class TestCase
    {
        public RejectionType RejectionType { get; set; }
        public double[] TestArray { get; set; }
        public double[] ExpectedArray { get; set; }

        public TestCase(RejectionType rejection, double[] test, double[] expected)
        {
            RejectionType = rejection;
            TestArray = test;
            ExpectedArray = expected; 
        }
    }
    


    [OneTimeSetUp]
    public void OneTimeSetup()
    {
        options = new();
        options.SetDefaultValues();
        options.MaxSigmaValue = 1.5;
        options.MinSigmaValue = 1.5;
        options.BinSize = 1.0d;
        options.PerformNormalization = true;
        options.Percentile = 0.9;
        options.SpectrumMergingType = SpectrumMergingType.MrsNoiseEstimate;
    }


    [Test]
    [TestCaseSource("arguments")]
    public void TestRejectionMethods(TestCase testCase)
    {
        var testArray = testCase.TestArray;
        var expectedArray = testCase.ExpectedArray;
        var rejection = testCase.RejectionType; 

        double[] testXarray = Enumerable.Repeat(1.0, testArray.Length).ToArray(); 

        PixelStack pixelStack = new(testXarray, testArray); 
        options.RejectionType = rejection;
        pixelStack.PerformRejection(options);
        
        Assert.That(pixelStack.UnrejectedIntensities, Is.EqualTo(expectedArray));
    }
}