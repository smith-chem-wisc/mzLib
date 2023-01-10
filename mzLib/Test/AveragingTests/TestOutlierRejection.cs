using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;
using NUnit.Framework.Internal;

namespace Test.AveragingTests;
[ExcludeFromCodeCoverage]
public class TestOutlierRejection
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
    private static double[] expectedThreshold = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

    private static readonly object[] arguments =
    {
        new object[] { new TestCase(OutlierRejectionType.PercentileClipping, percentileTest, percentileExpected) },
        new object[] { new TestCase(OutlierRejectionType.WinsorizedSigmaClipping, winsorizedTest, winsorizedExpected) },
        new object[] { new TestCase(OutlierRejectionType.AveragedSigmaClipping, testAveragedSigma, expectedAveragedSigma) },
        new object[] { new TestCase(OutlierRejectionType.MinMaxClipping, minMaxTest, minMaxExpected) },
        new object[] { new TestCase(OutlierRejectionType.SigmaClipping, sigmaTest, sigmaExpected) },
        new object[] { new TestCase(OutlierRejectionType.NoRejection, minMaxTest, minMaxTest) }, 
        new object[] { new TestCase(OutlierRejectionType.BelowThresholdRejection, testThreshold, expectedThreshold)}
    };

    public class TestCase
    {
        public OutlierRejectionType OutlierRejectionType { get; set; }
        public double[] TestArray { get; set; }
        public double[] ExpectedArray { get; set; }

        public TestCase(OutlierRejectionType outlierRejection, double[] test, double[] expected)
        {
            OutlierRejectionType = outlierRejection;
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
        options.SpectraMergingType = SpectraMergingType.MrsNoiseEstimate;
    }

    [Test]
    [TestCaseSource("arguments")]
    public void TestRejectionMethods(TestCase testCase)
    {
        // test pixel stack entry point
        var testArray = testCase.TestArray;
        var expectedArray = testCase.ExpectedArray;
        var rejection = testCase.OutlierRejectionType; 

        double[] testXarray = Enumerable.Repeat(1.0, testArray.Length).ToArray(); 

        PeakBin peakBin = new(testXarray, testArray); 
        options.OutlierRejectionType = rejection;
        OutlierRejection.RejectOutliers(peakBin, options);
        
        Assert.That(peakBin.UnrejectedIntensities, Is.EqualTo(expectedArray));

        // test double array entry point
        var results = OutlierRejection.RejectOutliers(testArray, options);
        Assert.That(results, Is.EqualTo(expectedArray));

        // test binned spectra entry point
        double[][] yArrays = new double[testArray.Length][];
        double[][] xArrays = new double[testXarray.Length][];
        for (var i = 0; i < testArray.Length; i++)
        {
            yArrays[i] = new double[1] { testArray[i] };
            xArrays[i] = new double[1] { testXarray[i] };
        }
        BinnedSpectra binnedSpectra = new(xArrays, yArrays, options.BinSize);
        OutlierRejection.RejectOutliers(binnedSpectra, options);
        Assert.That(binnedSpectra.PeakBins.First().UnrejectedIntensities, Is.EqualTo(expectedArray));
    }

    [Test]
    public static void TestMinMaxClipping()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        { OutlierRejectionType = OutlierRejectionType.MinMaxClipping };
        double[] test = { 10, 9, 8, 7, 6, 5 };
        double[] expected = { 9, 8, 7, 6 };
        double[] minMaxClipped = OutlierRejection.RejectOutliers(test, options);
        Assert.That(minMaxClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        double[] testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        PeakBin stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));
    }

    [Test]
    public static void TestPercentileClipping()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        { OutlierRejectionType = OutlierRejectionType.PercentileClipping, Percentile = 0.9 };
        double[] test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
        double[] expected = new double[] { 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5 };
        double[] percentileClipped = OutlierRejection.RejectOutliers(test, options);
        Assert.That(percentileClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        double[] testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        PeakBin stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
        expected = new double[] { 60, 50, 40, 30, 20, 10 };
        percentileClipped = OutlierRejection.RejectOutliers(test, options);
        Assert.That(percentileClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));
    }

    [Test]
    public static void TestSigmaClipping()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        { OutlierRejectionType = OutlierRejectionType.SigmaClipping, MaxSigmaValue = 1.5, MinSigmaValue = 1.5 };
        var test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
        double[] sigmaClipped = OutlierRejection.RejectOutliers(test, options);
        var expected = new double[] { 60, 50, 40, 30, 20, 10, 0 };
        Assert.That(sigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        var testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        PeakBin stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1;
        options.MaxSigmaValue = 1;
        test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
        sigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 60, 50, 40 };
        Assert.That(sigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1.3;
        options.MaxSigmaValue = 1.3;
        sigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 70, 60, 50, 40, 30 };
        Assert.That(sigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1;
        options.MaxSigmaValue = 1.3;
        sigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 80, 70, 60 };
        Assert.That(sigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1.5;
        options.MaxSigmaValue = 1.5;
        sigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
        Assert.That(sigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));
    }

    [Test]
    public static void TestWinsorizedSigmaClipping()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        {
            OutlierRejectionType = OutlierRejectionType.WinsorizedSigmaClipping,
            MinSigmaValue = 1.5,
            MaxSigmaValue = 1.5
        };
        var test = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
        double[] windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, options);
        var expected = new double[] { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
        Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        var testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        PeakBin stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        test = new double[] { 15, 30, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 100 };
        expected = new double[] { 56d, 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
        windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, options);
        Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1.3;
        options.MaxSigmaValue = 1.3;
        windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d };
        Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MaxSigmaValue = 1.5;
        windsorizedSigmaClipped = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 57d, 58d, 59d, 60d, 61d, 62d, 63d, 64d, 65d };
        Assert.That(windsorizedSigmaClipped, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));
    }

    [Test]
    public static void TestAveragedSigmaClipping()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        {
            OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping,
            MaxSigmaValue = 2,
            MinSigmaValue = 2
        };
        var test = new double[] { 120, 65, 64, 63, 62, 61, 60, 59, 59, 58, 57, 56, 30, 15 };
        double[] averagedSigmaClipping = OutlierRejection.RejectOutliers(test, options);
        var expected = test[1..test.Length];
        Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

        // test pixelstack entry point
        var testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        PeakBin stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1;
        options.MaxSigmaValue = 1;
        averagedSigmaClipping = OutlierRejection.RejectOutliers(test, options);
        expected = test[1..^2];
        Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1.3;
        options.MaxSigmaValue = 1.3;
        averagedSigmaClipping = OutlierRejection.RejectOutliers(test, options);
        expected = test[1..^2];
        Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));


        // test array entry point
        options.MinSigmaValue = 1.5;
        options.MaxSigmaValue = 1.5;
        averagedSigmaClipping = OutlierRejection.RejectOutliers(test, options);
        expected = test[1..^1];
        Assert.That(averagedSigmaClipping, Is.EqualTo(expected));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, test.Length).ToArray();
        stack = new(testXarray, test);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(expected));
    }

    [Test]
    public static void TestBelowThresholdRejection()
    {
        // test array entry point
        SpectralAveragingOptions options = new()
        {
            OutlierRejectionType = OutlierRejectionType.BelowThresholdRejection
        };
        var arr = new double[] { 0, 10, 0, 0, 0, 0, 0 };
        var arr2 = new double[] { 0, 10, 0, 0, 0 };
        var arr3 = new double[] { 0, 10, 0, 0 };
        var copy3 = new double[arr3.Length];
        Array.Copy(arr3, copy3, arr3.Length);

        // test array entry point
        var result = OutlierRejection.RejectOutliers(arr, options);
        Assert.That(result.All(p => p == 0));

        // test pixelstack entry point
        var testXarray = Enumerable.Repeat(1.0, arr.Length).ToArray();
        PeakBin stack = new(testXarray, arr);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(!stack.UnrejectedIntensities.Any());


        // test array entry point
        result = OutlierRejection.RejectOutliers(arr2, options);
        Assert.That(result.All(p => p == 0));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, arr2.Length).ToArray();
        stack = new(testXarray, arr2);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(!stack.UnrejectedIntensities.Any());


        // test array entry point
        result = OutlierRejection.RejectOutliers(arr3, options);
        Assert.That(result.SequenceEqual(copy3));

        // test pixelstack entry point
        testXarray = Enumerable.Repeat(1.0, arr3.Length).ToArray();
        stack = new(testXarray, arr3);
        OutlierRejection.RejectOutliers(stack, options);
        Assert.That(stack.UnrejectedIntensities, Is.EqualTo(copy3));
    }


    [Test]
    public static void TestSetValuesAndRejectOutliersSwitch()
    {
        SpectralAveragingOptions options = new SpectralAveragingOptions();

        options.SetDefaultValues();
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.NoRejection);
        Assert.That(options.SpectralWeightingType == SpectraWeightingType.None);
        Assert.That(0.1, Is.EqualTo(options.Percentile));
        Assert.That(1.5, Is.EqualTo(options.MinSigmaValue));
        Assert.That(1.5, Is.EqualTo(options.MaxSigmaValue));

        options.SetValues(OutlierRejectionType.MinMaxClipping, SpectraWeightingType.None, SpectraMergingType.MzBinning,
            performNormalization: true, percentile: .8, minSigma: 2, maxSigma: 4);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.MinMaxClipping);
        Assert.That(0.8, Is.EqualTo(options.Percentile));
        Assert.That(2, Is.EqualTo(options.MinSigmaValue));
        Assert.That(4, Is.EqualTo(options.MaxSigmaValue));

        // no rejection
        options.SetDefaultValues();
        double[] test = new double[] { 10, 8, 6, 5, 4, 2, 0 };
        double[] output = OutlierRejection.RejectOutliers(test, options);
        double[] expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
        Assert.That(output, Is.EqualTo(expected));

        // min max
        options.SetValues(outlierRejectionType: OutlierRejectionType.MinMaxClipping);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.MinMaxClipping);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 8, 6, 5, 4, 2 };
        Assert.That(output, Is.EqualTo(expected));

        // percentile
        options.SetValues(outlierRejectionType: OutlierRejectionType.PercentileClipping);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.PercentileClipping);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 6, 5, 4 };
        Assert.That(output, Is.EqualTo(expected));

        // sigma
        options.SetValues(outlierRejectionType: OutlierRejectionType.SigmaClipping);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.SigmaClipping);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
        Assert.That(output, Is.EqualTo(expected));

        // winsorized sigma
        options.SetValues(outlierRejectionType: OutlierRejectionType.WinsorizedSigmaClipping);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.WinsorizedSigmaClipping);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
        Assert.That(output, Is.EqualTo(expected));
        output = OutlierRejection.RejectOutliers(new double[] { }, options);
        Assert.That(output.SequenceEqual(new double[] { }));

        // averaged sigma
        options.SetValues(outlierRejectionType: OutlierRejectionType.AveragedSigmaClipping, minSigma: 1, maxSigma:1);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.AveragedSigmaClipping);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 8, 6, 5, 4, 2 };
        Assert.That(output, Is.EqualTo(expected));

        // below threshold
        options.SetValues(outlierRejectionType: OutlierRejectionType.BelowThresholdRejection);
        Assert.That(options.OutlierRejectionType == OutlierRejectionType.BelowThresholdRejection);
        output = OutlierRejection.RejectOutliers(test, options);
        expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
        Assert.That(output, Is.EqualTo(expected));
    }
}