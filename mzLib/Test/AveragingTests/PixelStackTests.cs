using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests;
[ExcludeFromCodeCoverage]
public class PixelStackTests
{
    private double[] xArray;
    private double[] yArray; 
    private PeakBin _peakBin; 
    [OneTimeSetUp]
    public void OneTimeSetup()
    {
        xArray = new[] { 1.0, 1.1, 0.9, 0.99 };
        yArray = new[] { 10.0, 11.0, 20.0, 50.0 };
        _peakBin = new PeakBin(xArray, yArray);
        _peakBin.Reject(_peakBin.Length - 1);
    }

    [Test]
    public void TestGetIntensities()
    {
        Assert.That(_peakBin.Intensities, Is.EqualTo(yArray)); 
    }

    [Test]
    public void TestGetMzValues()
    {
        Assert.That(_peakBin.Mzs, Is.EqualTo(xArray));
    }

    [Test]
    public void TestNonRejectedLength()
    {
        Assert.True(_peakBin.NonRejectedLength == _peakBin.Length - 1);
    }

    [Test]
    public void TestGetUnrejectedValues()
    {
        Assert.That(_peakBin.UnrejectedIntensities, Is.EqualTo(yArray[..^1]));
        Assert.That(_peakBin.UnrejectedMzs, Is.EqualTo(xArray[..^1]));
    }

    [Test]
    public void TestIsIndexRejected()
    {
        Assert.That(_peakBin.IsIndexRejected(_peakBin.Length - 1));
    }

    [Test]
    public void TestGetIntensityAtIndex()
    {
        Assert.AreEqual(_peakBin.GetIntensityAtIndex(0), yArray[0]); 
    }

    [Test]
    public void TestGetMzAtIndex()
    {
        Assert.AreEqual(_peakBin.GetMzAtIndex(0), xArray[0]);
    }

    [Test]
    public void TestMzAverage()
    {
        Assert.True(_peakBin.UnrejectedMzAverage == xArray[..^1].Average()); 
    }

    [Test]
    public void TestAverage()
    {
        Dictionary<int, double> weights = new(); 
        weights.Add(0, 0.2);
        weights.Add(1, 0.4);
        weights.Add(2, 0.8);
        weights.Add(3, 1.0);
        _peakBin.Average(weights);
        Assert.That(_peakBin.MergedIntensityValue,
            Is.EqualTo(15.9999).Within(0.01));
    }

    [Test]
    [TestCase(OutlierRejectionType.WinsorizedSigmaClipping)]
    [TestCase(OutlierRejectionType.SigmaClipping)]
    [TestCase(OutlierRejectionType.AveragedSigmaClipping)]
    [TestCase(OutlierRejectionType.MinMaxClipping)]
    [TestCase(OutlierRejectionType.PercentileClipping)]
    [TestCase(OutlierRejectionType.BelowThresholdRejection)]
    [TestCase(OutlierRejectionType.NoRejection)]
    public void TestEmptyPixelStack(OutlierRejectionType outlierRejection)
    {
        SpectralAveragingOptions options = new(); 
        options.SetDefaultValues();
        options.OutlierRejectionType = outlierRejection; 
        List<double> emptyXarray = new();
        List<double> emptyYarray = new(); 
        PeakBin emptyStack = new(emptyXarray, emptyYarray);
        OutlierRejection.RejectOutliers(emptyStack, options);


    }

}