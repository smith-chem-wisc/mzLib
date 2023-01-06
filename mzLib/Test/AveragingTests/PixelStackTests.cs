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
    private PixelStack pixelStack; 
    [OneTimeSetUp]
    public void OneTimeSetup()
    {
        xArray = new[] { 1.0, 1.1, 0.9, 0.99 };
        yArray = new[] { 10.0, 11.0, 20.0, 50.0 };
        pixelStack = new PixelStack(xArray, yArray);
        pixelStack.Reject(pixelStack.Length - 1);
    }

    [Test]
    public void TestGetIntensities()
    {
        Assert.That(pixelStack.Intensity, Is.EqualTo(yArray)); 
    }

    [Test]
    public void TestGetMzValues()
    {
        Assert.That(pixelStack.Mzs, Is.EqualTo(xArray));
    }

    [Test]
    public void TestNonRejectedLength()
    {
        Assert.True(pixelStack.NonRejectedLength == pixelStack.Length - 1);
    }

    [Test]
    public void TestGetUnrejectedValues()
    {
        Assert.That(pixelStack.UnrejectedIntensities, Is.EqualTo(yArray[..^1]));
        Assert.That(pixelStack.UnrejectedMzs, Is.EqualTo(xArray[..^1]));
    }

    [Test]
    public void TestIsIndexRejected()
    {
        Assert.That(pixelStack.IsIndexRejected(pixelStack.Length - 1));
    }

    [Test]
    public void TestGetIntensityAtIndex()
    {
        Assert.AreEqual(pixelStack.GetIntensityAtIndex(0), yArray[0]); 
    }

    [Test]
    public void TestGetMzAtIndex()
    {
        Assert.AreEqual(pixelStack.GetMzAtIndex(0), xArray[0]);
    }

    [Test]
    public void TestMzAverage()
    {
        Assert.True(pixelStack.MergedMzValue == xArray[..^1].Average()); 
    }

    [Test]
    public void TestAverage()
    {
        Dictionary<int, double> weights = new(); 
        weights.Add(0, 0.2);
        weights.Add(1, 0.4);
        weights.Add(2, 0.8);
        weights.Add(3, 1.0);
        pixelStack.Average(weights);
        Assert.That(pixelStack.MergedIntensityValue,
            Is.EqualTo(15.9999).Within(0.01));
    }

    [Test]
    [TestCase(RejectionType.WinsorizedSigmaClipping)]
    [TestCase(RejectionType.SigmaClipping)]
    [TestCase(RejectionType.AveragedSigmaClipping)]
    [TestCase(RejectionType.MinMaxClipping)]
    [TestCase(RejectionType.PercentileClipping)]
    [TestCase(RejectionType.BelowThresholdRejection)]
    [TestCase(RejectionType.NoRejection)]
    public void TestEmptyPixelStack(RejectionType rejection)
    {
        SpectralAveragingOptions options = new(); 
        options.SetDefaultValues();
        options.RejectionType = rejection; 
        List<double> emptyXarray = new();
        List<double> emptyYarray = new(); 
        PixelStack emptyStack = new(emptyXarray, emptyYarray); 
        emptyStack.PerformRejection(options);


    }

}