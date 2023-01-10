using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests;

[ExcludeFromCodeCoverage]
public class TestBinnedSpectra
{
    private double[][] xArrays;
    private double[][] yArrays;
    private double[] tics;
    private int numSpectra;
    private double binSize; 

    [OneTimeSetUp]
    public void OneTimeSetUp()
    {
        xArrays = new[]
        {
            new double[] { 0, 1, 2, 3, 3.49, 4 },
            new double[] { 0, 1, 2, 3, 4 },
            new double[] { 0.1, 1.1, 2.1, 3.1, 4.1}
        };
        yArrays = new[]
        {
            new double[] { 10, 11, 12, 12, 13, 14 },
            new double[] { 11, 12, 13, 14, 15 },
            new double[] { 20, 25, 30, 35, 40 }
        };
        tics = new double[3];
        tics[0] = yArrays[0].Sum();
        tics[1] = yArrays[1].Sum();
        tics[2] = yArrays[2].Sum();
        numSpectra = 3;
        binSize = 1.0; 

        SpectralAveragingOptions options = new SpectralAveragingOptions();
    }

    [Test]
    public void TestConsumeSpectra()
    {
        BinnedSpectra binnedSpectra = new(xArrays, yArrays, binSize); 
        Assert.True(binnedSpectra.PeakBins.Count == 5);
        Assert.That(binnedSpectra.PeakBins[0].UnrejectedMzAverage, Is.EqualTo(0.033333d).Within(0.01));
        Assert.That(binnedSpectra.PeakBins[2].Intensities, 
            Is.EqualTo(new double[] {12, 13, 30}));
    }

    [Test]
    [TestCase(0, 9, new[] {1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1}, new[] { 10.0, 20.0, 10.0, 10.0, 20.0, 10.0, 10.0, 20.0, 10.0 })]
    [TestCase(1, 3, new[] {3.0, 3.0, 3.0}, new[] {30.0, 30.0, 30.0} )]
    [TestCase(2, 3, new[] { 4.0, 4.0, 4.0 }, new[] { 20.0, 20.0, 20.0 })]
    [TestCase(3, 3, new[] { 5.0, 5.0, 5.0 }, new[] { 10.0, 10.0, 10.0 })]
    [TestCase(4,3,  new[] { 6.0, 6.0, 6.0 }, new[] { 10.0, 10.0, 10.0 })]
    public void BinnedSpectraConstructorAndConsumeSpectra(int stackIndex, int stackCount, double[] expectedMz, double[] expectedIntensity)
    {
        double[][] multipleXvals = new double[3][];
        double[][] xValsInBin = new double[3][];
        double[][] yVals = new double[3][];
        
        double[] yAxis = new double[] { 10, 10, 20, 30, 20, 10, 10 };
        double[] xAxisSkip = new double[] { 1, 1.1, 1.2, 3.0, 4.0, 5.0, 6.0 };

        multipleXvals[0] = xAxisSkip; 
        multipleXvals[1] = xAxisSkip; 
        multipleXvals[2] = xAxisSkip;
        yVals[0] = yAxis;
        yVals[1] = yAxis;
        yVals[2] = yAxis;

        BinnedSpectra bs = new(multipleXvals, yVals, 1.0);
        Assert.That(bs.PeakBins.Count == 5);


        var stack = bs.PeakBins[stackIndex];
        Assert.That(stack.Intensities.Count() == stackCount);
        Assert.That(stack.Mzs.SequenceEqual(expectedMz));
        Assert.That(stack.UnrejectedMzs.SequenceEqual(expectedMz));
        Assert.That(stack.Intensities.SequenceEqual(expectedIntensity));
        Assert.That(stack.UnrejectedIntensities.SequenceEqual(expectedIntensity));
        Assert.That(Math.Abs(stack.UnrejectedMzAverage - stack.Mzs.Average()) < 0.001);
    }
    
    [Test]
    [TestCase(0, 3, new[] {0, 0, 0.1 }, new[] {10.0, 11.0, 20.0})]
    [TestCase(1, 3, new[] {1, 1, 1.1 }, new[] {11.0, 12.0, 25.0})]
    [TestCase(2, 3, new[] {2, 2, 2.1 }, new[] {12.0, 13.0, 30.0})]
    [TestCase(3, 4, new[] {3.245, 3.245, 3, 3.1 }, new[] {13.0, 12.0, 14.0, 35.0})]
    [TestCase(4, 3, new[] {4, 4, 4.1}, new[] {14.0, 15.0, 40.0})]
    public void TestConsumeSpectraUnequalArrayLength(int stackIndex, int stackCount, double[] expectedMz, double[] expectedIntensity)
    {
        BinnedSpectra binnedSpectra = new(xArrays, yArrays, binSize);
        var stack = binnedSpectra.PeakBins[stackIndex];
        Assert.That(stack.Intensities.Count() == stackCount);
        Assert.That(stack.Mzs.SequenceEqual(expectedMz));
        Assert.That(stack.UnrejectedMzs.SequenceEqual(expectedMz));
        Assert.That(stack.Intensities.SequenceEqual(expectedIntensity));
        Assert.That(stack.UnrejectedIntensities.SequenceEqual(expectedIntensity));
        Assert.That(Math.Abs(stack.UnrejectedMzAverage - stack.Mzs.Average()) < 0.001);
    }

    [Test]
    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    public void TestPerformNormalization(int spectraId)
    {
        BinnedSpectra bs = new(xArrays, yArrays, binSize);
        bs.PerformNormalization();
        List<Peak> peaksFromSpectra = new();

        bs.PeakBins.ForEach(p => peaksFromSpectra.AddRange(p.Peaks.Where(m => m.SpectraId == spectraId)));
        double ticAfterNorm = peaksFromSpectra.Select(p => p.Intensity).Sum();
        Assert.That(ticAfterNorm, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void TestCalculateTics()
    {
        BinnedSpectra bs = new(xArrays, yArrays, binSize);
        double[] expected = new double[] { 72d, 65d, 150d };
        Assert.That(bs.Tics, Is.EqualTo(expected));
    }
    
    [Test]
    public void TestCalculateWeights()
    {
        BinnedSpectra bs = new(xArrays, yArrays, binSize);
        bs.PerformNormalization();
        SpectralWeighting.CalculateSpectraWeights(bs, new SpectralAveragingOptions() {SpectralWeightingType = SpectraWeightingType.MrsNoiseEstimation});
        double[] expectedWeights = new[]
        {
            //1539.23913, 
            //1636.63560,
            //755.37045
            2073.60000,
            3785.99353127,
            889.4099838
        }; 
        Assert.That(bs.Weights.Values.ToArray(), 
            Is.EqualTo(expectedWeights).Within(0.01));
    }


}
