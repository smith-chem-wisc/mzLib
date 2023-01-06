using System;
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
        BinnedSpectra binnedSpectra = new(this.numSpectra);

        int numSpectra = 3;
        double binSize = 1.0;

        binnedSpectra.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize); 
        Assert.True(binnedSpectra.PixelStacks.Count == 5);
        Assert.That(binnedSpectra.PixelStacks[0].MergedMzValue, Is.EqualTo(0.033333d).Within(0.01));
        Assert.That(binnedSpectra.PixelStacks[2].Intensity, 
            Is.EqualTo(new double[] {12, 13, 30}));
    }

    [Test]
    public void TestConsumeSpectraBinSizes()
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

        BinnedSpectra bs = new(3); 

        bs.ConsumeSpectra(multipleXvals, yVals, 3, 1.0);
        Assert.That(bs.PixelStacks[1].Intensity[0], Is.EqualTo(0d));
        Assert.That(bs.PixelStacks[0].MergedMzValue, Is.EqualTo(1.15));
    }
    
    [Test]
    public void TestConsumeSpectraUnequalArrayLength()
    {
        BinnedSpectra binnedSpectra = new(this.numSpectra);
        int numSpectra = 3;
        double binSize = 1.0;
        binnedSpectra.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        Assert.True(binnedSpectra.PixelStacks.Count == 5);
        Assert.True(binnedSpectra.PixelStacks.All(i => i.Length == 3));
        // 12.5 is the correct answer if the binning and averaging is correct. 
        // Value will be exactly 12.5, so okay to use 12.5 explicitly.  
        // ReSharper disable once CompareOfFloatsByEqualityOperator
        Assert.True(binnedSpectra.PixelStacks[3].Intensity[0] == 12.5d);
    }

    [Test]
    public void TestPerformNormalization()
    {
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        bs.PerformNormalization();
        double ticAfterNorm = bs.PixelStacks.Select(i => i.Intensity[0]).Sum(); 
        Assert.That(ticAfterNorm, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void TestCalculateNoiseEstimate()
    {
        // This is just testing to see if the noise estimates work.
        // The values produces by the NoiseEstimation are tested elsewhere. 
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        bs.PerformNormalization();
        bs.CalculateNoiseEstimates();
        foreach (var estimate in bs.NoiseEstimates)
        {
            Console.WriteLine("Spectra {0}. Noise estimate: {1}", 
                estimate.Key, estimate.Value);
        }
    }

    [Test]
    public void TestScaleEstimate()
    {
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        bs.PerformNormalization();
        bs.CalculateScaleEstimates();
        double[] expected = new[]
        {
            Math.Sqrt(0.024077), 
            Math.Sqrt(0.023317), 
            Math.Sqrt(0.050520)
        };
        Assert.That(bs.ScaleEstimates.Values.ToArray(), 
            Is.EqualTo(expected).Within(0.0001));
    }
    
    [Test]
    public void TestRecalculateTics()
    {
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        double[] expected = new double[] { 59.5, 65d, 150d };
        Assert.That(bs.Tics, Is.EqualTo(expected));
    }
    
    [Test]
    public void TestCalculateWeights()
    {
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        bs.PerformNormalization();
        bs.CalculateNoiseEstimates();
        bs.CalculateScaleEstimates();
        bs.CalculateWeights();
        double[] expectedWeights = new[]
        {
            1539.23913, 
            1636.63560,
            755.37045
        }; 
        Assert.That(bs.Weights.Values.ToArray(), 
            Is.EqualTo(expectedWeights).Within(0.01));
    }

    [Test]
    public void TestProcessPixelStacks()
    {
        SpectralAveragingOptions options = new(); 
        options.SetDefaultValues();
        BinnedSpectra bs = new(this.numSpectra);
        bs.ConsumeSpectra(xArrays, yArrays, numSpectra, binSize);
        bs.RecalculateTics();
        bs.PerformNormalization();
        bs.CalculateNoiseEstimates();
        bs.CalculateScaleEstimates();
        bs.CalculateWeights();
        bs.ProcessPixelStacks(options);
        bs.MergeSpectra();
        bs.GetMergedSpectrum();
    }
}
