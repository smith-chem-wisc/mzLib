using System;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using System.Linq;
using MathNet;
using System.Collections.Generic; 

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSpectrumSmoothing
    {
        [Test]
        // Standard case. 
        [TestCase(3, 2)]
        // These next two cases are to show what happens when window is as long or
        // longer than signal. They effectively return the average of the whole scan. 
        [TestCase(50, 2)]
        [TestCase(51, 2)] 
        public void TestKZ1D(int windowSize, int iterations)
        {
            // fuzzy square wave. i = 1 from 25 to 34 and 0 at all other points. 
            double[] testArray = new double[50];
            double[] squareWave = new double[10];
            // populate squareWave with 1's.
            for(int i = 0; i < squareWave.Length; i++)
            {
                squareWave[i] = 1D; 
            }
            Buffer.BlockCopy(squareWave, 0, testArray, sizeof(double)*24, sizeof(double) * squareWave.Length);

            // add noise to the square wave:
            // use a set seed to make things reproducible
            Random rnd = new(1551);
            for(int i = 0; i < testArray.Length; i++)
            {
                testArray[i] += rnd.NextDouble()/10; 
            }
            double[] smoothedData = MzSpectrum.KZ1D(testArray, windowSize, iterations);

            // Testing length of output. 
            Assert.AreEqual(50, smoothedData.Length);

            // A better test might be to perform a cross-correlation with the original 
            // square wave and see a high cross-correlation value from indices 24-35. 
        }
        [Test]
        [TestCase(5, 2)]
        public void TestSmoothSpectra(int windowSize, int iterations)
        {
            // construct an mzspectrum object that is just a square wave.
            // square wave construction for the intensity values
            double[] testArray = new double[50];
            double[] squareWave = new double[10];
            // populate squareWave with 1's.
            for (int i = 0; i < squareWave.Length; i++)
            {
                squareWave[i] = 1D;
            }
            Buffer.BlockCopy(squareWave, 0, testArray, sizeof(double) * 24, sizeof(double) * squareWave.Length);
            // m/z values will just be the sequence from 0 to 49.
            double[] mzVals = new double[testArray.Length];
            for(int i = 0; i < mzVals.Length; i++)
            {
                mzVals[i] = (double)i; 
            }

            MzSpectrum testMzSpectrum = new(mzVals, testArray, true);

            double[,] result = testMzSpectrum.SmoothSpectrumKZ(windowSize, iterations);
            // scan average should be zero based on window
            double shouldBeZero = result[1,0];
            // original scan has a zero at 23, but the filtered scan should
            // have a non-zero value. 
            double shouldNotBeZero = result[1,23];
            Assert.AreEqual(50, result.GetLength(1)); 
            Assert.IsTrue(shouldBeZero == 0);
            Assert.IsTrue(shouldNotBeZero > 0); 
        }
    }
    
}
