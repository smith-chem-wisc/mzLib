﻿using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Interpolation;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MassSpectrometry;
using Microsoft.ML.Transforms;
using MathNet.Numerics.Distributions;
using System.Collections;
using Proteomics.AminoAcidPolymer;

namespace Test
{
    public class TestXic
    {
        public static MsDataScan[] testScans { get; set; }
        public static IsotopicDistribution Dist { get; set; }
        public static List<IIndexedPeak> PeakList { get; set; } 

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            testScans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            Dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < testScans.Length; s++)
            {
                double[] mz = Dist.Masses.SelectMany(v => new List<double> { v.ToMz(1), (v + 0.0001).ToMz(1) }).ToArray();
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * intensity * intensityMultipliers[s], v * intensity }).ToArray();

                // add the scan
                testScans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Example XIC
            var peakIndexEngine = PeakIndexingEngine.InitializeIndexingEngine(testScans);
            PeakList = peakIndexEngine.GetXicByScanIndex(Dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
        }

        [Test]
        public static void TestXicConstruction()
        {
            var xic = new ExtractedIonChromatogram(PeakList);

            //Test XIC properties
            Assert.That(xic.Peaks.Count, Is.EqualTo(10));
            Assert.That(xic.ApexRT, Is.EqualTo(1.6f));
            Assert.That(xic.ApexScanIndex, Is.EqualTo(6));
            Assert.That(xic.StartRT, Is.EqualTo(1.0f));
            Assert.That(xic.EndRT, Is.EqualTo(1.9f));
            var mass = Dist.Masses.First().ToMz(1);
            Assert.That(xic.AveragedMassOrMz, Is.EqualTo(Dist.Masses.First().ToMz(1)).Within(0.0001));
            Assert.That(xic.ApexPeak.Intensity, Is.EqualTo(Dist.Intensities.First() * 1e6 * 10).Within(1));

            //Test normalized peak intensities
            xic.SetNormalizedPeakIntensities();
            Assert.That(xic.NormalizedPeakIntensities.Length, Is.EqualTo(10));
            Assert.That(xic.NormalizedPeakIntensities.Sum(), Is.EqualTo(100f).Within(0.0001));
        }

        [Test]
        public static void TestXicSpline()
        {
            var xic = new ExtractedIonChromatogram(PeakList);

            //Test XIC spline
            //in retention time
            var cubicSpline = new XicCubicSpline(0.05);
            var linearSpline = new XicLinearSpline(0.05);
            cubicSpline.SetXicSplineXYData(xic);
            Assert.That(xic.XYData.Length, Is.EqualTo(19)); 
            linearSpline.SetXicSplineXYData(xic);  
            Assert.That(xic.XYData.Length, Is.EqualTo(19));
            //in scan cycle
            cubicSpline.SetXicSplineXYData(xic, cycle: true);
            Assert.That(xic.XYData.Length, Is.EqualTo(181));

            //Test adding peaks
            var cubicSpline2 = new XicCubicSpline(0.05, 1, 0.1);
            cubicSpline2.SetXicSplineXYData(xic);
            Assert.That(xic.XYData.Min(xy => xy.Item1), Is.EqualTo(0.9).Within(0.0000001));
            Assert.That(xic.XYData.Max(xy => xy.Item1), Is.EqualTo(2.0).Within(0.0000001)); 

            //ensure that add peaks works correctly
            var peakList1 = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 3, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList1.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i + 5, mz: 500.0));
            }
            //The xic contains three original peaks, and we want to add two more peaks at the begining and the end. This setting of spline only adds peaks, no spline.
            var xic1 = new ExtractedIonChromatogram(peakList1);
            int numberOfPeaksToAdd = 2;
            var linearSpline2 = new XicLinearSpline(1, numberOfPeaksToAdd, 1);
            linearSpline2.SetXicSplineXYData(xic1, cycle: true);
            //the xic length after adding the peaks should be 3 + 2*2
            //the first two peaks and the last two peaks should have intensity 0
            Assert.That(xic1.XYData.Length, Is.EqualTo(7));
            for (int i = 0; i < numberOfPeaksToAdd; i++)
            {
                Assert.That(xic1.XYData[i].Item2, Is.EqualTo(0)); 
            }
            for (int i = xic1.XYData.Length - 1; i > xic1.Peaks.Count + numberOfPeaksToAdd - 1; i--)
            {
                Assert.That(xic1.XYData[i].Item2, Is.EqualTo(0));
            }
            //the orginal peaks should be present in the XYData after the 0 peaks are added
            foreach (var peak in xic1.Peaks)
            {
                Assert.That(xic1.XYData.First(xy => xy.Item1 == peak.ZeroBasedScanIndex).Item2 == peak.Intensity, Is.True);
            }
        }

        [Test]
        public static void TestXicSplineExceptionHandling()
        {
            var cubicSpline = new XicCubicSpline(0.05);
            var linearSpline = new XicLinearSpline(0.05);
            var rtArray = new float[] { 1.0f, 1.1f, 1.2f };
            var intensityArray = new float[] { 100, 200, 300 };
            var ex = Assert.Throws<MzLibException>(() => cubicSpline.GetXicSplineData(rtArray, intensityArray, 1.0, 1.2));
            Assert.That(ex.Message, Is.EqualTo("Input arrays must contain at least 5 points."));
            var intensityArray2 = new float[] { 100, 200, 300, 400, 500 };
            var ex2 = Assert.Throws<MzLibException>(() => cubicSpline.GetXicSplineData(rtArray, intensityArray2, 1.0, 1.2));
            Assert.That(ex2.Message, Is.EqualTo("Input arrays must have the same length."));

            var Bspline = new Bspline(4, 10);
            var ex3 = Assert.Throws<MzLibException>(() => Bspline.GetXicSplineData(rtArray, intensityArray));
            Assert.That(ex3.Message, Is.EqualTo("The number of points in the input array must be greater than the degree of the Bspline."));
        }

        [Test]
        public static void TestGetAllXics()
        {
            //Test case using FakeScans
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testScans);
            var xics = indexingEngine.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(xics.Count, Is.EqualTo(20));
            foreach(var xic in xics)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(10));
            }

            //Test with massIndexingEngine
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var allMasses = Deconvoluter.Deconvolute(testScans[0].MassSpectrum, deconParameters);
            var massIndexingEngine = MassIndexingEngine.InitializeMassIndexingEngine(testScans, deconParameters);
            var massXics = massIndexingEngine.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(massXics.Count, Is.EqualTo(2));
            foreach (var mass in allMasses)
            {
                Assert.That(massXics.Any(x => x.Peaks.Select(p => (IndexedMass)p).First().IsotopicEnvelope.MonoisotopicMass == mass.MonoisotopicMass));
            }
            foreach (var xic in massXics)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(10));
            }

            //Test if there are three missed scans in the middle
            var fakeScans2 = (MsDataScan[])testScans.Clone();
            var missedIndices = new List<int> { 3, 4, 5 }; //the ten scnas would be 1,1,1,0,0,0,1,1,1,1
            foreach(var s in missedIndices)
            {
                // replace the scan at index s with an empty scan
                fakeScans2[s] = new MsDataScan(massSpectrum: new MzSpectrum(new double[] {  }, new double[] {  }, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
            var indexingEngine2 = PeakIndexingEngine.InitializeIndexingEngine(fakeScans2);
            var xics2 = indexingEngine2.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(xics2.Count, Is.EqualTo(40)); //the first three scans and the last four scans will each contain two XICs
            Assert.That(xics2.SequenceEqual(xics2.OrderByDescending(x => x.ApexPeak.Intensity)));//make sure the XICs are ordered by descending apex intensity
            //Test with massIndexingEngine
            var massIndexingEngine2 = MassIndexingEngine.InitializeMassIndexingEngine(fakeScans2, deconParameters);
            var massXics2 = massIndexingEngine2.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(massXics2.Count, Is.EqualTo(4));

            var fakeScans3 = (MsDataScan[])testScans.Clone();
            var missedIndices2 = new List<int> { 2, 3, 4 }; //the ten scnas would be 1,1,0,0,0,1,1,1,1,1
            foreach (var s in missedIndices2)
            {
                fakeScans3[s] = new MsDataScan(massSpectrum: new MzSpectrum(new double[] {  }, new double[] {  }, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
            var indexingEngine3 = PeakIndexingEngine.InitializeIndexingEngine(fakeScans3);
            var xics3 = indexingEngine3.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(xics3.Count, Is.EqualTo(20)); // Because the minumum number of peaks required is set to 3, the first two scans do not contain any XICs with only two consecutive peaks.
            foreach (var xic in xics3)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(5));
            }
            //Test with massIndexingEngine
            var massIndexingEngine3 = MassIndexingEngine.InitializeMassIndexingEngine(fakeScans3, deconParameters);
            var massXics3 = massIndexingEngine3.GetAllXics(new PpmTolerance(20), 2, 2, 3);
            Assert.That(massXics3.Count, Is.EqualTo(2));
            foreach (var xic in massXics3)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(5));
            }
        }   

        [Test]
        public static void TestCutPeak()
        {
            var normal = new Normal(10, 1);
            var RTs = Enumerable.Range(0, 10).Select(i => 7.5 + i * 0.5).ToArray();
            var intensities = RTs.Select(r => normal.Density(r)).ToArray();

            //generate two sets of peaks and put them together so the integrated xic shows double apex
            var peak1 = new List<IIndexedPeak>();
            for (int i = 0; i < RTs.Length; i++)
            {
                peak1.Add(new IndexedMassSpectralPeak(intensity: intensities[i] * 10, retentionTime: RTs[i], zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic1 = new ExtractedIonChromatogram(peak1);
            var peak2 = new List<IIndexedPeak>();
            for (int i = 0; i < RTs.Length - 1; i++)
            {
                peak2.Add(new IndexedMassSpectralPeak(intensity: intensities[i], retentionTime: RTs[RTs.Length - 1] + (i + 1) * 0.5, zeroBasedScanIndex: i + RTs.Length, mz: 500.0));
            }
            var xic = new ExtractedIonChromatogram(peak1.Concat(peak2).OrderBy(p => p.RetentionTime).ToList());
            Assert.That(xic.Peaks.Count, Is.EqualTo(peak1.Count + peak2.Count));
            xic.CutPeak();
            Assert.That(xic.Peaks.Count, Is.EqualTo(peak1.Count));
            Assert.That(xic.ApexRT, Is.EqualTo(xic1.ApexRT));
            Assert.That(xic.StartRT, Is.EqualTo(xic1.StartRT));
            Assert.That(xic.EndRT, Is.EqualTo(xic1.EndRT));

            //if there is only one apex, it should not be cut
            xic1.CutPeak();
            Assert.That(xic1.Peaks.Count, Is.EqualTo(peak1.Count));

            //if the number of peaks is smaller than 5, it should not be cut
            var xic2 = new ExtractedIonChromatogram(peak1.Take(3).ToList());
            xic2.CutPeak();
            Assert.That(xic2.Peaks.Count, Is.EqualTo(3));

            //if the intensity difference does not exceed the discrimination factor, it should not be cut
            var peak3 = new List<IIndexedPeak>();
            var intensityIncrement = intensities[intensities.Length - 1] * 0.01;
            for (int i = 0; i < 4; i++)
            {
                peak3.Add(new IndexedMassSpectralPeak(intensity: intensities[intensities.Length - 1] + intensityIncrement * i, retentionTime: RTs[RTs.Length - 1] + (i + 1) * 0.5, zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic3 = new ExtractedIonChromatogram(peak1.Concat(peak3).OrderBy(p => p.RetentionTime).ToList());
            Assert.That(xic3.Peaks.Count, Is.EqualTo(peak1.Count + peak3.Count));
            xic3.CutPeak();
            Assert.That(xic3.Peaks.Count, Is.EqualTo(peak1.Count + peak3.Count));
        }

        [Test]
        public static void TestGetAllXicsWithCutPeak()
        {
            var testScans = new MsDataScan[10];
            //make all peak traces to have two apex at the third and the eighth scan
            double[] intensityMultipliers = { 1, 2, 6, 2, 1, 1.1, 2, 8, 2, 1.1 };
            var intensity = 1e6;

            for (int s = 0; s < testScans.Length; s++)
            {
                double[] mz = Dist.Masses.SelectMany(v => new List<double> { v.ToMz(1) }).ToArray();
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * intensity * intensityMultipliers[s] }).ToArray();
                testScans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(testScans);
            var xics = indexingEngine.GetAllXics(new PpmTolerance(5), 2, 3, 3);
            // 10 XICs should be found without peak cutting, each containing 10 peaks
            Assert.That(xics.Count(), Is.EqualTo(10)); 
            Assert.That(xics.All(x => x.Peaks.Count == 10), Is.True);
            //The XICs should have an apex at the eighth scan
            Assert.That(xics.All(x => x.ApexScanIndex == 7), Is.True);

            var xicsWithCutPeak = indexingEngine.GetAllXics(new PpmTolerance(5), 2, 3, 3, 0.6);
            // since each original XIC contain two apex, cutting peaks should double the number of XICs found and each XIC should contain 5 peaks now
            Assert.That(xicsWithCutPeak.Count(), Is.EqualTo(20));
            Assert.That(xicsWithCutPeak.All(x => x.Peaks.Count == 5), Is.True);
            // The XICs should have an apex at either the third or the eighth scan
            Assert.That(xicsWithCutPeak.All(x => x.ApexScanIndex == 2 || x.ApexScanIndex == 7), Is.True);
        }

        [Test]
        public static void TestMassXicExceptionHandling()
        {
            var peakIndexEngine = PeakIndexingEngine.InitializeIndexingEngine(testScans);
            var ex = Assert.Throws<MzLibException>(() => peakIndexEngine.GetXicByScanIndex(Dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1, 10, 1));
            Assert.That(ex.Message, Is.EqualTo("Error: Attempted to access a peak using a charge parameter, but the peaks do not have charge information available."));
        }

        [Test]
        public static void TestXicOnePeak()
        {
            var peaks = new List<IIndexedPeak> { new IndexedMassSpectralPeak(intensity: 10, retentionTime: 1, zeroBasedScanIndex: 0, mz: 500.0) };
            var xic = new ExtractedIonChromatogram(peaks);
            Assert.That(xic.Peaks.Count, Is.EqualTo(1));
            Assert.That(xic.ApexRT, Is.EqualTo(peaks.First().RetentionTime));
            Assert.That(xic.ApexScanIndex, Is.EqualTo(peaks.First().ZeroBasedScanIndex));
            Assert.That(xic.StartRT, Is.EqualTo(peaks.First().RetentionTime));
            Assert.That(xic.EndRT, Is.EqualTo(peaks.First().RetentionTime));
            Assert.That(xic.AveragedMassOrMz, Is.EqualTo(peaks.First().M).Within(0.0001));
            Assert.That(xic.ApexPeak.Intensity, Is.EqualTo(peaks.First().Intensity));
        }

        [Test]
        public static void TestBspline()
        {
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var peakList = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10.0, zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic = new ExtractedIonChromatogram(peakList);
            var bspline1 = new Bspline(1, 50);
            bspline1.SetXicSplineXYData(xic);
            //When the smoothing degree is 1, it is equivalent to linear spline, so the apex of the interpolated points should be the same as the original xic
            Assert.That(xic.XYData.MaxBy(p => p.Item2).Item1, Is.EqualTo(xic.ApexRT).Within(0.0001));
            Assert.That(xic.XYData.Max(p => p.Item2), Is.EqualTo(xic.Peaks.Max(p => p.Intensity)).Within(0.0001));

            var bspline2 = new Bspline(2, 150);
            bspline2.SetXicSplineXYData(xic);
            //When the smoothing degree is 2 the apex should still be around the original apex RT
            Assert.That(xic.XYData.MaxBy(p => p.Item2).Item1, Is.EqualTo(xic.ApexRT).Within(0.001));

            //The start and end RT should always stay the same
            Assert.That(xic.XYData.First().Item1, Is.EqualTo(xic.StartRT).Within(0.0001));
            Assert.That(xic.XYData.Last().Item1, Is.EqualTo(xic.EndRT).Within(0.0001));
        }
    }
}
