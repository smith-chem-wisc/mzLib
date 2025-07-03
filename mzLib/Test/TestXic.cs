using System;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using FlashLFQ;
using MzLibUtil;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MassSpectrometry;
using Microsoft.ML.Transforms;
using MathNet.Numerics.Distributions;
using System.Collections;

namespace Test
{
    public class TestXic
    {
        public static MsDataScan[] FakeScans { get; set; }
        public static IsotopicDistribution Dist { get; set; }
        public static List<IIndexedPeak> PeakList { get; set; } 

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            FakeScans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            Dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            // Create mzSpectra where two peaks appear very close together
            for (int s = 0; s < FakeScans.Length; s++)
            {
                double[] mz = Dist.Masses.SelectMany(v => new List<double> { v.ToMz(1), (v + 0.0001).ToMz(1) }).ToArray();
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * intensity * intensityMultipliers[s], v * intensity }).ToArray();

                // add the scan
                FakeScans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Example XIC
            var peakIndexEngine = PeakIndexingEngine.InitializeIndexingEngine(FakeScans);
            PeakList = peakIndexEngine.GetXic(Dist.Masses.First().ToMz(1), zeroBasedStartIndex: 4, new PpmTolerance(20), 1);
        }

        [Test]
        public static void TestXicConstruction()
        {
            var xic = new ExtractedIonChromatogram(PeakList);

            //Test XIC properties
            Assert.That(xic.Peaks.Count, Is.EqualTo(10));
            Assert.That(xic.ApexRT, Is.EqualTo(1.6));
            Assert.That(xic.ApexScanIndex, Is.EqualTo(6));
            Assert.That(xic.StartRT, Is.EqualTo(1.0));
            Assert.That(xic.EndRT, Is.EqualTo(1.9));
            Assert.That(xic.AveragedM, Is.EqualTo(Dist.Masses.First().ToMz(1)).Within(0.0000001));

            //Test normalized peak intensities
            xic.SetNormalizedPeakIntensities();
            Assert.That(xic.NormalizedPeakIntensities.Length, Is.EqualTo(10));
            Assert.That(xic.NormalizedPeakIntensities.Sum(), Is.EqualTo(100).Within(0.0001));
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
            Assert.That(xic.XYData.Max(xy => xy.Item1), Is.EqualTo(2).Within(0.0000001));
        }

        [Test]
        public static void TestExceptionHandling()
        {
            var cubicSpline = new XicCubicSpline(0.05);
            var linearSpline = new XicLinearSpline(0.05);
            var rtArray = new double[] { 1.0, 1.1, 1.2 };
            var intensityArray = new double[] { 100, 200, 300 };
            var ex = Assert.Throws<MzLibException>(() => cubicSpline.GetXicSplineData(rtArray, intensityArray, 1.0, 1.2));
            Assert.That(ex.Message, Is.EqualTo("Input arrays must contain at least 5 points."));
            var intensityArray2 = new double[] { 100, 200, 300, 400, 500 };
            var ex2 = Assert.Throws<MzLibException>(() => cubicSpline.GetXicSplineData(rtArray, intensityArray2, 1.0, 1.2));
            Assert.That(ex2.Message, Is.EqualTo("Input arrays must have the same length."));
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

            //if there is only one peak, it should not be cut
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
    }
}
