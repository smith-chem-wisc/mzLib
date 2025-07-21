using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

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
            Assert.That(xic.ApexRT, Is.EqualTo(1.6f));
            Assert.That(xic.ApexScanIndex, Is.EqualTo(6));
            Assert.That(xic.StartRT, Is.EqualTo(1.0f));
            Assert.That(xic.EndRT, Is.EqualTo(1.9f));
            var mass = Dist.Masses.First().ToMz(1);
            Assert.That(xic.AveragedM, Is.EqualTo(Dist.Masses.First().ToMz(1)).Within(0.0001));

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
            Assert.That(xic.XYData.Length, Is.EqualTo(18)); // Because the last time point will be stored as 18.999999 (origin 19) while convert the float to double. 
            linearSpline.SetXicSplineXYData(xic);                            // Then we will lose one numPoint (19 to 18) in the XYData. 
            Assert.That(xic.XYData.Length, Is.EqualTo(18));
            //in scan cycle
            cubicSpline.SetXicSplineXYData(xic, cycle: true);
            Assert.That(xic.XYData.Length, Is.EqualTo(181));

            //Test adding peaks
            var cubicSpline2 = new XicCubicSpline(0.05, 1, 0.1);
            cubicSpline2.SetXicSplineXYData(xic);
            Assert.That(xic.XYData.Min(xy => xy.Item1), Is.EqualTo(0.9).Within(0.0000001));
            Assert.That(xic.XYData.Max(xy => xy.Item1), Is.EqualTo(1.95).Within(0.0000001)); // Because we lost one numPoint, then last point will be 1.95 instead of original value 2.0.

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
        public static void TestExceptionHandling()
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
        }
    }
}
