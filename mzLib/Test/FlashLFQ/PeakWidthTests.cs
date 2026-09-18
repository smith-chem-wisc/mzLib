using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;
using NUnit.Framework;

namespace Test.FlashLFQ
{
    [TestFixture]
    public class PeakWidthTests
    {
        /// <summary>
        /// Builds a peak from (retentionTime, intensity) points on one charge state. Intensity is
        /// pre-multiplied by the charge because IsotopicEnvelope's constructor divides by it.
        /// </summary>
        private static ChromatographicPeak MakePeak(IEnumerable<(double Rt, double Intensity)> points, int charge = 2)
        {
            var identification = new Identification(
                new SpectraFileInfo("f.mzML", "c", 0, 0, 0),
                "PEPTIDE", "PEPTIDE", 799.36, 10.0, charge, new List<ProteinGroup>());

            var peak = new ChromatographicPeak(identification, new SpectraFileInfo("f.mzML", "c", 0, 0, 0));

            peak.IsotopicEnvelopes = points
                .Select((p, i) => new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(799.36, p.Intensity, i, p.Rt),
                    charge, p.Intensity * charge, 1.0))
                .ToList();

            peak.CalculateIntensityForThisFeature(false);
            return peak;
        }

        /// <summary>
        /// A symmetric triangle rising 0 to 100 over 1 minute and falling back. Half maximum is 50,
        /// crossed exactly halfway up each flank, so FWHM is exactly 1 minute. Chosen because the
        /// answer is arithmetic rather than a fit.
        /// </summary>
        [Test]
        public void TriangleHasExactlyKnownWidth()
        {
            var peak = MakePeak(new[]
            {
                (10.0, 0.0), (10.5, 50.0), (11.0, 100.0), (11.5, 50.0), (12.0, 0.0)
            });

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Measured));
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(width.HalfMaximumStart, Is.EqualTo(10.5).Within(1e-9));
            Assert.That(width.HalfMaximumEnd, Is.EqualTo(11.5).Within(1e-9));
        }

        /// <summary>
        /// A triangle sampled so that no point lands on half maximum: apex 200, flanks at 60, so half
        /// maximum (100) falls between samples. The crossing is 100/140 of the way from the apex down to
        /// the 60 point, i.e. 5/14 min either side, so FWHM is 5/7 min. Deliberately not a multiple of
        /// the 0.5 min sampling interval - without interpolation the answer could only be 0.5 or 1.0.
        /// </summary>
        [Test]
        public void CrossingsAreInterpolatedRatherThanSnappedToSamplePoints()
        {
            var peak = MakePeak(new[]
            {
                (10.0, 0.0), (10.5, 60.0), (11.0, 200.0), (11.5, 60.0), (12.0, 0.0)
            });

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Measured));
            Assert.That(width.HalfMaximumStart, Is.EqualTo(11.0 - (5.0 / 14.0)).Within(1e-9));
            Assert.That(width.HalfMaximumEnd, Is.EqualTo(11.0 + (5.0 / 14.0)).Within(1e-9));
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(5.0 / 7.0).Within(1e-9));

            // The whole point: the answer is not reachable by snapping to sample points.
            Assert.That(width.FullWidthAtHalfMaximum, Is.Not.EqualTo(0.5).Within(1e-6));
            Assert.That(width.FullWidthAtHalfMaximum, Is.Not.EqualTo(1.0).Within(1e-6));
        }

        /// <summary>
        /// A gaussian, where FWHM is 2*sqrt(2*ln2)*sigma analytically. Sampled densely enough that
        /// linear interpolation across one interval is a good approximation. This is the check that the
        /// measurement means what the name says on a realistic shape, not just on a triangle.
        /// </summary>
        [Test]
        public void GaussianMatchesTheAnalyticWidth()
        {
            const double sigma = 0.25;
            const double centre = 20.0;
            double expected = 2.0 * Math.Sqrt(2.0 * Math.Log(2.0)) * sigma;

            var points = new List<(double, double)>();
            for (double rt = centre - 1.5; rt <= centre + 1.5 + 1e-9; rt += 0.02)
            {
                points.Add((rt, 1000.0 * Math.Exp(-Math.Pow(rt - centre, 2) / (2 * sigma * sigma))));
            }

            PeakWidth width = MakePeak(points).PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Measured));
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(expected).Within(0.005),
                $"expected {expected:F4} min for sigma {sigma}");
        }

        /// <summary>
        /// A trace that is still above half maximum where it ends is censored, not measured. Reporting a
        /// number here would report the truncation rather than the peak.
        /// </summary>
        [Test]
        public void TraceEndingAboveHalfMaximumIsCensored()
        {
            var peak = MakePeak(new[]
            {
                (10.0, 60.0), (10.5, 80.0), (11.0, 100.0), (11.5, 90.0), (12.0, 70.0)
            });

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Censored));
            Assert.That(width.FullWidthAtHalfMaximum, Is.NaN);
        }

        /// <summary>
        /// Censoring on one side only is still censoring: half of a width is not a width.
        /// </summary>
        [Test]
        public void OneSidedCensoringIsStillCensored()
        {
            var peak = MakePeak(new[]
            {
                (10.0, 60.0), (10.5, 80.0), (11.0, 100.0), (11.5, 50.0), (12.0, 10.0)
            });

            Assert.That(peak.PeakWidth.Status, Is.EqualTo(PeakWidthStatus.Censored));
        }

        /// <summary>
        /// Under the point floor nothing is reported. Matches CutPeak, which declines to look for
        /// valleys below five envelopes.
        /// </summary>
        [Test]
        public void TooFewPointsIsNotMeasured()
        {
            var peak = MakePeak(new[] { (10.0, 0.0), (10.5, 50.0), (11.0, 100.0), (11.5, 0.0) });

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.TooFewPoints));
            Assert.That(width.TimePointsOnApexCharge, Is.EqualTo(4));
        }

        /// <summary>
        /// The measurement must use the apex charge state alone. Here +2 is a clean 1-minute triangle
        /// and +3 is interleaved noise at the same retention times. Walking every envelope in retention
        /// time order would describe a sawtooth; the answer must be unchanged by the +3 points.
        /// </summary>
        [Test]
        public void InterleavedSecondChargeStateDoesNotChangeTheWidth()
        {
            var apexChargePoints = new[]
            {
                (10.0, 0.0), (10.5, 50.0), (11.0, 100.0), (11.5, 50.0), (12.0, 0.0)
            };

            ChromatographicPeak peak = MakePeak(apexChargePoints, charge: 2);
            double widthWithoutInterference = peak.PeakWidth.FullWidthAtHalfMaximum;

            // Same retention times, a different charge, intensities that would wreck a naive walk.
            var interference = apexChargePoints
                .Select((p, i) => new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(533.24, 5.0, i, p.Item1), 3, 5.0 * 3, 1.0));

            peak.IsotopicEnvelopes.AddRange(interference);
            peak.CalculateIntensityForThisFeature(false);

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Measured));
            Assert.That(width.TimePointsOnApexCharge, Is.EqualTo(5), "the +3 envelopes were not filtered out");
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(widthWithoutInterference).Within(1e-9));
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(1.0).Within(1e-9));
        }

        /// <summary>
        /// A shoulder that dips below half maximum and comes back must not be swallowed. Taking the
        /// first and last crossing of the whole trace would report 3.0 minutes here; the nearest
        /// crossings either side of the apex give 1.0.
        /// </summary>
        [Test]
        public void ShoulderIsNotIncludedInTheWidth()
        {
            var peak = MakePeak(new[]
            {
                (10.0, 0.0), (10.5, 50.0), (11.0, 100.0), (11.5, 50.0), (12.0, 0.0),
                (12.5, 40.0), (13.0, 60.0), (13.5, 40.0), (14.0, 0.0)
            });

            PeakWidth width = peak.PeakWidth;

            Assert.That(width.Status, Is.EqualTo(PeakWidthStatus.Measured));
            Assert.That(width.FullWidthAtHalfMaximum, Is.EqualTo(1.0).Within(1e-9));
        }

        /// <summary>
        /// Every row must carry as many fields as the header, for both the base peak and the MBR subclass
        /// that duplicates the whole row rather than extending it. Adding the two width columns to
        /// ChromatographicPeak alone shifted every MBR row two fields out of step with the header, which
        /// is invisible in a spot check of one row and corrupts the whole file.
        /// </summary>
        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public void EveryRowHasAsManyFieldsAsTheHeader(bool withApex)
        {
            int headerFields = ChromatographicPeak.TabSeparatedHeader.Split('\t').Length;

            var file = new SpectraFileInfo("f.mzML", "c", 0, 0, 0);
            var identification = new Identification(
                file, "PEPTIDE", "PEPTIDE", 799.36, 10.0, 2, new List<ProteinGroup>());

            var peaks = new ChromatographicPeak[]
            {
                new ChromatographicPeak(identification, file),
                new MbrChromatographicPeak(identification, file, 1, false)
            };

            foreach (ChromatographicPeak peak in peaks)
            {
                peak.IsotopicEnvelopes = new[]
                    {
                        (10.0, 0.0), (10.5, 50.0), (11.0, 100.0), (11.5, 50.0), (12.0, 0.0)
                    }
                    .Select((p, i) => new IsotopicEnvelope(
                        new IndexedMassSpectralPeak(799.36, p.Item2, i, p.Item1), 2, p.Item2 * 2, 1.0))
                    .ToList();

                // Apex is only set by CalculateIntensityForThisFeature, so skipping it exercises the
                // no-apex branch, which writes a different set of literals.
                if (withApex)
                {
                    peak.CalculateIntensityForThisFeature(false);
                }

                Assert.That(peak.ToString().Split('\t').Length, Is.EqualTo(headerFields),
                    $"{peak.GetType().Name} wrote a row of a different width than the header " +
                    $"(apex {(withApex ? "set" : "null")})");
            }
        }

        /// <summary>
        /// A peak with no apex reports NoApex rather than throwing.
        /// </summary>
        [Test]
        public void EmptyPeakReportsNoApex()
        {
            var identification = new Identification(
                new SpectraFileInfo("f.mzML", "c", 0, 0, 0),
                "PEPTIDE", "PEPTIDE", 799.36, 10.0, 2, new List<ProteinGroup>());
            var peak = new ChromatographicPeak(identification, new SpectraFileInfo("f.mzML", "c", 0, 0, 0));
            peak.CalculateIntensityForThisFeature(false);

            Assert.That(peak.PeakWidth.Status, Is.EqualTo(PeakWidthStatus.NoApex));
        }
    }
}
