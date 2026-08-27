using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

namespace Test.FlashLFQ
{
    [TestFixture]
    public class PeakWidthSummaryTests
    {
        /// <summary>
        /// Widths are derived from IIndexedPeak.RetentionTime, which is a float, so a width built from
        /// retention times that are not exactly representable carries about 1e-6 min of representation
        /// error. Irrelevant against peak widths of tenths of a minute, but it rules out exact equality.
        /// </summary>
        private const double FloatRetentionTimeTolerance = 1e-5;

        /// <summary>
        /// A symmetric triangle whose half-maximum crossings sit exactly on sample points, so its width
        /// is the requested value with no interpolation involved. Apex is placed at apexRetentionTime.
        /// </summary>
        private static ChromatographicPeak PeakOfWidth(double width, double apexRetentionTime, SpectraFileInfo file)
        {
            const int charge = 2;
            var identification = new Identification(
                file, "PEPTIDE", "PEPTIDE", 799.36, apexRetentionTime, charge, new List<ProteinGroup>());
            var peak = new ChromatographicPeak(identification, file);

            // Half maximum is crossed at +/- width/2, so the flanking points carry half the apex intensity.
            var points = new (double Rt, double Intensity)[]
            {
                (apexRetentionTime - width, 0.0),
                (apexRetentionTime - (width / 2), 50.0),
                (apexRetentionTime, 100.0),
                (apexRetentionTime + (width / 2), 50.0),
                (apexRetentionTime + width, 0.0)
            };

            peak.IsotopicEnvelopes = points
                .Select((p, i) => new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(799.36, p.Intensity, i, p.Rt), charge, p.Intensity * charge, 1.0))
                .ToList();
            peak.CalculateIntensityForThisFeature(false);

            return peak;
        }

        /// <summary>A peak that is still above half maximum where it ends, so it has no measurable width.</summary>
        private static ChromatographicPeak CensoredPeak(double apexRetentionTime, SpectraFileInfo file)
        {
            const int charge = 2;
            var identification = new Identification(
                file, "PEPTIDE", "PEPTIDE", 799.36, apexRetentionTime, charge, new List<ProteinGroup>());
            var peak = new ChromatographicPeak(identification, file);

            peak.IsotopicEnvelopes = new[] { 60.0, 80.0, 100.0, 90.0, 70.0 }
                .Select((intensity, i) => new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(799.36, intensity, i, apexRetentionTime + (i * 0.1)),
                    charge, intensity * charge, 1.0))
                .ToList();
            peak.CalculateIntensityForThisFeature(false);

            return peak;
        }

        private static SpectraFileInfo File(string name, int order) => new SpectraFileInfo(name + ".mzML", "c", order, 0, 0);

        /// <summary>A peak with no isotopic envelopes, so no apex and ApexRetentionTime of -1.</summary>
        private static ChromatographicPeak ApexlessPeak(SpectraFileInfo file)
        {
            var identification = new Identification(
                file, "PEPTIDE", "PEPTIDE", 799.36, 15.0, 2, new List<ProteinGroup>());
            var peak = new ChromatographicPeak(identification, file);
            peak.IsotopicEnvelopes = new List<IsotopicEnvelope>();
            return peak;
        }

        [Test]
        public void FilesSharingABaseNameAreSummarisedSeparately()
        {
            // Same file name, different directories -- two runs, not one.
            var runA = new SpectraFileInfo(@"C:\Data\RunA\sample.mzML", "c", 0, 0, 0);
            var runB = new SpectraFileInfo(@"C:\Data\RunB\sample.mzML", "c", 1, 0, 0);

            var peaks = new List<ChromatographicPeak>
            {
                PeakOfWidth(0.2, 10, runA), PeakOfWidth(0.2, 20, runA),
                PeakOfWidth(0.8, 10, runB), PeakOfWidth(0.8, 20, runB)
            };

            List<PeakWidthSummary> summaries = PeakWidthSummary.Summarize(peaks);

            Assert.Multiple(() =>
            {
                // Grouped by base name, these would be one bucket of four peaks with a median between
                // the two -- a number describing neither run.
                Assert.That(summaries, Has.Count.EqualTo(2));
                Assert.That(summaries.Select(x => x.TotalPeakCount), Is.EqualTo(new[] { 2, 2 }));
                Assert.That(summaries[0].MedianFullWidthAtHalfMaximum,
                    Is.EqualTo(0.2).Within(FloatRetentionTimeTolerance));
                Assert.That(summaries[1].MedianFullWidthAtHalfMaximum,
                    Is.EqualTo(0.8).Within(FloatRetentionTimeTolerance));
            });
        }

        [Test]
        public void PeaksWithNoApexAreLeftOutEntirely()
        {
            var file = File("a", 0);

            var peaks = new List<ChromatographicPeak>
            {
                ApexlessPeak(file),
                PeakOfWidth(0.4, 10, file),
                PeakOfWidth(0.4, 20, file)
            };

            List<PeakWidthSummary> summaries = PeakWidthSummary.Summarize(peaks);

            Assert.Multiple(() =>
            {
                // ApexRetentionTime is -1 without an apex, so the peak would sort first, drag
                // RetentionTimeStart to -1 and count toward the bin it cannot contribute a width to.
                Assert.That(summaries, Has.Count.EqualTo(1));
                Assert.That(summaries[0].TotalPeakCount, Is.EqualTo(2));
                Assert.That(summaries[0].RetentionTimeStart, Is.EqualTo(10.0).Within(FloatRetentionTimeTolerance));
            });
        }

        [Test]
        public void OneSummaryPerFileWithTheMedianOfThatFile()
        {
            var a = File("a", 0);
            var b = File("b", 1);

            var peaks = new List<ChromatographicPeak>
            {
                // median of 0.2, 0.4, 0.6 is 0.4
                PeakOfWidth(0.2, 10, a), PeakOfWidth(0.6, 20, a), PeakOfWidth(0.4, 30, a),
                // median of 1.0, 2.0 is 1.5
                PeakOfWidth(1.0, 10, b), PeakOfWidth(2.0, 20, b)
            };

            List<PeakWidthSummary> summaries = PeakWidthSummary.Summarize(peaks);

            Assert.That(summaries.Count, Is.EqualTo(2));
            Assert.That(summaries[0].FileName, Is.EqualTo("a"));
            Assert.That(summaries[0].MedianFullWidthAtHalfMaximum, Is.EqualTo(0.4).Within(FloatRetentionTimeTolerance));
            Assert.That(summaries[0].MeasuredPeakCount, Is.EqualTo(3));
            Assert.That(summaries[1].FileName, Is.EqualTo("b"));
            Assert.That(summaries[1].MedianFullWidthAtHalfMaximum, Is.EqualTo(1.5).Within(FloatRetentionTimeTolerance));
        }

        /// <summary>
        /// The point of binning: a file whose width drifts must not report one flat median. Ten narrow
        /// peaks early and ten wide peaks late average to something that describes neither.
        /// </summary>
        [Test]
        public void BinningSeparatesADriftingFile()
        {
            var file = File("drifting", 0);
            var peaks = new List<ChromatographicPeak>();

            for (int i = 0; i < 10; i++) peaks.Add(PeakOfWidth(0.2, 10 + i, file));
            for (int i = 0; i < 10; i++) peaks.Add(PeakOfWidth(0.8, 100 + i, file));

            List<PeakWidthSummary> unbinned = PeakWidthSummary.Summarize(peaks);
            Assert.That(unbinned.Count, Is.EqualTo(1));
            Assert.That(unbinned[0].MedianFullWidthAtHalfMaximum, Is.EqualTo(0.5).Within(FloatRetentionTimeTolerance),
                "the flat median sits between the two populations and describes neither");

            List<PeakWidthSummary> binned = PeakWidthSummary.Summarize(peaks, retentionTimeBins: 2);

            Assert.That(binned.Count, Is.EqualTo(2));
            Assert.That(binned[0].MedianFullWidthAtHalfMaximum, Is.EqualTo(0.2).Within(FloatRetentionTimeTolerance));
            Assert.That(binned[1].MedianFullWidthAtHalfMaximum, Is.EqualTo(0.8).Within(FloatRetentionTimeTolerance));
            Assert.That(binned[0].RetentionTimeEnd, Is.LessThan(binned[1].RetentionTimeStart));
        }

        /// <summary>
        /// Bins hold equal numbers of peaks, and a remainder joins the last bin rather than forming a
        /// short one whose median rests on a couple of peaks.
        /// </summary>
        [Test]
        public void RemainderGoesToTheLastBinRatherThanFormingAShortOne()
        {
            var file = File("seven", 0);
            var peaks = Enumerable.Range(0, 7).Select(i => PeakOfWidth(0.2, 10 + i, file)).ToList();

            List<PeakWidthSummary> binned = PeakWidthSummary.Summarize(peaks, retentionTimeBins: 3);

            Assert.That(binned.Count, Is.EqualTo(3));
            Assert.That(binned.Select(s => s.TotalPeakCount), Is.EqualTo(new[] { 2, 2, 3 }));
            Assert.That(binned.Sum(s => s.TotalPeakCount), Is.EqualTo(7));
        }

        /// <summary>
        /// Unmeasurable peaks are excluded from the median but still counted, so a caller can tell a
        /// median resting on most of the file from one resting on a fraction of it.
        /// </summary>
        [Test]
        public void CensoredPeaksAreCountedButDoNotMoveTheMedian()
        {
            var file = File("mixed", 0);
            var peaks = new List<ChromatographicPeak>
            {
                PeakOfWidth(0.2, 10, file), PeakOfWidth(0.4, 20, file), PeakOfWidth(0.6, 30, file),
                CensoredPeak(40, file), CensoredPeak(50, file)
            };

            PeakWidthSummary summary = PeakWidthSummary.Summarize(peaks).Single();

            Assert.That(summary.MedianFullWidthAtHalfMaximum, Is.EqualTo(0.4).Within(FloatRetentionTimeTolerance));
            Assert.That(summary.MeasuredPeakCount, Is.EqualTo(3));
            Assert.That(summary.TotalPeakCount, Is.EqualTo(5));
        }

        /// <summary>
        /// A bin in which nothing could be measured reports NaN, not zero. Zero would read as "these
        /// peaks are infinitely sharp" rather than "there is no measurement here".
        /// </summary>
        [Test]
        public void BinWithNothingMeasurableReportsNaNRatherThanZero()
        {
            var file = File("allCensored", 0);
            var peaks = new List<ChromatographicPeak> { CensoredPeak(10, file), CensoredPeak(20, file) };

            PeakWidthSummary summary = PeakWidthSummary.Summarize(peaks).Single();

            Assert.That(summary.MedianFullWidthAtHalfMaximum, Is.NaN);
            Assert.That(summary.MeasuredPeakCount, Is.EqualTo(0));
            Assert.That(summary.TotalPeakCount, Is.EqualTo(2));
        }

        [Test]
        public void QuartilesBracketTheMedian()
        {
            var file = File("quartiles", 0);
            var peaks = Enumerable.Range(1, 9).Select(i => PeakOfWidth(i * 0.1, 10 + i, file)).ToList();

            PeakWidthSummary summary = PeakWidthSummary.Summarize(peaks).Single();

            // widths 0.1 .. 0.9: median 0.5, quartiles 0.3 and 0.7
            Assert.That(summary.MedianFullWidthAtHalfMaximum, Is.EqualTo(0.5).Within(FloatRetentionTimeTolerance));
            Assert.That(summary.LowerQuartile, Is.EqualTo(0.3).Within(FloatRetentionTimeTolerance));
            Assert.That(summary.UpperQuartile, Is.EqualTo(0.7).Within(FloatRetentionTimeTolerance));
        }

        [Test]
        public void MoreBinsThanPeaksDoesNotProduceEmptyBins()
        {
            var file = File("two", 0);
            var peaks = new List<ChromatographicPeak> { PeakOfWidth(0.2, 10, file), PeakOfWidth(0.4, 20, file) };

            List<PeakWidthSummary> binned = PeakWidthSummary.Summarize(peaks, retentionTimeBins: 10);

            Assert.That(binned.Count, Is.EqualTo(2));
            Assert.That(binned.All(s => s.TotalPeakCount > 0));
        }

        [Test]
        public void EmptyInputGivesNoSummariesAndBadArgumentsThrow()
        {
            Assert.That(PeakWidthSummary.Summarize(new List<ChromatographicPeak>()), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => PeakWidthSummary.Summarize(null));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                PeakWidthSummary.Summarize(new List<ChromatographicPeak>(), retentionTimeBins: 0));
        }
    }
}
