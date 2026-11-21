using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using SpectralAveraging;

namespace Test.AveragingTests
{
    /// <summary>
    /// The averaging code was rewritten for speed and must give bit-identical results. These tests run it side by side
    /// with verbatim copies of the code it replaced (<see cref="LegacySpectraAveraging"/>), over randomized spectra built to
    /// hit the awkward cases: NaN, +0 and -0, repeated intensities, bins holding one peak or several from one spectrum,
    /// every rejection, weighting and normalization type, and one thread against several.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestAveragingMatchesLegacy
    {
        private static readonly double[] AwkwardIntensities = { 0.0, -0.0, 1.0, 2.0, 2.0, 5.0, 1e6, 1e-300, double.NaN };

        private static double RandomIntensity(Random random)
        {
            // mostly ordinary values, often repeated, sometimes zero, -0 or NaN
            int pick = random.Next(20);
            if (pick < 4)
                return AwkwardIntensities[random.Next(AwkwardIntensities.Length - 1)];
            if (pick == 4)
                return random.Next(50) == 0 ? double.NaN : 3.0;
            return Math.Round(random.NextDouble() * 1000, random.Next(4));
        }

        private static double[] RandomValues(Random random, int count)
        {
            var values = new double[count];
            for (int i = 0; i < count; i++)
                values[i] = RandomIntensity(random);
            return values;
        }

        /// <summary>
        /// A group of spectra whose peaks fall on a shared set of m/z values with jitter, so that bins get one peak, a
        /// peak from every spectrum, or several from the same spectrum; occasionally a NaN m/z.
        /// </summary>
        private static (double[][] x, double[][] y) RandomSpectra(Random random, int spectra, double binSize)
        {
            double[] pool = Enumerable.Range(0, random.Next(1, 40)).Select(_ => 100 + random.NextDouble() * 1900).ToArray();
            var x = new double[spectra][];
            var y = new double[spectra][];
            for (int s = 0; s < spectra; s++)
            {
                var mzs = new List<double>();
                foreach (double mz in pool)
                {
                    if (random.NextDouble() < 0.7)
                        mzs.Add(mz + (random.NextDouble() - 0.5) * binSize * 0.6);
                    if (random.NextDouble() < 0.1)
                        mzs.Add(mz + (random.NextDouble() - 0.5) * binSize * 0.6);   // a second peak in the same bin
                }

                if (mzs.Count == 0 || random.NextDouble() < 0.2)
                    mzs.Add(100 + random.NextDouble() * 1900);                         // a peak no other spectrum has
                mzs.Sort();
                if (random.Next(100) == 0)
                    mzs[random.Next(mzs.Count)] = double.NaN;

                x[s] = mzs.ToArray();
                y[s] = mzs.Select(_ => RandomIntensity(random)).ToArray();
            }

            // make sure the awkward values turn up in a good share of groups, not once in a few thousand peaks
            double[] someSpectrum() => y[random.Next(spectra)];
            if (random.NextDouble() < 0.25) { var a = someSpectrum(); a[random.Next(a.Length)] = double.NaN; }
            if (random.NextDouble() < 0.4) { var a = someSpectrum(); a[random.Next(a.Length)] = -0.0; }
            if (random.NextDouble() < 0.5) { var a = someSpectrum(); var b = someSpectrum(); b[random.Next(b.Length)] = a[random.Next(a.Length)]; }

            return (x, y);
        }

        private static double[][] Copy(double[][] arrays) => arrays.Select(a => (double[])a.Clone()).ToArray();

        /// <summary>
        /// The same double, bit for bit - so 0 and -0 differ - except that any NaN equals any NaN. A NaN's sign and payload
        /// are not something either version controls: when two NaNs meet in one operation (a NaN intensity times a NaN
        /// weight, say), x64 passes on whichever is the first operand, and the JIT is free to swap the operands of a
        /// commutative operation, differently in different methods and tiers. The comparison caught exactly that, once in
        /// four runs: the same NaN with its sign bit set in one version and not in the other.
        /// </summary>
        private static bool SameDouble(double a, double b) =>
            BitConverter.DoubleToInt64Bits(a) == BitConverter.DoubleToInt64Bits(b) || (double.IsNaN(a) && double.IsNaN(b));

        private static bool BitIdentical(double[] a, double[] b) => FirstDifference(a, b) is null;

        /// <summary>Where two arrays first differ (see <see cref="SameDouble"/>), or null when they do not.</summary>
        private static string FirstDifference(double[] a, double[] b)
        {
            if (a.Length != b.Length)
                return $"lengths {a.Length} and {b.Length}";
            for (int i = 0; i < a.Length; i++)
                if (!SameDouble(a[i], b[i]))
                    return $"[{i}]: {a[i]} ({BitConverter.DoubleToInt64Bits(a[i]):X16}) and {b[i]} ({BitConverter.DoubleToInt64Bits(b[i]):X16})";
            return null;
        }

        /// <summary>Runs both, and requires the same values bit for bit, or the same exception type.</summary>
        private static void AssertSame<T>(Func<T> legacy, Func<T> current, Func<T, T, bool> equal, string what,
            Func<T, T, string> describe = null)
        {
            T expected = default, actual = default;
            Exception expectedException = null, actualException = null;
            try { expected = legacy(); } catch (Exception e) { expectedException = e; }
            try { actual = current(); } catch (Exception e) { actualException = e; }

            if (expectedException != null || actualException != null)
            {
                Assert.That(actualException?.GetType(), Is.EqualTo(expectedException?.GetType()), what);
                return;
            }

            if (!equal(expected, actual))
                Assert.Fail($"{what}: {describe?.Invoke(expected, actual)}");
        }

        [Test]
        public static void StreamingStandardDeviationMatchesMathNet()
        {
            var random = new Random(1);
            for (int trial = 0; trial < 20000; trial++)
            {
                double[] values = RandomValues(random, random.Next(0, 40));
                if (random.Next(10) == 0)
                    values = values.Select(v => v * 1e-3 + 0.1).ToArray();

                double expected = values.ToList().StandardDeviation();   // the list path the old code took
                double actual = OutlierRejection.StreamingStandardDeviation(values);

                Assert.That(SameDouble(actual, expected), Is.True,
                    $"trial {trial}: {actual} and {expected} for [{string.Join(", ", values)}]");
            }
        }

        [Test]
        public static void RejectOutliersOfAnArrayMatchesLegacy()
        {
            var random = new Random(2);
            foreach (OutlierRejectionType type in Enum.GetValues<OutlierRejectionType>())
            {
                var parameters = new SpectralAveragingParameters
                {
                    OutlierRejectionType = type, MinSigmaValue = 0.5 + random.NextDouble(), MaxSigmaValue = 1 + 2 * random.NextDouble(),
                    Percentile = 0.1 + 0.8 * random.NextDouble()
                };
                for (int trial = 0; trial < 2000; trial++)
                {
                    double[] values = RandomValues(random, random.Next(1, 30));
                    AssertSame(() => LegacyOutlierRejection.RejectOutliers((double[])values.Clone(), parameters),
                        () => OutlierRejection.RejectOutliers((double[])values.Clone(), parameters),
                        BitIdentical, $"{type}, trial {trial}: [{string.Join(", ", values)}]");
                }
            }
        }

        [Test]
        public static void RejectOutliersOfBinnedPeaksMatchesLegacy()
        {
            var random = new Random(3);
            foreach (OutlierRejectionType type in Enum.GetValues<OutlierRejectionType>())
            {
                var parameters = new SpectralAveragingParameters { OutlierRejectionType = type, MinSigmaValue = 0.5, MaxSigmaValue = 3 };
                for (int trial = 0; trial < 2000; trial++)
                {
                    var peaks = RandomValues(random, random.Next(1, 30))
                        .Select((intensity, i) => new BinnedPeak(7, 500 + i * 1e-3, intensity, random.Next(5))).ToList();

                    AssertSame(() => LegacyOutlierRejection.RejectOutliers(peaks.ToList(), parameters),
                        () => OutlierRejection.RejectOutliers(peaks.ToList(), parameters),
                        (a, b) => a.Count == b.Count && a.Zip(b).All(pair => SamePeak(pair.First, pair.Second)), $"{type}, trial {trial}");
                }
            }
        }

        private static bool SamePeak(BinnedPeak a, BinnedPeak b) =>
            a.Bin == b.Bin && a.SpectraId == b.SpectraId && SameDouble(a.Mz, b.Mz) && SameDouble(a.Intensity, b.Intensity);

        [Test]
        public static void GetBinsMatchesLegacy()
        {
            var random = new Random(5);
            for (int group = 0; group < 300; group++)
            {
                double binSize = new[] { 1e-5, 0.01, 0.1, 1.0 }[random.Next(4)];
                var (x, y) = RandomSpectra(random, random.Next(1, 8), binSize);

                var expected = LegacySpectraAveraging.GetBins(x, y, binSize).ToList();
                using var actual = SpectraAveraging.GetBins(x, y, binSize);

                // the bins in the order the dictionary enumerated them, each with the same peaks in the same order
                Assert.That(actual.Count, Is.EqualTo(expected.Count), $"group {group}");
                for (int b = 0; b < expected.Count; b++)
                {
                    var peaks = actual[b].ToArray();
                    Assert.That(peaks[0].Bin, Is.EqualTo(expected[b].Key), $"group {group}, bin {b}");
                    Assert.That(peaks.Length, Is.EqualTo(expected[b].Value.Count), $"group {group}, bin {b}");
                    for (int p = 0; p < peaks.Length; p++)
                        Assert.That(SamePeak(peaks[p], expected[b].Value[p]), Is.True, $"group {group}, bin {b}, peak {p}");
                }
            }
        }



        /// <summary>
        /// A DDA-like run: MS1 scans, some at another compensation voltage, each followed by a few MS2 scans; sometimes a
        /// scan number used twice, which is the one case where the order the unaveraged scans are added in matters.
        /// </summary>
        private static List<MsDataScan> RandomDdaRun(Random random)
        {
            var scans = new List<MsDataScan>();
            int scanNumber = 1;
            int ms1Count = random.Next(1, 14);
            for (int m = 0; m < ms1Count; m++)
            {
                var (x, y) = RandomSpectra(random, 1, 0.01);
                double? cv = random.Next(4) == 0 ? -45.0 : null;
                int ms1Number = scanNumber++;
                double[] mzs = x[0].Where(v => !double.IsNaN(v)).DefaultIfEmpty(500).ToArray();
                var spectrum = new MzSpectrum(mzs, mzs.Select(_ => RandomIntensity(random)).ToArray(), false);
                scans.Add(new MsDataScan(spectrum, ms1Number, 1, true, Polarity.Positive, m, spectrum.Range, null,
                    MZAnalyzerType.Orbitrap, 10, 1, null, $"controllerType=0 controllerNumber=1 scan={ms1Number}",
                    compensationVoltage: cv));

                for (int k = random.Next(0, 4); k > 0; k--)
                {
                    int number = random.Next(15) == 0 ? ms1Number : scanNumber++;   // now and then a repeated number
                    var ms2 = new MzSpectrum(new[] { 200.0, 300.0 + k }, new[] { 1.0, 2.0 }, false);
                    scans.Add(new MsDataScan(ms2, number, 2, true, Polarity.Positive, m + 0.1 * k, ms2.Range, null,
                        MZAnalyzerType.Orbitrap, 3, 1, null, $"controllerType=0 controllerNumber=1 scan={number}",
                        oneBasedPrecursorScanNumber: ms1Number, compensationVoltage: cv));
                }
            }

            return scans;
        }

        [Test]
        public static void AverageDdaScansMatchesLegacy()
        {
            var random = new Random(6);
            for (int run = 0; run < 300; run++)
            {
                var scans = RandomDdaRun(random);
                var parameters = new SpectralAveragingParameters
                {
                    SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScans,
                    OutlierRejectionType = OutlierRejectionType.SigmaClipping, MinSigmaValue = 0.5, MaxSigmaValue = 3,
                    NumberOfScansToAverage = random.Next(1, 6), ScanOverlap = 0, BinSize = 0.01
                };

                MsDataScan[] expected = null, actual = null;
                Exception expectedException = null, actualException = null;
                try { expected = LegacySpectraFileAveraging.AverageDdaScans(scans, parameters); } catch (Exception e) { expectedException = e; }
                try { actual = SpectraFileAveraging.AverageSpectraFile(scans, parameters); } catch (Exception e) { actualException = e; }
                Assert.That(actualException?.GetType(), Is.EqualTo(expectedException?.GetType()), $"run {run}");
                if (expectedException != null)
                    continue;

                Assert.That(actual.Length, Is.EqualTo(expected.Length), $"run {run}");
                for (int i = 0; i < expected.Length; i++)
                {
                    // the scans that were not averaged must be the very same objects, in the same places
                    Assert.That(scans.Contains(actual[i]), Is.EqualTo(scans.Contains(expected[i])), $"run {run}, scan {i}");
                    if (scans.Contains(expected[i]))
                    {
                        Assert.That(actual[i], Is.SameAs(expected[i]), $"run {run}, scan {i}");
                        continue;
                    }

                    Assert.That(actual[i].OneBasedScanNumber, Is.EqualTo(expected[i].OneBasedScanNumber), $"run {run}, scan {i}");
                    Assert.That(actual[i].NativeId, Is.EqualTo(expected[i].NativeId), $"run {run}, scan {i}");
                    Assert.That(BitIdentical(actual[i].MassSpectrum.XArray, expected[i].MassSpectrum.XArray), $"run {run}, scan {i}");
                    Assert.That(BitIdentical(actual[i].MassSpectrum.YArray, expected[i].MassSpectrum.YArray), $"run {run}, scan {i}");
                }
            }
        }

        [Test]
        public static void InMzOrderMatchesOrderBy()
        {
            var random = new Random(7);
            for (int trial = 0; trial < 20000; trial++)
            {
                // runs of bins, as GetBins numbers them: usually each run already in m/z order, with ties between runs
                int runs = random.Next(10) == 0 ? random.Next(17, 25) : random.Next(1, 8);
                var firstBinOfSpectrum = new int[runs + 1];
                var peaks = new List<(double mz, double intensity)>();
                double[] shared = Enumerable.Range(0, 6).Select(_ => Math.Round(100 + random.NextDouble() * 10, 1)).ToArray();
                for (int run = 0; run < runs; run++)
                {
                    firstBinOfSpectrum[run] = peaks.Count;
                    var mzs = Enumerable.Range(0, random.Next(0, 8))
                        .Select(_ => random.Next(3) == 0 ? shared[random.Next(shared.Length)] : Math.Round(100 + random.NextDouble() * 10, 2))
                        .ToList();
                    if (random.Next(50) == 0) mzs.Add(double.NaN);
                    if (random.Next(30) == 0) mzs.Add(-0.0);
                    if (random.Next(4) != 0) mzs.Sort();            // a run out of order now and then
                    foreach (double mz in mzs)
                    {
                        double intensity = random.Next(6) == 0 ? (random.Next(2) == 0 ? 0.0 : -0.0)
                            : random.Next(40) == 0 ? double.NaN : random.NextDouble();
                        peaks.Add((mz, intensity));
                    }
                }

                firstBinOfSpectrum[runs] = peaks.Count;
                var array = peaks.ToArray();

                var expected = array.Where(p => p.intensity != 0).OrderBy(p => p.mz).ToList();
                var actual = SpectraAveraging.InMzOrder(array, firstBinOfSpectrum);

                Assert.That(actual.Count, Is.EqualTo(expected.Count), $"trial {trial}");
                for (int i = 0; i < expected.Count; i++)
                {
                    Assert.That(BitConverter.DoubleToInt64Bits(actual[i].mz), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected[i].mz)), $"trial {trial}, peak {i}");
                    Assert.That(BitConverter.DoubleToInt64Bits(actual[i].intensity), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected[i].intensity)), $"trial {trial}, peak {i}");
                }
            }
        }

        [Test]
        public static void AverageSpectraMatchesLegacy()
        {
            var random = new Random(4);
            int compared = 0;
            for (int group = 0; group < 30; group++)
            {
                double binSize = new[] { 0.01, 0.1, 1.0 }[random.Next(3)];
                var (x, y) = RandomSpectra(random, random.Next(1, 8), binSize);

                foreach (OutlierRejectionType rejection in Enum.GetValues<OutlierRejectionType>())
                foreach (SpectraWeightingType weighting in Enum.GetValues<SpectraWeightingType>())
                foreach (NormalizationType normalization in Enum.GetValues<NormalizationType>())
                foreach (int threads in new[] { 1, 4 })
                {
                    var parameters = new SpectralAveragingParameters
                    {
                        OutlierRejectionType = rejection, SpectralWeightingType = weighting, NormalizationType = normalization,
                        SpectralAveragingType = SpectralAveragingType.MzBinning, BinSize = binSize,
                        MinSigmaValue = 0.5, MaxSigmaValue = 3, Percentile = 0.1, MaxThreadsToUsePerFile = threads
                    };

                    AssertSame(() => LegacySpectraAveraging.AverageSpectra(x, Copy(y), parameters),
                        () => SpectraAveraging.AverageSpectra(x, Copy(y), parameters),
                        (a, b) => BitIdentical(a[0], b[0]) && BitIdentical(a[1], b[1]),
                        $"group {group}, {rejection}/{weighting}/{normalization}, {threads} threads",
                        (a, b) => $"m/z {FirstDifference(a[0], b[0])}; intensity {FirstDifference(a[1], b[1])}");
                    compared++;
                }
            }

            Assert.That(compared, Is.EqualTo(30 * 7 * 3 * 4 * 2));
        }
    }
}
