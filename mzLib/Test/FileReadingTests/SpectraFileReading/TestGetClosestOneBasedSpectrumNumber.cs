using System;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// GetClosestOneBasedSpectrumNumber binary searches rather than walking every scan (#156). These pin
    /// the answers it must give, because a binary search is easy to write in a way that is right for
    /// interior values and wrong at the edges, on exact hits, or where scans share a retention time.
    ///
    /// Note the variables below are typed as MsDataFile, not FakeMsDataFile. FakeMsDataFile declares its
    /// own GetClosestOneBasedSpectrumNumber which HIDES rather than overrides the real one, and returns an
    /// insertion point instead of the closest scan. A test that held the derived type would quietly
    /// exercise the fake instead of the implementation.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestGetClosestOneBasedSpectrumNumber
    {
        private static MsDataFile FileWithRetentionTimes(params double[] retentionTimes)
        {
            var scans = new MsDataScan[retentionTimes.Length];
            for (int i = 0; i < retentionTimes.Length; i++)
            {
                var spectrum = new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false);
                scans[i] = new MsDataScan(spectrum, i + 1, 1, true, Polarity.Positive, retentionTimes[i],
                    new MzRange(50, 2000), "f", MZAnalyzerType.Orbitrap, 1.0, null, null, "scan=" + (i + 1));
            }
            return new FakeMsDataFile(scans);
        }

        [Test]
        public static void ExactRetentionTimeReturnsThatScan()
        {
            MsDataFile file = FileWithRetentionTimes(1.0, 2.0, 3.0, 4.0, 5.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(1.0), Is.EqualTo(1));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(3.0), Is.EqualTo(3));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(5.0), Is.EqualTo(5));
            });
        }

        /// <summary>
        /// Queries outside the acquired range clamp to the ends rather than running off them, which is
        /// where an off-by-one in a binary search usually shows up first.
        /// </summary>
        [Test]
        public static void RetentionTimesOutsideTheRangeClampToTheEnds()
        {
            MsDataFile file = FileWithRetentionTimes(10.0, 20.0, 30.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(-100), Is.EqualTo(1));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(9.99), Is.EqualTo(1));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(1000), Is.EqualTo(3));
            });
        }

        [Test]
        public static void NearerNeighbourWinsOnEitherSide()
        {
            MsDataFile file = FileWithRetentionTimes(10.0, 20.0, 30.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(12), Is.EqualTo(1), "12 is nearer 10 than 20");
                Assert.That(file.GetClosestOneBasedSpectrumNumber(18), Is.EqualTo(2), "18 is nearer 20 than 10");
                Assert.That(file.GetClosestOneBasedSpectrumNumber(29), Is.EqualTo(3));
            });
        }

        /// <summary>
        /// A retention time falling exactly between two scans resolves to the later one. This is not an
        /// arbitrary choice: it is what the linear implementation did, and callers such as
        /// GetMsScansInTimeRange were written against it.
        /// </summary>
        [Test]
        public static void ExactMidpointResolvesToTheLaterScan()
        {
            MsDataFile file = FileWithRetentionTimes(10.0, 20.0, 30.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(15.0), Is.EqualTo(2));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(25.0), Is.EqualTo(3));
            });
        }

        /// <summary>
        /// Scans can share a retention time, and the linear implementation kept walking while the distance
        /// stayed equal, so it came to rest on the last scan of the run. A plain binary search lands on the
        /// first instead. This is the case that caught the first version of the change: a differential test
        /// over random data disagreed on 19 of 33,056 queries, all of them runs of equal retention times.
        /// </summary>
        [Test]
        public static void RunsOfEqualRetentionTimesReturnTheLastScanInTheRun()
        {
            MsDataFile file = FileWithRetentionTimes(1.0, 5.0, 5.0, 5.0, 9.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(5.0), Is.EqualTo(4), "last scan of the run");
                Assert.That(file.GetClosestOneBasedSpectrumNumber(4.9), Is.EqualTo(4));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(5.1), Is.EqualTo(4));
            });
        }

        [Test]
        public static void SingleScanFileAlwaysReturnsThatScan()
        {
            MsDataFile file = FileWithRetentionTimes(7.5);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(0), Is.EqualTo(1));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(7.5), Is.EqualTo(1));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(100), Is.EqualTo(1));
            });
        }

        /// <summary>
        /// A file with no scans returns 0 rather than 1 or an exception, matching what the linear version
        /// returned when its loop never ran.
        /// </summary>
        [Test]
        public static void EmptyFileReturnsZero()
        {
            MsDataFile file = FileWithRetentionTimes();

            Assert.That(file.GetClosestOneBasedSpectrumNumber(5), Is.EqualTo(0));
        }
    }
}
