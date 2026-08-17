using System.Diagnostics.CodeAnalysis;
using System.Linq;
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
    /// These deliberately do not use FakeMsDataFile. It declares its own GetClosestOneBasedSpectrumNumber
    /// which HIDES rather than overrides the real one and returns an insertion point instead of the closest
    /// scan, so a test holding that type would quietly exercise the double instead of the implementation.
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
            return new AlreadyLoadedFile(scans);
        }

        private static MsDataFile FileWithGappedScanNumbers(params (int scanNumber, double retentionTime)[] scans)
        {
            var dense = new MsDataScan[scans.Length];
            for (int i = 0; i < scans.Length; i++)
            {
                var spectrum = new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false);
                dense[i] = new MsDataScan(spectrum, scans[i].scanNumber, 1, true, Polarity.Positive,
                    scans[i].retentionTime, new MzRange(50, 2000), "f", MZAnalyzerType.Orbitrap, 1.0, null, null,
                    "scan=" + scans[i].scanNumber);
            }
            return new SparselyIndexedFile(dense);
        }

        /// <summary>
        /// A file whose scans are already in memory, so LoadAllStaticData is a no-op rather than a throw.
        /// FakeMsDataFile cannot be used here for two reasons: it hides GetClosestOneBasedSpectrumNumber
        /// with its own version, and its LoadAllStaticData throws. The second matters for the zero-scan
        /// case, because CheckIfScansLoaded is "Scans != null &amp;&amp; Scans.Length > 0" - a file with no
        /// scans always counts as not loaded, so the method tries to load before it can return anything.
        /// </summary>
        private class AlreadyLoadedFile : MsDataFile
        {
            public AlreadyLoadedFile(MsDataScan[] scans)
                : base(scans, new SourceFile(@"scan number only nativeID format", "mzML format", null,
                    "SHA-1", @"C:\fake.mzML", null))
            {
            }

            public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1) => this;
            public override SourceFile GetSourceFile() => SourceFile;
            public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber,
                IFilteringParams filterParams = null) => GetOneBasedScan(oneBasedScanNumber);
            public override void CloseDynamicConnection() { }
            public override void InitiateDynamicConnection() { }
        }

        /// <summary>
        /// Stores scans the way Mgf and MsAlign do: Scans is dense, and GetOneBasedScan reads a separate
        /// array sized to the largest scan number with nulls wherever a number is missing.
        /// </summary>
        private sealed class SparselyIndexedFile : AlreadyLoadedFile
        {
            private readonly MsDataScan[] _indexedScans;

            public SparselyIndexedFile(MsDataScan[] dense) : base(dense)
            {
                _indexedScans = new MsDataScan[dense[^1].OneBasedScanNumber];
                foreach (var scan in dense)
                    _indexedScans[scan.OneBasedScanNumber - 1] = scan;
            }

            public override MsDataScan GetOneBasedScan(int scanNumber) => _indexedScans[scanNumber - 1];
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
        /// A file that loads successfully but contains no scans returns 0 rather than 1 or an exception,
        /// matching what the linear version returned when its loop never ran. Reaching this requires a
        /// LoadAllStaticData that succeeds without producing scans, because CheckIfScansLoaded treats a
        /// zero-scan file as not loaded and the method therefore attempts a load first.
        /// </summary>
        [Test]
        public static void FileThatLoadsWithNoScansReturnsZero()
        {
            MsDataFile file = FileWithRetentionTimes();

            Assert.That(file.GetClosestOneBasedSpectrumNumber(5), Is.EqualTo(0));
        }

        /// <summary>
        /// A run of equal retention times that ends at the last scan. RunsOfEqualRetentionTimes... keeps its
        /// run strictly interior, so the end of the search range is never the answer there and the bound on
        /// the run search is never what stops it. Here it is: the run must not be walked off the end.
        /// </summary>
        [Test]
        public static void RunOfEqualRetentionTimesAtTheEndOfTheFileReturnsTheLastScan()
        {
            MsDataFile file = FileWithRetentionTimes(1.0, 5.0, 9.0, 9.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(9.0), Is.EqualTo(4), "last scan of a trailing run");
                Assert.That(file.GetClosestOneBasedSpectrumNumber(100), Is.EqualTo(4), "target past the end");
            });
        }

        /// <summary>
        /// A run of equal retention times that starts at the first scan, so the predecessor comparison is
        /// skipped entirely and the run search begins from the first index. The interior fixture never takes
        /// that path, because there is always an earlier scan to compare against.
        /// </summary>
        [Test]
        public static void RunOfEqualRetentionTimesAtTheStartOfTheFileReturnsTheLastScanInTheRun()
        {
            MsDataFile file = FileWithRetentionTimes(5.0, 5.0, 9.0);

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(0), Is.EqualTo(2), "target before a leading run");
                Assert.That(file.GetClosestOneBasedSpectrumNumber(5.0), Is.EqualTo(2), "exact hit on a leading run");
            });
        }

        /// <summary>
        /// GetMsScansInTimeRange seeds from GetClosestOneBasedSpectrumNumber, so the tie-break and the
        /// run answer above reach it directly. Pinning the scan numbers it yields is what makes "range
        /// queries are unaffected" a checked claim rather than an argument.
        ///
        /// Note the third case: the window covers scans 2, 3 and 4, but only 4 and 5 come back. The seed is
        /// the last scan of the equal-retention-time run and the walk only goes forward, so the earlier
        /// members are skipped. That is pre-existing, is not changed by the binary search, and is tracked
        /// as #1172 - it is asserted here so that fixing it has to be a deliberate act.
        /// </summary>
        /// <summary>
        /// Mgf and MsAlign keep a dense Scans array but override GetOneBasedScan to read a separate array
        /// sized to the largest scan number, null in the gaps (Mgf.cs:62-67). NumSpectra counts the dense
        /// one, so on a file whose scan numbers skip, the holes land inside the searched range [1,
        /// NumSpectra] and probing one throws. Both assertions here throw NullReferenceException if the
        /// search goes back through that accessor; the second also pins that the number handed back is
        /// the one the reader's own accessor wants, not a position in Scans.
        /// </summary>
        [Test]
        public static void GappedScanNumbersAreSearchableAndTheAnswerRoundTrips()
        {
            MsDataFile file = FileWithGappedScanNumbers((1, 1.0), (2, 2.0), (3, 3.0), (10, 10.0), (11, 11.0));

            Assert.Multiple(() =>
            {
                Assert.That(file.GetClosestOneBasedSpectrumNumber(2.4), Is.EqualTo(2));
                Assert.That(file.GetClosestOneBasedSpectrumNumber(10.6), Is.EqualTo(11), "across the gap");
                Assert.That(file.GetOneBasedScan(file.GetClosestOneBasedSpectrumNumber(10.6)).RetentionTime,
                    Is.EqualTo(11.0), "the answer round-trips through the reader's accessor");
            });
        }

        [Test]
        public static void RangeQueriesSeededByThisMethodYieldTheExpectedScans()
        {
            MsDataFile spread = FileWithRetentionTimes(1.0, 3.0, 5.0, 7.0, 9.0);
            MsDataFile withRun = FileWithRetentionTimes(1.0, 5.0, 5.0, 5.0, 9.0);

            Assert.Multiple(() =>
            {
                Assert.That(spread.GetMsScansInTimeRange(3.0, 7.0).Select(s => s.OneBasedScanNumber),
                    Is.EqualTo(new[] { 2, 3, 4 }));
                Assert.That(spread.GetMsScansInTimeRange(4.0, 6.0).Select(s => s.OneBasedScanNumber),
                    Is.EqualTo(new[] { 3 }), "midpoint seed resolves to the later scan");
                Assert.That(withRun.GetMsScansInTimeRange(5.0, 10.0).Select(s => s.OneBasedScanNumber),
                    Is.EqualTo(new[] { 4, 5 }), "#1172: scans 2 and 3 are in the window but are skipped");
            });
        }
    }
}
