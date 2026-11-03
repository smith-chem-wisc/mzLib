using System;
using System.Collections.Concurrent;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// A parallel timsTOF load used to read every frame through one Bruker handle under one lock, so its workers queued
    /// behind each other. Each worker now reads through a handle of its own. These tests hold the new load to the old one
    /// - a verbatim copy of its parallel loops, which read through the shared handle - scan for scan and bit for bit, at
    /// one thread and at several.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    [ExcludeFromCodeCoverage]
    public sealed class TestTimsTofParallelHandles
    {
        private static readonly FilteringParams Filtering = new FilteringParams(numberOfPeaksToKeepPerWindow: 200, minimumAllowedIntensityRatioToBasePeak: 0.01);

        private static string DataFile(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", name);

        [TestCase("timsTOF_snippet.d", true)]
        [TestCase("timsTOF_snippet.d", false)]
        [TestCase("timsTOF_MRM.d", true)]
        public void LoadIsIdenticalToTheSharedHandleLoadAtOneAndSeveralThreads(string folder, bool filter)
        {
            string path = DataFile(folder);
            FilteringParams filteringParams = filter ? Filtering : null;
            string asBefore = Dump(LoadAsBefore(path, filteringParams, maxThreads: 10));
            Assert.That(asBefore.Length, Is.GreaterThan(1000), "nothing was loaded");

            foreach (int threads in new[] { 1, 4, 10 })
            {
                var reader = new TimsTofFileReader(path);
                reader.LoadAllStaticData(filteringParams, threads);
                Assert.That(Dump(reader), Is.EqualTo(asBefore), $"{folder} at {threads} threads");
            }
        }

        [Test]
        public void AFactoryWithItsOwnHandleSharesTheTablesAndReadsTheSameFrames()
        {
            var reader = new TimsTofFileReader(DataFile("timsTOF_snippet.d"));
            reader.InitiateDynamicConnection();
            ulong ownHandle = TimsTofFileReader.tims_open(BrukerFileReader.ConvertStringToUTF8ByteArray(reader.FilePath), 0);
            try
            {
                Assert.That(ownHandle, Is.Not.EqualTo(0UL));
                FrameProxyFactory shared = reader.FrameProxyFactory;
                FrameProxyFactory own = shared.WithOwnHandle(ownHandle);

                Assert.That(own.FileHandle, Is.EqualTo(ownHandle));
                Assert.That(own.FileLock, Is.Not.SameAs(shared.FileLock));
                Assert.That(own.FramesTable, Is.SameAs(shared.FramesTable));
                Assert.That(own.MzLookupArray, Is.SameAs(shared.MzLookupArray));
                Assert.That(own.OneOverK0LookupArray, Is.SameAs(shared.OneOverK0LookupArray));

                // every frame, read through both handles at once from several threads. Only the part of the buffer the
                // reader filled is compared - peak counts, then indices and intensities; GetScanRawData copies the whole
                // native buffer, whose tail past the data is whatever the allocation held.
                long[] frameIds = shared.FramesTable.OneBasedFrameIndex.ToArray();
                var differing = new ConcurrentBag<long>();
                Parallel.ForEach(frameIds, new ParallelOptions { MaxDegreeOfParallelism = 4 }, frameId =>
                {
                    if (!FilledPart(shared.GetFrameProxy(frameId)).SequenceEqual(FilledPart(own.GetFrameProxy(frameId))))
                    {
                        differing.Add(frameId);
                    }
                });
                Assert.That(differing, Is.Empty);
            }
            finally
            {
                TimsTofFileReader.tims_close(ownHandle);
                reader.CloseDynamicConnection();
            }
        }

        private static ReadOnlySpan<uint> FilledPart(FrameProxy frame)
        {
            long peaks = 0;
            for (int scan = 0; scan < frame.NumberOfScans; scan++)
            {
                peaks += frame._rawData[scan];
            }
            return frame._rawData.AsSpan(0, (int)(frame.NumberOfScans + 2 * peaks));
        }

        /// <summary>
        /// The parallel load as it was, verbatim: every worker builds its scans through the reader's shared handle, which
        /// is what BuildDDAScans and BuildMrmScan still do when given no frame reader of their own.
        /// </summary>
        private static TimsTofFileReader LoadAsBefore(string path, FilteringParams filteringParams, int maxThreads)
        {
            var reader = new TimsTofFileReader(path);
            reader.InitiateDynamicConnection();
            switch (reader.ScanMode)
            {
                case ScanMode.PASEF:
                    reader.Ms1ScansNoPrecursorsBag = new();
                    Parallel.ForEach(
                        Partitioner.Create(0, reader.Ms1FrameIds.Count),
                        new ParallelOptions() { MaxDegreeOfParallelism = maxThreads },
                        (range) =>
                        {
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                reader.BuildDDAScans(reader.Ms1FrameIds[i], filteringParams);
                            }
                        });
                    reader.AssignOneBasedPrecursorsToPasefScans();
                    break;
                case ScanMode.MRM:
                    // The implicit assumption here is that in MRM mode, no MS1 scans are collected
                    reader.MrmScanArray = new TimsDataScan[reader.NumberOfFrames];
                    Parallel.ForEach(
                        Partitioner.Create(0, reader.NumberOfFrames),
                        new ParallelOptions() { MaxDegreeOfParallelism = maxThreads },
                        (range) =>
                        {
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                reader.BuildMrmScan(i + 1, filteringParams); // i is zero-based, frame ids are one-based
                            }
                        });
                    reader.AssignScanNumbersToMrmScans();
                    break;
                default:
                    throw new InvalidOperationException("unexpected scan mode " + reader.ScanMode);
            }
            reader.CloseDynamicConnection();
            return reader;
        }

        /// <summary>
        /// Every field of every scan, with each double as its bits, so equal dumps mean bit-identical scans.
        /// </summary>
        internal static string Dump(MsDataFile file)
        {
            static string D(double? x) => x.HasValue ? BitConverter.DoubleToInt64Bits(x.Value).ToString("X16") : "null";
            var text = new StringBuilder();
            foreach (MsDataScan scan in file.Scans)
            {
                var tims = (TimsDataScan)scan;
                text.Append(scan.OneBasedScanNumber).Append('|').Append(scan.MsnOrder).Append('|').Append(scan.NativeId)
                    .Append('|').Append(D(scan.RetentionTime)).Append('|').Append(scan.Polarity).Append('|').Append(scan.MzAnalyzer)
                    .Append('|').Append(D(scan.ScanWindowRange.Minimum)).Append('|').Append(D(scan.ScanWindowRange.Maximum))
                    .Append('|').Append(scan.ScanFilter).Append('|').Append(scan.IsCentroid).Append('|').Append(D(scan.TotalIonCurrent))
                    .Append('|').Append(D(scan.InjectionTime)).Append('|').Append(D(scan.IsolationMz)).Append('|').Append(D(scan.IsolationWidth))
                    .Append('|').Append(scan.SelectedIonChargeStateGuess).Append('|').Append(D(scan.SelectedIonIntensity))
                    .Append('|').Append(D(scan.SelectedIonMZ)).Append('|').Append(scan.DissociationType)
                    .Append('|').Append(scan.OneBasedPrecursorScanNumber).Append('|').Append(D(scan.SelectedIonMonoisotopicGuessMz))
                    .Append('|').Append(D(scan.SelectedIonMonoisotopicGuessIntensity)).Append('|').Append(scan.HcdEnergy)
                    .Append('|').Append(tims.FrameId).Append('|').Append(tims.ScanNumberStart).Append('|').Append(tims.ScanNumberEnd)
                    .Append('|').Append(D(tims.OneOverK0)).Append('|').Append(tims.PrecursorId)
                    .Append('|').Append(tims.FrameIds == null ? "" : string.Join(",", tims.FrameIds)).Append('\n');
                foreach (double mz in scan.MassSpectrum.XArray)
                {
                    text.Append(D(mz)).Append(' ');
                }
                text.Append('\n');
                foreach (double intensity in scan.MassSpectrum.YArray)
                {
                    text.Append(D(intensity)).Append(' ');
                }
                text.Append('\n');
            }
            return text.ToString();
        }
    }
}
