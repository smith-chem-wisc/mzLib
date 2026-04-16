using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MassSpectrometry
{
    public sealed class ImspExportService : IImspExportService
    {
        public const int DefaultBinsPerDalton = 100;
        public const double DefaultIntensityThreshold = 10000;

        public int WriteFile(IEnumerable<MsDataScan> ms1Scans, string outputPath, int binsPerDalton = DefaultBinsPerDalton,
            double intensityThreshold = DefaultIntensityThreshold)
        {
            if (string.IsNullOrWhiteSpace(outputPath))
            {
                throw new ArgumentException("Output path is required.", nameof(outputPath));
            }

            using (var fs = File.Create(outputPath))
            {
                return WriteStream(ms1Scans, fs, binsPerDalton, intensityThreshold);
            }
        }

        public byte[] WriteBytes(IEnumerable<MsDataScan> ms1Scans, int binsPerDalton = DefaultBinsPerDalton,
            double intensityThreshold = DefaultIntensityThreshold)
        {
            using (var stream = new MemoryStream())
            {
                WriteStream(ms1Scans, stream, binsPerDalton, intensityThreshold);
                return stream.ToArray();
            }
        }

        public int WriteStream(IEnumerable<MsDataScan> ms1Scans, Stream outputStream, int binsPerDalton = DefaultBinsPerDalton,
            double intensityThreshold = DefaultIntensityThreshold)
        {
            if (ms1Scans == null)
            {
                throw new ArgumentNullException(nameof(ms1Scans));
            }

            if (outputStream == null)
            {
                throw new ArgumentNullException(nameof(outputStream));
            }

            if (!outputStream.CanWrite)
            {
                throw new ArgumentException("Output stream must be writable.", nameof(outputStream));
            }

            if (binsPerDalton <= 0)
            {
                throw new ArgumentOutOfRangeException(nameof(binsPerDalton), "Bins per dalton must be positive.");
            }

            var scans = ms1Scans.ToArray();
            if (scans.Length == 0)
            {
                throw new ArgumentException("At least one MS1 scan is required.", nameof(ms1Scans));
            }

            double maxMz = scans
                .Where(s => s.MassSpectrum.LastX.HasValue)
                .Select(s => s.MassSpectrum.LastX.Value)
                .DefaultIfEmpty(double.NaN)
                .Max();

            if (double.IsNaN(maxMz))
            {
                throw new ArgumentException("At least one scan must contain m/z values.", nameof(ms1Scans));
            }

            int binArrayLength = (int)Math.Ceiling(maxMz * binsPerDalton) + 1;
            var bins = new List<ImspPeak>[binArrayLength];

            for (int scanIndex = 0; scanIndex < scans.Length; scanIndex++)
            {
                var scan = scans[scanIndex];
                for (int peakIndex = 0; peakIndex < scan.MassSpectrum.XArray.Length; peakIndex++)
                {
                    double intensity = scan.MassSpectrum.YArray[peakIndex];
                    if (intensity < intensityThreshold)
                    {
                        continue;
                    }

                    double mz = scan.MassSpectrum.XArray[peakIndex];
                    int binIndex = (int)Math.Round(mz * binsPerDalton, 0);
                    if (bins[binIndex] == null)
                    {
                        bins[binIndex] = new List<ImspPeak>();
                    }

                    bins[binIndex].Add(new ImspPeak(mz, intensity, scanIndex));
                }
            }

            var nonEmptyBins = Enumerable.Range(0, binArrayLength)
                .Where(i => bins[i] != null)
                .Select(i => new ImspBin(i, bins[i]))
                .ToArray();

            int totalPeakCount = nonEmptyBins.Sum(b => b.Peaks.Count);

            var ticPerScan = new double[scans.Length];
            foreach (var bin in nonEmptyBins)
            {
                foreach (var peak in bin.Peaks)
                {
                    ticPerScan[peak.ScanIndex] += peak.Intensity;
                }
            }

            using (var bw = new BinaryWriter(outputStream, System.Text.Encoding.UTF8, leaveOpen: true))
            {
                bw.Write(new[] { (byte)'I', (byte)'M', (byte)'S', (byte)'P' });
                bw.Write((uint)1);
                bw.Write((uint)binsPerDalton);
                bw.Write((uint)nonEmptyBins.Length);
                bw.Write((uint)totalPeakCount);
                bw.Write((uint)scans.Length);

                for (int scanIndex = 0; scanIndex < scans.Length; scanIndex++)
                {
                    bw.Write((uint)scans[scanIndex].OneBasedScanNumber);
                    bw.Write(scans[scanIndex].RetentionTime);
                    bw.Write((float)ticPerScan[scanIndex]);
                }

                uint peakOffset = 0;
                foreach (var bin in nonEmptyBins)
                {
                    bw.Write((uint)bin.BinIndex);
                    bw.Write(peakOffset);
                    bw.Write((uint)bin.Peaks.Count);
                    peakOffset += (uint)bin.Peaks.Count;
                }

                foreach (var bin in nonEmptyBins)
                {
                    foreach (var peak in bin.Peaks)
                    {
                        bw.Write((uint)Math.Round(peak.Mz * 10000));
                        bw.Write((float)peak.Intensity);
                        bw.Write((uint)peak.ScanIndex);
                    }
                }
            }

            return totalPeakCount;
        }

        private sealed class ImspBin
        {
            public ImspBin(int binIndex, List<ImspPeak> peaks)
            {
                BinIndex = binIndex;
                Peaks = peaks;
            }

            public int BinIndex { get; }
            public List<ImspPeak> Peaks { get; }
        }

        private readonly struct ImspPeak
        {
            public ImspPeak(double mz, double intensity, int scanIndex)
            {
                Mz = mz;
                Intensity = intensity;
                ScanIndex = scanIndex;
            }

            public double Mz { get; }
            public double Intensity { get; }
            public int ScanIndex { get; }
        }
    }
}
