using MassSpectrometry;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;

namespace FlashLFQ
{
    /// <summary>
    /// Everything the <see cref="Ms1FeatureData"/> of one spectra file have in common, written once ahead of them:
    /// what the file is, and one entry per MS1 scan those features were extracted from, giving the scan's one based
    /// number, its retention time and its total ion current.
    ///
    /// This is what lets a window hold nothing but m/z, intensity and a scan number. A retention time repeated
    /// beside every peak of every window is the same handful of values written thousands of times over - a scan
    /// covered by fifty features carries its retention time fifty times, once per peak in each - where a scan table
    /// per file writes it once and lets a scan number stand for it everywhere else.
    ///
    /// Only the scans the file's features actually cover are listed. The point is to resolve the scan numbers the
    /// windows refer to, and a file's features usually cover a small part of a run, so listing every MS1 scan of an
    /// hour long gradient would be mostly rows nothing points at.
    /// </summary>
    public class SpectraFileHeaderData
    {
        public SpectraFileInfo SpectraFileInfo { get; }

        /// <summary>
        /// The mass analyzer the covered scans were acquired on, or null if the reader did not say. Taken from a
        /// scan rather than from the source file, which describes the format the data was written in rather than
        /// the instrument it came off. Not every reader fills it in when a scan is read on its own: mzML holds the
        /// instrument configuration in its head, which a reader seeking to one scan of it never passes.
        /// </summary>
        public MZAnalyzerType? MzAnalyzer { get; }

        /// <summary>
        /// Every scan of the file, of every MS order, or null when the file was read scan by scan and so was never
        /// asked how many scans it holds. Readers that stream a file do not all know its length without walking it,
        /// and walking it to fill in a metadata field would undo the reason for streaming.
        /// </summary>
        public int? TotalScanCount { get; }

        /// <summary>
        /// Every MS1 scan of the file, not just the ones listed in Scans.
        /// </summary>
        public int Ms1ScanCount { get; }

        /// <summary>
        /// The retention time of the file's last MS1 scan, which is the length of the run as this output sees it.
        /// </summary>
        public double LastMs1RetentionTime { get; }

        /// <summary>
        /// How many features follow this header.
        /// </summary>
        public int FeatureCount { get; }

        /// <summary>
        /// The MS1 scans the features cover, ascending by scan number. Every scan number appearing in a window of
        /// a following feature appears here exactly once.
        /// </summary>
        public ScanInfo[] Scans { get; }

        private SpectraFileHeaderData(SpectraFileInfo spectraFileInfo, MZAnalyzerType? mzAnalyzer,
            int? totalScanCount, int ms1ScanCount, double lastMs1RetentionTime, int featureCount, ScanInfo[] scans)
        {
            SpectraFileInfo = spectraFileInfo;
            MzAnalyzer = mzAnalyzer;
            TotalScanCount = totalScanCount;
            Ms1ScanCount = ms1ScanCount;
            LastMs1RetentionTime = lastMs1RetentionTime;
            FeatureCount = featureCount;
            Scans = scans;
        }

        /// <summary>
        /// Summarizes the file and the scans its features will cover. Returns null when none of the peaks can be
        /// placed in the file at all, in which case no feature will be written for it either and a header would
        /// describe an empty section.
        ///
        /// Which scans the features cover is worked out from the peaks rather than from the features themselves,
        /// so the header can be written before the first feature is extracted, and the features can then be
        /// streamed out one at a time instead of being held until the header is known.
        /// </summary>
        /// <param name="spectraFileInfo"> the file the peaks were quantified from </param>
        /// <param name="dataFile"> that same file, opened </param>
        /// <param name="ms1ScanInfo"> the MS1 scans of the file, as Ms1FeatureData numbers them </param>
        /// <param name="peaks"> the peaks whose features are about to be written </param>
        /// <param name="useDynamicConnection"> read scans through a dynamic connection the caller has already
        /// opened, rather than out of a statically loaded file </param>
        public static SpectraFileHeaderData Create(SpectraFileInfo spectraFileInfo, MsDataFile dataFile,
            ScanInfo[] ms1ScanInfo, IEnumerable<ChromatographicPeak> peaks, bool useDynamicConnection = false)
        {
            if (spectraFileInfo == null || dataFile == null || ms1ScanInfo == null || ms1ScanInfo.Length == 0
                || peaks == null)
            {
                return null;
            }

            var covered = new bool[ms1ScanInfo.Length];
            int featureCount = 0;

            foreach (ChromatographicPeak peak in peaks)
            {
                bool peakIsCovered = false;
                foreach ((int firstScanIndex, int lastScanIndex) in Ms1FeatureData.GetScanIndexRanges(peak))
                {
                    // The same bounds check the window extraction makes, so that the scans listed here are exactly
                    // the scans the features end up referring to.
                    if (firstScanIndex < 0 || lastScanIndex >= ms1ScanInfo.Length)
                        continue;

                    for (int i = firstScanIndex; i <= lastScanIndex; i++)
                    {
                        covered[i] = true;
                    }
                    peakIsCovered = true;
                }

                if (peakIsCovered)
                {
                    featureCount++;
                }
            }

            ScanInfo[] scans = ms1ScanInfo.Where((_, i) => covered[i]).ToArray();
            if (scans.Length == 0)
                return null;

            // One scan is read to name the mass analyzer, which no scan metadata carries. It is one of the scans
            // the features are about to read anyway.
            MZAnalyzerType? mzAnalyzer = null;
            MsDataScan firstScan = PeakWindowData.GetScan(dataFile, scans[0].OneBasedScanNumber, useDynamicConnection);
            if (firstScan != null && firstScan.MzAnalyzer != MZAnalyzerType.Unknown)
            {
                mzAnalyzer = firstScan.MzAnalyzer;
            }

            return new SpectraFileHeaderData(spectraFileInfo, mzAnalyzer,
                dataFile.NumSpectra > 0 ? dataFile.NumSpectra : (int?)null, ms1ScanInfo.Length,
                ms1ScanInfo[ms1ScanInfo.Length - 1].RetentionTime, featureCount, scans);
        }

        /// <summary>
        /// Writes the header as a single JSON object: the file's metadata, followed by its covered scans as three
        /// parallel arrays.
        /// </summary>
        public void WriteTo(Utf8JsonWriter writer)
        {
            writer.WriteStartObject();

            writer.WriteString("type", "file");
            writer.WriteString("fileName", SpectraFileInfo.FilenameWithoutExtension);
            writer.WriteString("filePath", SpectraFileInfo.FullFilePathWithExtension);

            if (MzAnalyzer == null)
            {
                writer.WriteNull("mzAnalyzer");
            }
            else
            {
                writer.WriteString("mzAnalyzer", MzAnalyzer.ToString());
            }

            if (TotalScanCount == null)
            {
                writer.WriteNull("totalScanCount");
            }
            else
            {
                writer.WriteNumber("totalScanCount", TotalScanCount.Value);
            }

            writer.WriteNumber("ms1ScanCount", Ms1ScanCount);
            Ms1FeatureData.WriteNumberOrNull(writer, "lastMs1RetentionTime", LastMs1RetentionTime);
            writer.WriteNumber("featureCount", FeatureCount);
            writer.WriteNumber("scanCount", Scans.Length);

            writer.WriteStartArray("scanNumbers");
            foreach (ScanInfo scan in Scans)
            {
                writer.WriteNumberValue(scan.OneBasedScanNumber);
            }
            writer.WriteEndArray();

            writer.WriteStartArray("retentionTimes");
            foreach (ScanInfo scan in Scans)
            {
                writer.WriteNumberValue(scan.RetentionTime);
            }
            writer.WriteEndArray();

            // A scan whose total ion current was never recorded is written as null rather than as NaN, which is
            // not a JSON number and which Utf8JsonWriter refuses outright.
            writer.WriteStartArray("totalIonCurrents");
            foreach (ScanInfo scan in Scans)
            {
                if (double.IsNaN(scan.TotalIonCurrent) || double.IsInfinity(scan.TotalIonCurrent))
                {
                    writer.WriteNullValue();
                }
                else
                {
                    writer.WriteNumberValue(scan.TotalIonCurrent);
                }
            }
            writer.WriteEndArray();

            writer.WriteEndObject();
        }

        public override string ToString()
        {
            return SpectraFileInfo.FilenameWithoutExtension + "; " + FeatureCount + " features over "
                + Scans.Length + " of " + Ms1ScanCount + " MS1 scans";
        }
    }
}
