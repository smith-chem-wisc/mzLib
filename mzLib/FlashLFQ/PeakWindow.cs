using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;

namespace FlashLFQ
{
    /// <summary>
    /// A rectangular slice of the MS1 data surrounding a single ChromatographicPeak. The slice spans the scans the
    /// peak was observed in, and runs from mzExpansion below the lowest m/z observed for the peak to mzExpansion
    /// above the highest. Every MS1 peak in that region is retained, including peaks that belong to species other
    /// than the quantified precursor, which is the point of the type: it captures the co-eluting and interfering
    /// signal that quantification discards.
    ///
    /// The peaks are held as one m/z array and one intensity array per scan, not as peak objects. A window can
    /// cover tens of thousands of peaks and exists to be handed to a consumer or written to disk, so there is
    /// nothing to gain from an object per peak, and m/z and intensity are kept at the float precision the raw data
    /// is stored at anyway.
    ///
    /// The peaks are read straight out of the spectra file rather than out of a PeakIndexingEngine. The index is
    /// organized m/z-major and scan-minor, so a contiguous m/z slice of a handful of consecutive scans costs one
    /// binary search per m/z bin - 100 bins per dalton - where reading the scans costs one binary search per scan.
    /// </summary>
    public class PeakWindow
    {
        public const double DefaultMzExpansion = 0.5;

        public ChromatographicPeak Peak { get; }
        public float MinMz { get; }
        public float MaxMz { get; }
        public double ApexRetentionTime { get; }

        /// <summary>
        /// The one based scan numbers of the raw file, one per scan the peak was observed in, ascending. Not to be
        /// confused with the zero based indices carried by indexed peaks, which are numbered against an MS1-only
        /// array and address a different scan as soon as MS2 scans are interleaved.
        /// </summary>
        public int[] ScanNumbers { get; }

        /// <summary>
        /// Parallel to ScanNumbers.
        /// </summary>
        public double[] RetentionTimes { get; }

        /// <summary>
        /// One m/z array per entry in ScanNumbers, ascending within each scan. A scan with nothing inside the m/z
        /// window is kept, holding an empty array.
        /// </summary>
        public float[][] Mz { get; }

        /// <summary>
        /// One intensity array per entry in ScanNumbers, parallel to the matching array in Mz.
        /// </summary>
        public float[][] Intensity { get; }

        /// <summary>
        /// Indices into the matching array in Mz, marking the peaks FlashLFQ used to trace the feature: the most
        /// abundant isotope peak of each isotopic envelope, one per scan per charge state. Stored as indices rather
        /// than a flag per peak because there are only a handful per scan. The remaining isotope peaks of each
        /// envelope are not retained by ChromatographicPeak, so they appear in the window unmarked.
        /// </summary>
        public int[][] PeakfindingIndices { get; }

        public double MinRetentionTime => RetentionTimes[0];
        public double MaxRetentionTime => RetentionTimes[RetentionTimes.Length - 1];

        /// <summary>
        /// The total number of MS1 peaks inside the window, across all of its scans.
        /// </summary>
        public int PeakCount { get; }

        private PeakWindow(ChromatographicPeak peak, float minMz, float maxMz, double apexRetentionTime,
            int[] scanNumbers, double[] retentionTimes, float[][] mz, float[][] intensity, int[][] peakfindingIndices)
        {
            Peak = peak;
            MinMz = minMz;
            MaxMz = maxMz;
            ApexRetentionTime = apexRetentionTime;
            ScanNumbers = scanNumbers;
            RetentionTimes = retentionTimes;
            Mz = mz;
            Intensity = intensity;
            PeakfindingIndices = peakfindingIndices;
            PeakCount = mz.Sum(a => a.Length);
        }

        /// <summary>
        /// Maps the zero based scan indices carried by indexed peaks onto the one based scan numbers of the raw
        /// file. The MS1 scans are selected and ordered exactly as PeakIndexingEngine.InitializeIndexingEngine
        /// selects them, so that the indices recorded during quantification line up with this map.
        /// </summary>
        public static int[] BuildMs1ScanNumberMap(MsDataFile dataFile)
        {
            return dataFile.GetMS1Scans()
                .Where(s => s != null && s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .Select(s => s.OneBasedScanNumber)
                .ToArray();
        }

        /// <summary>
        /// Extracts the window of MS1 data surrounding a chromatographic peak. Returns null for peaks that have no
        /// isotopic envelopes, as those have no observed m/z or scan to build a window around.
        /// </summary>
        /// <param name="peak"> the quantified peak to build a window around </param>
        /// <param name="dataFile"> the loaded spectra file the peak was quantified from </param>
        /// <param name="ms1ScanNumberMap"> the map built by BuildMs1ScanNumberMap for that same file </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z the window extends </param>
        public static PeakWindow Create(ChromatographicPeak peak, MsDataFile dataFile, int[] ms1ScanNumberMap,
            double mzExpansion = DefaultMzExpansion)
        {
            if (peak == null || dataFile == null || ms1ScanNumberMap == null || !peak.IsotopicEnvelopes.Any())
                return null;

            float minMz = (float)(peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.M) - mzExpansion);
            float maxMz = (float)(peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.M) + mzExpansion);

            int firstMs1Index = peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            if (firstMs1Index < 0 || lastMs1Index >= ms1ScanNumberMap.Length)
                return null;

            // The m/z values FlashLFQ traced, grouped by the scan they were found in, so that the peaks it actually
            // used can be marked as each scan is sliced.
            var peakfindingMzByScanNumber = new Dictionary<int, HashSet<float>>();
            foreach (IsotopicEnvelope envelope in peak.IsotopicEnvelopes)
            {
                int scanNumber = ms1ScanNumberMap[envelope.IndexedPeak.ZeroBasedScanIndex];
                if (!peakfindingMzByScanNumber.TryGetValue(scanNumber, out HashSet<float> mzValues))
                {
                    mzValues = new HashSet<float>();
                    peakfindingMzByScanNumber[scanNumber] = mzValues;
                }
                mzValues.Add(envelope.IndexedPeak.M);
            }

            int scanCount = lastMs1Index - firstMs1Index + 1;
            var scanNumbers = new int[scanCount];
            var retentionTimes = new double[scanCount];
            var mz = new float[scanCount][];
            var intensity = new float[scanCount][];
            var peakfindingIndices = new int[scanCount][];

            for (int i = 0; i < scanCount; i++)
            {
                int oneBasedScanNumber = ms1ScanNumberMap[firstMs1Index + i];
                MsDataScan scan = dataFile.GetOneBasedScan(oneBasedScanNumber);
                peakfindingMzByScanNumber.TryGetValue(oneBasedScanNumber, out HashSet<float> peakfindingMz);

                scanNumbers[i] = oneBasedScanNumber;
                retentionTimes[i] = scan.RetentionTime;
                SliceScan(scan, minMz, maxMz, peakfindingMz, out mz[i], out intensity[i], out peakfindingIndices[i]);
            }

            double apexRetentionTime = double.NaN;
            if (peak.Apex != null)
            {
                int apexMs1Index = peak.Apex.IndexedPeak.ZeroBasedScanIndex;
                if (apexMs1Index >= 0 && apexMs1Index < ms1ScanNumberMap.Length)
                {
                    apexRetentionTime = dataFile.GetOneBasedScan(ms1ScanNumberMap[apexMs1Index]).RetentionTime;
                }
            }

            return new PeakWindow(peak, minMz, maxMz, apexRetentionTime, scanNumbers, retentionTimes, mz, intensity,
                peakfindingIndices);
        }

        /// <summary>
        /// Copies the peaks between minMz and maxMz, inclusive, out of a scan. The scan's m/z array is sorted
        /// ascending, so the slice is found with one binary search and read off contiguously.
        /// </summary>
        /// <param name="peakfindingMz"> m/z values FlashLFQ traced in this scan, or null if it traced none </param>
        private static void SliceScan(MsDataScan scan, float minMz, float maxMz, HashSet<float> peakfindingMz,
            out float[] mz, out float[] intensity, out int[] peakfindingIndices)
        {
            double[] mzArray = scan.MassSpectrum.XArray;
            double[] intensityArray = scan.MassSpectrum.YArray;

            int start = 0;
            int count = 0;
            if (mzArray.Length > 0)
            {
                // GetClosestIndex clamps to the ends of the array when the target falls outside it, so the returned
                // index is only a starting guess and has to be confirmed to be inside the window.
                start = mzArray.GetClosestIndex(minMz, ArraySearchOption.Next);
                while (start < mzArray.Length && mzArray[start] < minMz)
                {
                    start++;
                }

                int end = start;
                while (end < mzArray.Length && mzArray[end] <= maxMz)
                {
                    end++;
                }
                count = end - start;
            }

            if (count == 0)
            {
                mz = Array.Empty<float>();
                intensity = Array.Empty<float>();
                peakfindingIndices = Array.Empty<int>();
                return;
            }

            mz = new float[count];
            intensity = new float[count];
            List<int> foundPeakfindingIndices = null;

            for (int i = 0; i < count; i++)
            {
                mz[i] = (float)mzArray[start + i];
                intensity[i] = (float)intensityArray[start + i];

                // Indexed peaks store m/z narrowed from this same array, so comparing the narrowed value recovers
                // an exact match rather than needing a tolerance.
                if (peakfindingMz != null && peakfindingMz.Contains(mz[i]))
                {
                    foundPeakfindingIndices ??= new List<int>();
                    foundPeakfindingIndices.Add(i);
                }
            }

            peakfindingIndices = foundPeakfindingIndices?.ToArray() ?? Array.Empty<int>();
        }

        /// <summary>
        /// Writes the window as a single JSON object. The peak's metadata is written once, followed by the scans
        /// as parallel arrays, so the per-peak overhead does not scale with the number of MS1 peaks in the window.
        /// </summary>
        /// <param name="writer"> the writer to emit the object into </param>
        /// <param name="peakId"> an identifier distinguishing this window's peak from the file's other peaks </param>
        public void WriteTo(Utf8JsonWriter writer, int peakId)
        {
            writer.WriteStartObject();

            writer.WriteString("fileName", Peak.SpectraFileInfo.FilenameWithoutExtension);
            writer.WriteNumber("peakId", peakId);
            writer.WriteString("baseSequence", string.Join("|", Peak.Identifications.Select(p => p.BaseSequence).Distinct()));
            writer.WriteString("fullSequence", string.Join("|", Peak.Identifications.Select(p => p.ModifiedSequence).Distinct()));
            writer.WriteString("detectionType", Peak.DetectionType.ToString());
            writer.WriteNumber("charge", Peak.Apex?.ChargeState ?? 0);
            writer.WriteNumber("rtStart", MinRetentionTime);
            WriteNumberOrNull(writer, "rtApex", ApexRetentionTime);
            writer.WriteNumber("rtEnd", MaxRetentionTime);
            writer.WriteNumber("minMz", MinMz);
            writer.WriteNumber("maxMz", MaxMz);
            writer.WriteNumber("peakCount", PeakCount);

            writer.WriteStartArray("scanNumbers");
            foreach (int scanNumber in ScanNumbers)
            {
                writer.WriteNumberValue(scanNumber);
            }
            writer.WriteEndArray();

            writer.WriteStartArray("retentionTimes");
            foreach (double retentionTime in RetentionTimes)
            {
                writer.WriteNumberValue(retentionTime);
            }
            writer.WriteEndArray();

            WriteJaggedFloatArray(writer, "mz", Mz);
            WriteJaggedFloatArray(writer, "intensity", Intensity);

            writer.WriteStartArray("peakfindingIndices");
            foreach (int[] indices in PeakfindingIndices)
            {
                writer.WriteStartArray();
                foreach (int index in indices)
                {
                    writer.WriteNumberValue(index);
                }
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            writer.WriteEndObject();
        }

        private static void WriteJaggedFloatArray(Utf8JsonWriter writer, string propertyName, float[][] values)
        {
            writer.WriteStartArray(propertyName);
            foreach (float[] scanValues in values)
            {
                writer.WriteStartArray();
                foreach (float value in scanValues)
                {
                    writer.WriteNumberValue(value);
                }
                writer.WriteEndArray();
            }
            writer.WriteEndArray();
        }

        /// <summary>
        /// NaN is not representable in JSON, so a peak with no apex is written as null rather than as a value the
        /// reader would have to reject.
        /// </summary>
        private static void WriteNumberOrNull(Utf8JsonWriter writer, string propertyName, double value)
        {
            if (double.IsNaN(value) || double.IsInfinity(value))
            {
                writer.WriteNull(propertyName);
            }
            else
            {
                writer.WriteNumber(propertyName, value);
            }
        }

        public override string ToString()
        {
            return Peak.Identifications.First().ModifiedSequence + "; " + PeakCount + " peaks in " + ScanNumbers.Length
                + " scans; " + MinMz.ToString("F2") + "-" + MaxMz.ToString("F2") + " m/z; "
                + MinRetentionTime.ToString("F2") + "-" + MaxRetentionTime.ToString("F2") + " min";
        }
    }
}
