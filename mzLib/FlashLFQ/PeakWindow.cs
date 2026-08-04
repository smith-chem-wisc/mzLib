using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;

namespace FlashLFQ
{
    /// <summary>
    /// A rectangular slice of the MS1 data surrounding one charge state of a single ChromatographicPeak. The slice
    /// spans the scans that charge state was observed in, and runs from mzExpansion below the lowest m/z observed
    /// for it to mzExpansion above the highest. Every MS1 peak in that region is retained, including peaks that
    /// belong to species other than the quantified precursor, which is the point of the type: it captures the
    /// co-eluting and interfering signal that quantification discards.
    ///
    /// A window covers one charge state and not the whole peak because the charge states of a peptide sit at
    /// unrelated m/z values. A window spanning all of them would either have to be a rectangle wide enough to
    /// enclose the gaps between them - hundreds of daltons of mostly irrelevant signal for a peak that was
    /// observed at 2+ and 3+ - or stop being a rectangle. The charge states are recombined by
    /// <see cref="QuantifiedMs1Feature"/>, which holds one window per charge state of a peak.
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

        /// <summary>
        /// The charge state every isotopic envelope in this window was traced at.
        /// </summary>
        public int ChargeState { get; }

        public float MinMz { get; }
        public float MaxMz { get; }

        /// <summary>
        /// The retention time of the most intense isotopic envelope observed at this charge state. This is the
        /// apex of the charge state, which is not the apex of the peak unless this is also its most intense
        /// charge state.
        /// </summary>
        public double ApexRetentionTime { get; }

        /// <summary>
        /// The intensity of the most intense isotopic envelope observed at this charge state.
        /// </summary>
        public double ApexIntensity { get; }

        /// <summary>
        /// The one based scan numbers of the raw file, one per scan the charge state was observed across,
        /// ascending. Not to be confused with the zero based indices carried by indexed peaks, which are numbered
        /// against an MS1-only array and address a different scan as soon as MS2 scans are interleaved.
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
        /// Indices into the matching array in Mz, marking the peaks FlashLFQ used to trace this charge state: the
        /// most abundant isotope peak of each isotopic envelope, at most one per scan. Stored as indices rather
        /// than a flag per peak because there is only a handful per window. The remaining isotope peaks of each
        /// envelope are not retained by ChromatographicPeak, so they appear in the window unmarked.
        /// </summary>
        public int[][] PeakfindingIndices { get; }

        public double MinRetentionTime => RetentionTimes[0];
        public double MaxRetentionTime => RetentionTimes[RetentionTimes.Length - 1];

        /// <summary>
        /// The total number of MS1 peaks inside the window, across all of its scans.
        /// </summary>
        public int PeakCount { get; }

        private PeakWindow(int chargeState, float minMz, float maxMz, double apexRetentionTime, double apexIntensity,
            int[] scanNumbers, double[] retentionTimes, float[][] mz, float[][] intensity, int[][] peakfindingIndices)
        {
            ChargeState = chargeState;
            MinMz = minMz;
            MaxMz = maxMz;
            ApexRetentionTime = apexRetentionTime;
            ApexIntensity = apexIntensity;
            ScanNumbers = scanNumbers;
            RetentionTimes = retentionTimes;
            Mz = mz;
            Intensity = intensity;
            PeakfindingIndices = peakfindingIndices;
            PeakCount = mz.Sum(a => a.Length);
        }

        /// <summary>
        /// Extracts the window of MS1 data surrounding one charge state of a chromatographic peak. Returns null if
        /// there are no envelopes to build a window around, or if one of them sits in a scan the map does not
        /// cover.
        /// </summary>
        /// <param name="envelopes"> the isotopic envelopes of a single charge state of one peak </param>
        /// <param name="dataFile"> the spectra file the peak was quantified from </param>
        /// <param name="ms1ScanNumberMap"> the map built by QuantifiedMs1Feature.BuildMs1ScanNumberMap for that same file </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z the window extends </param>
        /// <param name="useDynamicConnection"> read the scans through a dynamic connection the caller has already
        /// opened, rather than out of a statically loaded file </param>
        /// <exception cref="ArgumentException"> if the envelopes were not all traced at the same charge state </exception>
        public static PeakWindow Create(IReadOnlyList<IsotopicEnvelope> envelopes, MsDataFile dataFile,
            int[] ms1ScanNumberMap, double mzExpansion = DefaultMzExpansion, bool useDynamicConnection = false)
        {
            if (envelopes == null || envelopes.Count == 0 || dataFile == null || ms1ScanNumberMap == null)
                return null;

            int chargeState = envelopes[0].ChargeState;
            if (envelopes.Any(e => e.ChargeState != chargeState))
            {
                throw new ArgumentException("A peak window covers a single charge state, but envelopes from charge "
                    + "states " + string.Join(", ", envelopes.Select(e => e.ChargeState).Distinct().OrderBy(z => z))
                    + " were passed to it");
            }

            float minMz = (float)(envelopes.Min(e => e.IndexedPeak.M) - mzExpansion);
            float maxMz = (float)(envelopes.Max(e => e.IndexedPeak.M) + mzExpansion);

            int firstMs1Index = envelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = envelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            if (firstMs1Index < 0 || lastMs1Index >= ms1ScanNumberMap.Length)
                return null;

            // The m/z values FlashLFQ traced, grouped by the scan they were found in, so that the peaks it actually
            // used can be marked as each scan is sliced.
            var peakfindingMzByScanNumber = new Dictionary<int, HashSet<float>>();
            foreach (IsotopicEnvelope envelope in envelopes)
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

            double apexRetentionTime = double.NaN;
            int apexScanNumber = ms1ScanNumberMap[envelopes.MaxBy(e => e.Intensity).IndexedPeak.ZeroBasedScanIndex];

            for (int i = 0; i < scanCount; i++)
            {
                int oneBasedScanNumber = ms1ScanNumberMap[firstMs1Index + i];
                MsDataScan scan = GetScan(dataFile, oneBasedScanNumber, useDynamicConnection);
                peakfindingMzByScanNumber.TryGetValue(oneBasedScanNumber, out HashSet<float> peakfindingMz);

                scanNumbers[i] = oneBasedScanNumber;
                retentionTimes[i] = scan.RetentionTime;
                SliceScan(scan, minMz, maxMz, peakfindingMz, out mz[i], out intensity[i], out peakfindingIndices[i]);

                // The apex scan is one of the scans being read anyway, so its retention time is taken as it goes
                // past rather than by fetching the scan a second time.
                if (oneBasedScanNumber == apexScanNumber)
                {
                    apexRetentionTime = scan.RetentionTime;
                }
            }

            return new PeakWindow(chargeState, minMz, maxMz, apexRetentionTime,
                envelopes.Max(e => e.Intensity), scanNumbers, retentionTimes, mz, intensity, peakfindingIndices);
        }

        /// <summary>
        /// Reads one scan of a spectra file, through whichever kind of connection the caller opened. A dynamic
        /// connection seeks to the scan and reads it alone; the static path expects the whole file to have been
        /// loaded already, and quietly loads it if not.
        /// </summary>
        internal static MsDataScan GetScan(MsDataFile dataFile, int oneBasedScanNumber, bool useDynamicConnection)
        {
            return useDynamicConnection
                ? dataFile.GetOneBasedScanFromDynamicConnection(oneBasedScanNumber)
                : dataFile.GetOneBasedScan(oneBasedScanNumber);
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
        /// Writes the window as a single JSON object: the bounds of the charge state, followed by its scans as
        /// parallel arrays. The peak level metadata is not repeated here, as the window is written as part of the
        /// <see cref="QuantifiedMs1Feature"/> that carries it.
        /// </summary>
        public void WriteTo(Utf8JsonWriter writer)
        {
            writer.WriteStartObject();

            writer.WriteNumber("charge", ChargeState);
            writer.WriteNumber("minMz", MinMz);
            writer.WriteNumber("maxMz", MaxMz);
            writer.WriteNumber("rtStart", MinRetentionTime);
            writer.WriteNumber("rtApex", ApexRetentionTime);
            writer.WriteNumber("rtEnd", MaxRetentionTime);
            writer.WriteNumber("apexIntensity", ApexIntensity);
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

        public override string ToString()
        {
            return "+" + ChargeState + "; " + PeakCount + " peaks in " + ScanNumbers.Length + " scans; "
                + MinMz.ToString("F2") + "-" + MaxMz.ToString("F2") + " m/z; "
                + MinRetentionTime.ToString("F2") + "-" + MaxRetentionTime.ToString("F2") + " min";
        }
    }
}
