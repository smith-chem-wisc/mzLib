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
    /// <see cref="Ms1FeatureData"/>, which holds one window per charge state of a peak.
    ///
    /// The peaks are held as three parallel flat arrays - m/z, intensity, and the one based scan number each peak
    /// was found in - rather than as peak objects, or as one array per scan. A window can cover tens of thousands
    /// of peaks and exists to be handed to a consumer or written to disk, so there is nothing to gain from an
    /// object per peak, and m/z and intensity are kept at the float precision the raw data is stored at anyway.
    /// Everything a scan number implies - its retention time, its total ion current - is written once per scan in
    /// the <see cref="SpectraFileHeaderData"/> heading the file's features, so the window carries none of it.
    ///
    /// The peaks are read straight out of the spectra file rather than out of a PeakIndexingEngine. The index is
    /// organized m/z-major and scan-minor, so a contiguous m/z slice of a handful of consecutive scans costs one
    /// binary search per m/z bin - 100 bins per dalton - where reading the scans costs one binary search per scan.
    /// </summary>
    public class PeakWindowData
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
        /// The one based scan numbers of the raw file bounding the window, inclusive. Not to be confused with the
        /// zero based indices carried by indexed peaks, which are numbered against an MS1-only array and address a
        /// different scan as soon as MS2 scans are interleaved. The window covers every MS1 scan between the two,
        /// including any that turned out to hold nothing inside its m/z range and so appear in none of the arrays.
        /// </summary>
        public int FirstScanNumber { get; }
        public int LastScanNumber { get; }

        public double MinRetentionTime { get; }
        public double MaxRetentionTime { get; }

        /// <summary>
        /// The m/z of every MS1 peak in the window, ascending within each scan and grouped by scan in ascending
        /// scan order.
        /// </summary>
        public float[] Mz { get; }

        /// <summary>
        /// Parallel to Mz.
        /// </summary>
        public float[] Intensity { get; }

        /// <summary>
        /// Parallel to Mz: the one based scan number each peak was found in. Scan numbers repeat, once per peak
        /// the scan contributed, which is what lets the three arrays be flat rather than one array per scan.
        /// </summary>
        public int[] ScanNumbers { get; }

        /// <summary>
        /// Indices into Mz, marking the peaks FlashLFQ used to trace this charge state: the most abundant isotope
        /// peak of each isotopic envelope, at most one per scan. Stored as indices rather than a flag per peak
        /// because there is only a handful per window. The remaining isotope peaks of each envelope are not
        /// retained by ChromatographicPeak, so they appear in the window unmarked.
        /// </summary>
        public int[] PeakfindingIndices { get; }

        /// <summary>
        /// The total number of MS1 peaks inside the window, across all of its scans.
        /// </summary>
        public int PeakCount => Mz.Length;

        /// <summary>
        /// The number of MS1 scans the window spans, counting scans that held nothing inside its m/z range.
        /// </summary>
        public int ScanCount { get; }

        private PeakWindowData(int chargeState, float minMz, float maxMz, double apexRetentionTime,
            double apexIntensity, int firstScanNumber, int lastScanNumber, double minRetentionTime,
            double maxRetentionTime, int scanCount, float[] mz, float[] intensity, int[] scanNumbers,
            int[] peakfindingIndices)
        {
            ChargeState = chargeState;
            MinMz = minMz;
            MaxMz = maxMz;
            ApexRetentionTime = apexRetentionTime;
            ApexIntensity = apexIntensity;
            FirstScanNumber = firstScanNumber;
            LastScanNumber = lastScanNumber;
            MinRetentionTime = minRetentionTime;
            MaxRetentionTime = maxRetentionTime;
            ScanCount = scanCount;
            Mz = mz;
            Intensity = intensity;
            ScanNumbers = scanNumbers;
            PeakfindingIndices = peakfindingIndices;
        }

        /// <summary>
        /// Extracts the window of MS1 data surrounding one charge state of a chromatographic peak. Returns null if
        /// there are no envelopes to build a window around, or if one of them sits in a scan the scan metadata does
        /// not cover.
        /// </summary>
        /// <param name="envelopes"> the isotopic envelopes of a single charge state of one peak </param>
        /// <param name="dataFile"> the spectra file the peak was quantified from </param>
        /// <param name="ms1ScanInfo"> the MS1 scans of that same file, as built by Ms1FeatureData.BuildMs1ScanInfo
        /// or left behind by the run that quantified the peak </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z the window extends </param>
        /// <param name="useDynamicConnection"> read the scans through a dynamic connection the caller has already
        /// opened, rather than out of a statically loaded file </param>
        /// <exception cref="ArgumentException"> if the envelopes were not all traced at the same charge state </exception>
        public static PeakWindowData Create(IReadOnlyList<IsotopicEnvelope> envelopes, MsDataFile dataFile,
            ScanInfo[] ms1ScanInfo, double mzExpansion = DefaultMzExpansion, bool useDynamicConnection = false)
        {
            if (envelopes == null || envelopes.Count == 0 || dataFile == null || ms1ScanInfo == null)
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
            if (firstMs1Index < 0 || lastMs1Index >= ms1ScanInfo.Length)
                return null;

            // The m/z values FlashLFQ traced, grouped by the scan they were found in, so that the peaks it actually
            // used can be marked as each scan is sliced.
            var peakfindingMzByScanNumber = new Dictionary<int, HashSet<float>>();
            foreach (IsotopicEnvelope envelope in envelopes)
            {
                int scanNumber = ms1ScanInfo[envelope.IndexedPeak.ZeroBasedScanIndex].OneBasedScanNumber;
                if (!peakfindingMzByScanNumber.TryGetValue(scanNumber, out HashSet<float> mzValues))
                {
                    mzValues = new HashSet<float>();
                    peakfindingMzByScanNumber[scanNumber] = mzValues;
                }
                mzValues.Add(envelope.IndexedPeak.M);
            }

            int scanCount = lastMs1Index - firstMs1Index + 1;
            var mz = new List<float>();
            var intensity = new List<float>();
            var scanNumbers = new List<int>();
            var peakfindingIndices = new List<int>();

            for (int i = 0; i < scanCount; i++)
            {
                int oneBasedScanNumber = ms1ScanInfo[firstMs1Index + i].OneBasedScanNumber;
                MsDataScan scan = GetScan(dataFile, oneBasedScanNumber, useDynamicConnection);
                peakfindingMzByScanNumber.TryGetValue(oneBasedScanNumber, out HashSet<float> peakfindingMz);

                AppendScan(scan, oneBasedScanNumber, minMz, maxMz, peakfindingMz, mz, intensity, scanNumbers,
                    peakfindingIndices);
            }

            // Every retention time the window reports is one the scan metadata already holds, so resolving them
            // costs no reads of its own - not even for the apex, which need not be a scan the window covers when
            // the charge state was traced across a gap.
            int apexMs1Index = envelopes.MaxBy(e => e.Intensity).IndexedPeak.ZeroBasedScanIndex;

            return new PeakWindowData(chargeState, minMz, maxMz, ms1ScanInfo[apexMs1Index].RetentionTime,
                envelopes.Max(e => e.Intensity), ms1ScanInfo[firstMs1Index].OneBasedScanNumber,
                ms1ScanInfo[lastMs1Index].OneBasedScanNumber, ms1ScanInfo[firstMs1Index].RetentionTime,
                ms1ScanInfo[lastMs1Index].RetentionTime, scanCount, mz.ToArray(), intensity.ToArray(),
                scanNumbers.ToArray(), peakfindingIndices.ToArray());
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
        /// Appends the peaks between minMz and maxMz, inclusive, of one scan to the window's arrays. The scan's m/z
        /// array is sorted ascending, so the slice is found with one binary search and read off contiguously.
        /// </summary>
        /// <param name="peakfindingMz"> m/z values FlashLFQ traced in this scan, or null if it traced none </param>
        private static void AppendScan(MsDataScan scan, int oneBasedScanNumber, float minMz, float maxMz,
            HashSet<float> peakfindingMz, List<float> mz, List<float> intensity, List<int> scanNumbers,
            List<int> peakfindingIndices)
        {
            double[] mzArray = scan.MassSpectrum.XArray;
            double[] intensityArray = scan.MassSpectrum.YArray;

            if (mzArray.Length == 0)
                return;

            // GetClosestIndex clamps to the ends of the array when the target falls outside it, so the returned
            // index is only a starting guess and has to be confirmed to be inside the window.
            int start = mzArray.GetClosestIndex(minMz, ArraySearchOption.Next);
            while (start < mzArray.Length && mzArray[start] < minMz)
            {
                start++;
            }

            for (int i = start; i < mzArray.Length && mzArray[i] <= maxMz; i++)
            {
                float peakMz = (float)mzArray[i];

                // Indexed peaks store m/z narrowed from this same array, so comparing the narrowed value recovers
                // an exact match rather than needing a tolerance.
                if (peakfindingMz != null && peakfindingMz.Contains(peakMz))
                {
                    peakfindingIndices.Add(mz.Count);
                }

                mz.Add(peakMz);
                intensity.Add((float)intensityArray[i]);
                scanNumbers.Add(oneBasedScanNumber);
            }
        }

        /// <summary>
        /// Writes the window as a single JSON object: the bounds of the charge state, followed by its peaks as
        /// three parallel arrays. The peak level metadata is not repeated here, as the window is written as part of
        /// the <see cref="Ms1FeatureData"/> that carries it, and neither is anything a scan number stands for,
        /// which the file's <see cref="SpectraFileHeaderData"/> holds once per scan.
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
            writer.WriteNumber("scanStart", FirstScanNumber);
            writer.WriteNumber("scanEnd", LastScanNumber);
            writer.WriteNumber("scanCount", ScanCount);
            writer.WriteNumber("peakCount", PeakCount);

            WriteFloatArray(writer, "mz", Mz);
            WriteFloatArray(writer, "intensity", Intensity);
            WriteIntArray(writer, "scanNumbers", ScanNumbers);
            WriteIntArray(writer, "peakfindingIndices", PeakfindingIndices);

            writer.WriteEndObject();
        }

        private static void WriteFloatArray(Utf8JsonWriter writer, string propertyName, float[] values)
        {
            writer.WriteStartArray(propertyName);
            foreach (float value in values)
            {
                writer.WriteNumberValue(value);
            }
            writer.WriteEndArray();
        }

        private static void WriteIntArray(Utf8JsonWriter writer, string propertyName, int[] values)
        {
            writer.WriteStartArray(propertyName);
            foreach (int value in values)
            {
                writer.WriteNumberValue(value);
            }
            writer.WriteEndArray();
        }

        public override string ToString()
        {
            return "+" + ChargeState + "; " + PeakCount + " peaks in " + ScanCount + " scans; "
                + MinMz.ToString("F2") + "-" + MaxMz.ToString("F2") + " m/z; "
                + MinRetentionTime.ToString("F2") + "-" + MaxRetentionTime.ToString("F2") + " min";
        }
    }
}
