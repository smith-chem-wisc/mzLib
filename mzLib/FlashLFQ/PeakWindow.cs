using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    /// <summary>
    /// A rectangular slice of the MS1 data surrounding one charge state of a single ChromatographicPeak. The slice
    /// spans the scans that charge state was observed in, and runs from mzExpansion below the lowest m/z observed
    /// for it to mzExpansion above the highest. Every MS1 peak in that region is retained, including peaks that
    /// belong to species other than the quantified precursor, which is the point of the type: it captures the
    /// co-eluting and interfering signal that quantification discards.
    ///
    /// The window is per charge state rather than per peak because a ChromatographicPeak carries one isotopic
    /// envelope per scan per charge state, and the m/z of a peptide at one charge is nowhere near its m/z at
    /// another. Bounding all of them at once would produce a window hundreds of daltons wide spanning the gap
    /// between the charge states, whose contents are the rest of the spectrum rather than the peak's neighbours.
    ///
    /// The peaks are read straight out of the spectra file rather than out of a PeakIndexingEngine. The index is
    /// organized m/z-major and scan-minor, so a contiguous m/z slice of a handful of consecutive scans costs one
    /// binary search per m/z bin - 100 bins per dalton - where reading the scans costs one binary search per scan.
    /// </summary>
    public class PeakWindow
    {
        public const double DefaultMzExpansion = 0.5;

        public ChromatographicPeak Peak { get; }

        /// <summary>
        /// The charge state whose isotopic envelopes bound this window. A peak traced at more than one charge state
        /// produces one window per charge state, all of them describing the same peak.
        /// </summary>
        public int ChargeState { get; }

        public double MinMz { get; }
        public double MaxMz { get; }
        public double MinRetentionTime { get; }
        public double MaxRetentionTime { get; }
        public double ApexRetentionTime { get; }

        /// <summary>
        /// One entry per scan the peak was observed in, ordered by retention time. Scans in which nothing fell
        /// inside the m/z window are present but hold no peaks.
        /// </summary>
        public PeakWindowScan[] Scans { get; }

        /// <summary>
        /// The total number of MS1 peaks inside the window, across all of its scans.
        /// </summary>
        public int PeakCount { get; }

        private PeakWindow(ChromatographicPeak peak, int chargeState, double minMz, double maxMz,
            double apexRetentionTime, PeakWindowScan[] scans, int peakCount)
        {
            Peak = peak;
            ChargeState = chargeState;
            MinMz = minMz;
            MaxMz = maxMz;
            ApexRetentionTime = apexRetentionTime;
            Scans = scans;
            PeakCount = peakCount;
            MinRetentionTime = scans.First().RetentionTime;
            MaxRetentionTime = scans.Last().RetentionTime;
        }

        /// <summary>
        /// Maps the zero based scan indices carried by indexed peaks onto the one based scan numbers of the raw
        /// file. The two are not interchangeable: indexed peaks are numbered against an MS1-only array, so using a
        /// zero based index to address the data file lands on the wrong scan as soon as MS2 scans are interleaved.
        ///
        /// The MS1 scans are selected and ordered exactly as PeakIndexingEngine.InitializeIndexingEngine selects
        /// them, so that the indices recorded during quantification line up with this map.
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
        /// Extracts one window per charge state the peak was traced at, in ascending charge order. Peaks with no
        /// isotopic envelopes yield nothing, as those have no observed m/z or scan to build a window around.
        /// </summary>
        /// <param name="peak"> the quantified peak to build windows around </param>
        /// <param name="dataFile"> the loaded spectra file the peak was quantified from </param>
        /// <param name="ms1ScanNumberMap"> the map built by BuildMs1ScanNumberMap for that same file </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z each window extends </param>
        public static IEnumerable<PeakWindow> CreateForEachChargeState(ChromatographicPeak peak, MsDataFile dataFile,
            int[] ms1ScanNumberMap, double mzExpansion = DefaultMzExpansion)
        {
            if (peak == null || dataFile == null || ms1ScanNumberMap == null)
                yield break;

            foreach (int chargeState in peak.IsotopicEnvelopes.Select(e => e.ChargeState).Distinct().OrderBy(z => z))
            {
                PeakWindow window = Create(peak, dataFile, ms1ScanNumberMap, chargeState, mzExpansion);
                if (window != null)
                    yield return window;
            }
        }

        /// <summary>
        /// Extracts the window of MS1 data surrounding one charge state of a chromatographic peak. Returns null when
        /// the peak was not traced at that charge state, as there is then no observed m/z or scan to build a window
        /// around, and when the scans it was traced in fall outside the supplied map.
        /// </summary>
        /// <param name="peak"> the quantified peak to build a window around </param>
        /// <param name="dataFile"> the loaded spectra file the peak was quantified from </param>
        /// <param name="ms1ScanNumberMap"> the map built by BuildMs1ScanNumberMap for that same file </param>
        /// <param name="chargeState"> the charge state whose isotopic envelopes bound the window </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z the window extends </param>
        public static PeakWindow Create(ChromatographicPeak peak, MsDataFile dataFile, int[] ms1ScanNumberMap,
            int chargeState, double mzExpansion = DefaultMzExpansion)
        {
            if (peak == null || dataFile == null || ms1ScanNumberMap == null)
                return null;

            // Only this charge state's envelopes bound the window. A peptide's m/z at one charge is unrelated to its
            // m/z at another, so bounding every envelope at once would span the gap between the charge states.
            List<IsotopicEnvelope> envelopes = peak.IsotopicEnvelopes.Where(e => e.ChargeState == chargeState).ToList();
            if (envelopes.Count == 0)
                return null;

            double minMz = envelopes.Min(e => e.IndexedPeak.M) - mzExpansion;
            double maxMz = envelopes.Max(e => e.IndexedPeak.M) + mzExpansion;

            int firstMs1Index = envelopes.Min(e => e.IndexedPeak.ZeroBasedScanIndex);
            int lastMs1Index = envelopes.Max(e => e.IndexedPeak.ZeroBasedScanIndex);
            if (firstMs1Index < 0 || lastMs1Index >= ms1ScanNumberMap.Length)
                return null;

            // The m/z values FlashLFQ traced, grouped by the scan they were found in, so that the peaks it actually
            // used can be flagged as the window is written out.
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

            var scans = new PeakWindowScan[lastMs1Index - firstMs1Index + 1];
            int peakCount = 0;
            for (int ms1Index = firstMs1Index; ms1Index <= lastMs1Index; ms1Index++)
            {
                int oneBasedScanNumber = ms1ScanNumberMap[ms1Index];
                MsDataScan scan = dataFile.GetOneBasedScan(oneBasedScanNumber);
                peakfindingMzByScanNumber.TryGetValue(oneBasedScanNumber, out HashSet<float> peakfindingMz);

                PeakWindowScan windowScan = PeakWindowScan.Create(scan, minMz, maxMz, peakfindingMz);
                scans[ms1Index - firstMs1Index] = windowScan;
                peakCount += windowScan.Mz.Length;
            }

            // The apex of this charge state, not of the peak as a whole: the peak's apex can belong to a different
            // charge state, whose scan is not necessarily inside this window.
            double apexRetentionTime = double.NaN;
            int apexMs1Index = envelopes.MaxBy(e => e.Intensity).IndexedPeak.ZeroBasedScanIndex;
            if (apexMs1Index >= 0 && apexMs1Index < ms1ScanNumberMap.Length)
            {
                apexRetentionTime = dataFile.GetOneBasedScan(ms1ScanNumberMap[apexMs1Index]).RetentionTime;
            }

            return new PeakWindow(peak, chargeState, minMz, maxMz, apexRetentionTime, scans, peakCount);
        }

        public static string TabSeparatedHeader => string.Join("\t",
            "File Name",
            "Peak ID",
            "Base Sequence",
            "Full Sequence",
            "Peak Detection Type",
            "Peak Charge",
            "Window RT Start",
            "Window RT Apex",
            "Window RT End",
            "Window Min MZ",
            "Window Max MZ",
            "Scan Number",
            "Scan Retention Time",
            "MZ",
            "Intensity",
            "Is Peakfinding Peak");

        /// <summary>
        /// One tab separated line per MS1 peak in the window. Rows are yielded rather than returned as a list so
        /// that a caller writing many windows never holds more than one window's worth of text at a time.
        /// </summary>
        /// <param name="peakId"> an identifier distinguishing this window's peak from the file's other peaks. All of
        /// a peak's charge states share it, and are told apart by the Peak Charge column. </param>
        public IEnumerable<string> ToTsvRows(int peakId)
        {
            // Everything to the left of the scan columns is constant across the window, so it is built once.
            StringBuilder peakColumns = new StringBuilder();
            peakColumns.Append(Peak.SpectraFileInfo.FilenameWithoutExtension + "\t");
            peakColumns.Append(peakId.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(string.Join("|", Peak.Identifications.Select(p => p.BaseSequence).Distinct()) + "\t");
            peakColumns.Append(string.Join("|", Peak.Identifications.Select(p => p.ModifiedSequence).Distinct()) + "\t");
            peakColumns.Append(Peak.DetectionType.ToString() + "\t");
            peakColumns.Append(ChargeState.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(MinRetentionTime.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(ApexRetentionTime.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(MaxRetentionTime.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(MinMz.ToString(CultureInfo.InvariantCulture) + "\t");
            peakColumns.Append(MaxMz.ToString(CultureInfo.InvariantCulture) + "\t");
            string peakColumnString = peakColumns.ToString();

            foreach (PeakWindowScan scan in Scans)
            {
                string scanColumns = scan.OneBasedScanNumber.ToString(CultureInfo.InvariantCulture) + "\t"
                    + scan.RetentionTime.ToString(CultureInfo.InvariantCulture) + "\t";

                for (int i = 0; i < scan.Mz.Length; i++)
                {
                    StringBuilder sb = new StringBuilder(peakColumnString);
                    sb.Append(scanColumns);
                    sb.Append(scan.Mz[i].ToString(CultureInfo.InvariantCulture) + "\t");
                    sb.Append(scan.Intensity[i].ToString(CultureInfo.InvariantCulture) + "\t");
                    sb.Append(scan.IsPeakfindingPeak[i].ToString(CultureInfo.InvariantCulture));
                    yield return sb.ToString();
                }
            }
        }

        public override string ToString()
        {
            return Peak.Identifications.First().ModifiedSequence + "; z=" + ChargeState + "; " + PeakCount
                + " peaks in " + Scans.Length + " scans; " + MinMz.ToString("F2") + "-" + MaxMz.ToString("F2") + " m/z; "
                + MinRetentionTime.ToString("F2") + "-" + MaxRetentionTime.ToString("F2") + " min";
        }
    }

    /// <summary>
    /// The portion of one MS1 scan that falls inside a PeakWindow, held as parallel arrays rather than as a
    /// collection of peak objects. A window can cover thousands of peaks, and this output exists to be written
    /// straight to disk, so there is nothing to gain from boxing each one into an object.
    /// </summary>
    public class PeakWindowScan
    {
        public int OneBasedScanNumber { get; }
        public double RetentionTime { get; }
        public double[] Mz { get; }
        public double[] Intensity { get; }

        /// <summary>
        /// Parallel to Mz and Intensity. Flags the peaks FlashLFQ used to trace the feature: the most abundant
        /// isotope peak of each isotopic envelope, one per scan per charge state. The remaining isotope peaks of
        /// each envelope are not retained by ChromatographicPeak, so they appear in the window unflagged.
        /// </summary>
        public bool[] IsPeakfindingPeak { get; }

        private PeakWindowScan(int oneBasedScanNumber, double retentionTime, double[] mz, double[] intensity,
            bool[] isPeakfindingPeak)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            RetentionTime = retentionTime;
            Mz = mz;
            Intensity = intensity;
            IsPeakfindingPeak = isPeakfindingPeak;
        }

        /// <summary>
        /// Copies the peaks between minMz and maxMz, inclusive, out of a scan. The scan's m/z array is sorted
        /// ascending, so the slice is found with one binary search and copied out in bulk.
        /// </summary>
        /// <param name="peakfindingMz"> m/z values FlashLFQ traced in this scan, or null if it traced none </param>
        internal static PeakWindowScan Create(MsDataScan scan, double minMz, double maxMz, HashSet<float> peakfindingMz)
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
                return new PeakWindowScan(scan.OneBasedScanNumber, scan.RetentionTime, Array.Empty<double>(),
                    Array.Empty<double>(), Array.Empty<bool>());
            }

            double[] mz = new double[count];
            double[] intensity = new double[count];
            Array.Copy(mzArray, start, mz, 0, count);
            Array.Copy(intensityArray, start, intensity, 0, count);

            bool[] isPeakfindingPeak = new bool[count];
            if (peakfindingMz != null)
            {
                // Indexed peaks store m/z as a float narrowed from this same array, so comparing the narrowed value
                // recovers an exact match.
                for (int i = 0; i < count; i++)
                {
                    isPeakfindingPeak[i] = peakfindingMz.Contains((float)mz[i]);
                }
            }

            return new PeakWindowScan(scan.OneBasedScanNumber, scan.RetentionTime, mz, intensity, isPeakfindingPeak);
        }
    }
}
