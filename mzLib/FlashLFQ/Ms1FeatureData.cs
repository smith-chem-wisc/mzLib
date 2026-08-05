using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;

namespace FlashLFQ
{
    /// <summary>
    /// The raw MS1 data underlying one quantified ChromatographicPeak, held as one <see cref="PeakWindowData"/> per
    /// charge state the peak was observed at. A window is a rectangle in m/z and retention time, and the charge
    /// states of a peptide sit at unrelated m/z values, so one rectangle cannot cover them all without swallowing
    /// the gaps between them. This type is the join: the charge states are extracted separately and recombined
    /// here, under the peak they were all traced for.
    ///
    /// This class is designed to hold rich data representations of the MS1 signal produced by peptides or proteoforms
    /// as they elute from the LC column, and to write them out in a compact JSON representation. This data can be
    /// used for visualization, re-analysis, or training machine learning models for peptide detection and quantification.
    ///
    /// Named Ms1FeatureData rather than Ms1Feature because Readers.Ms1Feature already names a row of an
    /// external .ms1.feature file, and the two types are routinely in scope together.
    /// </summary>
    public class Ms1FeatureData
    {
        public ChromatographicPeak Peak { get; }

        /// <summary>
        /// One window per charge state the peak was observed at, ascending by charge state. A charge state whose
        /// envelopes could not be placed in the spectra file is dropped, so this can hold fewer charge states than
        /// the peak observed - but never zero, as a feature with no windows is not constructed at all.
        /// </summary>
        public IReadOnlyList<PeakWindowData> Windows { get; }

        public SpectraFileInfo SpectraFileInfo => Peak.SpectraFileInfo;

        /// <summary>
        /// The charge state of the peak's most intense isotopic envelope, or 0 for a peak whose intensity has not
        /// been calculated and so has no apex.
        /// </summary>
        public int ApexChargeState => Peak.Apex?.ChargeState ?? 0;

        /// <summary>
        /// The retention time of the peak's most intense isotopic envelope, across all of its charge states, or
        /// NaN for a peak with no apex.
        /// </summary>
        public double ApexRetentionTime { get; }

        public IEnumerable<int> ChargeStates => Windows.Select(w => w.ChargeState);

        public float MinMz => Windows.MaxBy(w => w.ChargeState).MinMz;
        public float MaxMz => Windows.MinBy(w => w.ChargeState).MaxMz;
        public double MinRetentionTime => Windows.Min(w => w.MinRetentionTime);
        public double MaxRetentionTime => Windows.Max(w => w.MaxRetentionTime);

        /// <summary>
        /// The total number of MS1 peaks across every window. Two charge states of the same peak can in principle
        /// land close enough in m/z for their windows to overlap, in which case the shared peaks are counted - and
        /// written - once per window.
        /// </summary>
        public int PeakCount => Windows.Sum(w => w.PeakCount);

        private Ms1FeatureData(ChromatographicPeak peak, IReadOnlyList<PeakWindowData> windows,
            double apexRetentionTime)
        {
            Peak = peak;
            Windows = windows;
            ApexRetentionTime = apexRetentionTime;
        }

        /// <summary>
        /// The window covering a given charge state, or null if the peak was not observed at it.
        /// </summary>
        public PeakWindowData GetWindow(int chargeState)
        {
            return Windows.FirstOrDefault(w => w.ChargeState == chargeState);
        }

        /// <summary>
        /// The MS1 scans Create would read for a peak, as one inclusive range of zero based scan indices per charge
        /// state: everything between that charge state's first and last isotopic envelope. Ranges from different
        /// charge states routinely overlap, since the charge states of a peptide elute together.
        ///
        /// This is what a peak covers, and it is derived from the envelopes alone, so a caller can size the work
        /// or collect the scans involved without reading anything. Both of them do, which is why it lives here
        /// rather than being spelled out again at each call site.
        /// </summary>
        public static IEnumerable<(int FirstScanIndex, int LastScanIndex)> GetScanIndexRanges(ChromatographicPeak peak)
        {
            if (peak == null)
                return Enumerable.Empty<(int, int)>();

            return peak.IsotopicEnvelopes
                .GroupBy(e => e.ChargeState)
                .Select(g => (g.Min(e => e.IndexedPeak.ZeroBasedScanIndex),
                    g.Max(e => e.IndexedPeak.ZeroBasedScanIndex)));
        }

        /// <summary>
        /// How many scans Create would read for a peak: every scan of every charge state's window, counting a scan
        /// once per window that covers it, since that is once per read. Lets a caller decide whether seeking to
        /// those scans is cheaper than loading the file they are in, without reading anything to find out.
        /// </summary>
        public static int CountScansToRead(ChromatographicPeak peak)
        {
            return GetScanIndexRanges(peak).Sum(r => r.LastScanIndex - r.FirstScanIndex + 1);
        }

        /// <summary>
        /// The MS1 scans of a file, in the order the zero based indices carried by indexed peaks number them, by
        /// asking PeakIndexingEngine which scans it would have indexed. Going through the engine is what makes the
        /// numbering correct: the indices being translated were assigned by that same selection, and a second copy
        /// of it here would line up only until one of them changed.
        ///
        /// This reads the file, so it is the fallback for callers with no run behind them. FlashLfqResults prefers
        /// the ScanInfo the indexing engine already built during quantification.
        /// </summary>
        public static ScanInfo[] BuildMs1ScanInfo(MsDataFile dataFile)
        {
            MsDataScan[] scans = PeakIndexingEngine.GetScansToIndex(dataFile);
            var scanInfo = new ScanInfo[scans.Length];
            for (int i = 0; i < scans.Length; i++)
            {
                scanInfo[i] = new ScanInfo(scans[i].OneBasedScanNumber, i, scans[i].RetentionTime,
                    scans[i].MsnOrder, scans[i].TotalIonCurrent);
            }
            return scanInfo;
        }

        /// <summary>
        /// Extracts the MS1 data surrounding a chromatographic peak, one window per charge state. Returns null for
        /// peaks that produced no windows at all: peaks with no isotopic envelopes, and peaks whose every charge
        /// state sits in scans the scan metadata does not cover.
        /// </summary>
        /// <param name="peak"> the quantified peak to build windows around </param>
        /// <param name="dataFile"> the spectra file the peak was quantified from </param>
        /// <param name="ms1ScanInfo"> the MS1 scans of that same file, as built by BuildMs1ScanInfo or left behind
        /// by the run that quantified the peak </param>
        /// <param name="mzExpansion"> how far below the lowest and above the highest observed m/z each window extends </param>
        /// <param name="useDynamicConnection"> read the scans through a dynamic connection the caller has already
        /// opened, rather than out of a statically loaded file </param>
        public static Ms1FeatureData Create(ChromatographicPeak peak, MsDataFile dataFile,
            ScanInfo[] ms1ScanInfo, double mzExpansion = PeakWindowData.DefaultMzExpansion,
            bool useDynamicConnection = false)
        {
            if (peak == null || dataFile == null || ms1ScanInfo == null || !peak.IsotopicEnvelopes.Any())
                return null;

            var windows = new List<PeakWindowData>();
            foreach (var envelopesForCharge in peak.IsotopicEnvelopes
                         .GroupBy(e => e.ChargeState)
                         .OrderBy(g => g.Key))
            {
                PeakWindowData window = PeakWindowData.Create(envelopesForCharge.ToList(), dataFile, ms1ScanInfo,
                    mzExpansion, useDynamicConnection);
                if (window != null)
                {
                    windows.Add(window);
                }
            }

            if (windows.Count == 0)
                return null;

            // The peak's apex is the apex of one of its charge states, and that charge state's window has already
            // resolved the retention time, so no scan is read a second time to find it here.
            double apexRetentionTime = double.NaN;
            if (peak.Apex != null)
            {
                PeakWindowData apexWindow = windows.FirstOrDefault(w => w.ChargeState == peak.Apex.ChargeState);
                if (apexWindow != null)
                {
                    apexRetentionTime = apexWindow.ApexRetentionTime;
                }
            }

            return new Ms1FeatureData(peak, windows, apexRetentionTime);
        }

        /// <summary>
        /// Writes the feature as a single JSON object. The peak's metadata is written once, followed by its
        /// windows, so the per-peak overhead does not scale with the number of MS1 peaks the windows cover.
        /// </summary>
        /// <param name="writer"> the writer to emit the object into </param>
        /// <param name="featureId"> an identifier distinguishing this feature from the file's other features </param>
        public void WriteTo(Utf8JsonWriter writer, int featureId)
        {
            writer.WriteStartObject();

            writer.WriteString("type", "feature");
            writer.WriteString("fileName", SpectraFileInfo.FilenameWithoutExtension);
            writer.WriteNumber("featureId", featureId);
            writer.WriteString("baseSequence", string.Join("|", Peak.Identifications.Select(p => p.BaseSequence).Distinct()));
            writer.WriteString("fullSequence", string.Join("|", Peak.Identifications.Select(p => p.ModifiedSequence).Distinct()));
            writer.WriteString("detectionType", Peak.DetectionType.ToString());
            writer.WriteNumber("apexCharge", ApexChargeState);
            writer.WriteNumber("rtStart", MinRetentionTime);
            WriteNumberOrNull(writer, "rtApex", ApexRetentionTime);
            writer.WriteNumber("rtEnd", MaxRetentionTime);
            writer.WriteNumber("minMz", MinMz);
            writer.WriteNumber("maxMz", MaxMz);
            writer.WriteNumber("peakCount", PeakCount);

            writer.WriteStartArray("charges");
            foreach (int chargeState in ChargeStates)
            {
                writer.WriteNumberValue(chargeState);
            }
            writer.WriteEndArray();

            writer.WriteStartArray("windows");
            foreach (PeakWindowData window in Windows)
            {
                window.WriteTo(writer);
            }
            writer.WriteEndArray();

            writer.WriteEndObject();
        }

        /// <summary>
        /// NaN is not representable in JSON, so a peak with no apex is written as null rather than as a value the
        /// reader would have to reject.
        /// </summary>
        internal static void WriteNumberOrNull(Utf8JsonWriter writer, string propertyName, double value)
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
            return Peak.Identifications.First().ModifiedSequence + "; "
                + Windows.Count + " charge states (" + string.Join(", ", ChargeStates.Select(z => "+" + z)) + "); "
                + PeakCount + " peaks; " + MinMz.ToString("F2") + "-" + MaxMz.ToString("F2") + " m/z; "
                + MinRetentionTime.ToString("F2") + "-" + MaxRetentionTime.ToString("F2") + " min";
        }
    }
}
