using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MassSpectrometry;

namespace Readers
{
    /// <summary>
    /// Writes an <see cref="MsDataFile"/> as Mascot Generic Format. Mirrors <see cref="MzmlMethods"/>,
    /// the equivalent entry point for mzML.
    /// </summary>
    /// <remarks>
    /// Format rules follow the Matrix Science specification: no whitespace around '=', all parameters
    /// precede the fragment peaks within a block, and numbers use the invariant culture so the decimal
    /// separator is always '.'.
    /// <para>
    /// MSLEVEL is not part of the core Mascot specification but is understood widely enough (SIRIUS,
    /// pyteomics, MsBackendMgf) to be the standard way to record MS level, and it is what lets MS1 and
    /// MS3 scans survive a round trip through this format.
    /// </para>
    /// <para>
    /// Scans with no peaks are skipped: MGF requires at least one fragment per block, and
    /// <see cref="Mgf"/> discards empty blocks on read.
    /// </para>
    /// </remarks>
    public static class MgfMethods
    {
        public static void WriteMgf(MsDataFile myMsDataFile, string outputFile)
        {
            if (myMsDataFile == null)
            {
                throw new ArgumentNullException(nameof(myMsDataFile));
            }

            if (myMsDataFile.Scans == null)
            {
                myMsDataFile.LoadAllStaticData();
            }

            string sourceName = string.IsNullOrEmpty(myMsDataFile.FilePath)
                ? Path.GetFileNameWithoutExtension(outputFile)
                : Path.GetFileNameWithoutExtension(myMsDataFile.FilePath);

            using StreamWriter output = new StreamWriter(outputFile);
            foreach (MsDataScan scan in myMsDataFile.GetAllScansList())
            {
                if (scan?.MassSpectrum == null || scan.MassSpectrum.Size == 0)
                {
                    continue;
                }

                WriteScan(output, scan, sourceName);
            }
        }

        private static void WriteScan(StreamWriter output, MsDataScan scan, string sourceName)
        {
            output.WriteLine("BEGIN IONS");
            output.WriteLine($"TITLE={BuildTitle(scan, sourceName)}");
            output.WriteLine($"MSLEVEL={scan.MsnOrder.ToString(CultureInfo.InvariantCulture)}");

            // MS1 has no precursor, so PEPMASS and CHARGE are omitted rather than written empty
            double? precursorMz = scan.SelectedIonMonoisotopicGuessMz ?? scan.SelectedIonMZ;
            if (precursorMz.HasValue)
            {
                string pepmass = precursorMz.Value.ToString(CultureInfo.InvariantCulture);
                if (scan.SelectedIonIntensity.HasValue)
                {
                    pepmass += " " + scan.SelectedIonIntensity.Value.ToString(CultureInfo.InvariantCulture);
                }
                output.WriteLine($"PEPMASS={pepmass}");
            }

            if (scan.SelectedIonChargeStateGuess.HasValue)
            {
                int magnitude = Math.Abs(scan.SelectedIonChargeStateGuess.Value);
                string sign = scan.Polarity == Polarity.Negative ? "-" : "+";
                output.WriteLine($"CHARGE={magnitude.ToString(CultureInfo.InvariantCulture)}{sign}");
            }

            output.WriteLine($"SCANS={scan.OneBasedScanNumber.ToString(CultureInfo.InvariantCulture)}");

            if (!double.IsNaN(scan.RetentionTime) && !double.IsInfinity(scan.RetentionTime))
            {
                double rtInSeconds = scan.RetentionTime * 60.0;
                output.WriteLine($"RTINSECONDS={rtInSeconds.ToString(CultureInfo.InvariantCulture)}");
            }

            double[] mzs = scan.MassSpectrum.XArray;
            double[] intensities = scan.MassSpectrum.YArray;
            for (int i = 0; i < mzs.Length; i++)
            {
                output.WriteLine($"{mzs[i].ToString(CultureInfo.InvariantCulture)} {intensities[i].ToString(CultureInfo.InvariantCulture)}");
            }

            output.WriteLine("END IONS");
        }

        /// <summary>
        /// Follows the widely parsed "basename.scan.scan.charge" convention, so tools that recover scan
        /// identity from TITLE rather than SCANS still work. Mgf itself ignores this line.
        /// </summary>
        private static string BuildTitle(MsDataScan scan, string sourceName)
        {
            int charge = scan.SelectedIonChargeStateGuess.HasValue
                ? Math.Abs(scan.SelectedIonChargeStateGuess.Value)
                : 0;

            return $"{sourceName}.{scan.OneBasedScanNumber}.{scan.OneBasedScanNumber}.{charge} " +
                   $"File:\"{sourceName}\", NativeID:\"{scan.NativeId ?? string.Empty}\"";
        }
    }
}
