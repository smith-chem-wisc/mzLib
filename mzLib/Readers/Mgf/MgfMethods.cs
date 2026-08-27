using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MassSpectrometry;
using MzLibUtil;

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
    /// <see cref="Mgf"/> discards empty blocks on read. Writing a file with no blocks at all throws,
    /// because the reader cannot load one.
    /// </para>
    /// <para>
    /// MS1 scans are written by default, which is what lets a precursor be re-deconvoluted after a round
    /// trip. They cost portability: the Matrix Science specification requires a PEPMASS in every block and
    /// an MS1 has none, so such a file is out of spec. ProteoWizard and OpenMS read it anyway, but
    /// MSToolkit -- and therefore Comet -- calls exit(-12) on the first block without one. Pass
    /// includeMs1Scans: false to emit an MS2-only file that every reader accepts.
    /// </para>
    /// <para>
    /// Polarity does not survive a round trip for any scan without a charge state guess. MGF has no
    /// polarity field, so CHARGE is the only carrier of the sign, and CHARGE is written only where a
    /// guess exists. Such scans read back as positive. This applies to MS1 scans and equally to an MS2
    /// with a null charge guess, and is a limit of the format rather than of this writer.
    /// </para>
    /// </remarks>
    public static class MgfMethods
    {
        /// <param name="includeMs1Scans">
        /// True to write MS1 scans as their own blocks, preserving the precursor spectra a later search
        /// needs. False to write only MS2 and above, which keeps every block PEPMASS-bearing and therefore
        /// readable by tools that enforce the specification.
        /// </param>
        public static void WriteMgf(MsDataFile myMsDataFile, string outputFile, bool includeMs1Scans = true)
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

            List<MsDataScan> scansToWrite = new List<MsDataScan>();
            foreach (MsDataScan scan in myMsDataFile.GetAllScansList())
            {
                if (scan?.MassSpectrum != null && scan.MassSpectrum.Size > 0
                    && (includeMs1Scans || scan.MsnOrder > 1))
                {
                    scansToWrite.Add(scan);
                }
            }

            // Checked before opening the writer: a blockless mgf is not readable, so producing one would
            // hand the caller a file that throws on the way back in.
            if (scansToWrite.Count == 0)
            {
                throw new MzLibException(
                    "Cannot write an mgf with no spectra: every scan was empty or filtered out, and mgf "
                    + "requires at least one fragment peak per block.");
            }

            using StreamWriter output = new StreamWriter(outputFile);
            foreach (MsDataScan scan in scansToWrite)
            {
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

            // Four fields a search needs that the Matrix Science local parameters have no home for.
            // Extension keys rather than dropped data: the spec does not say what a reader must do with
            // an unknown key, but ProteoWizard (SpectrumList_MGF, "else { continue; // ignored
            // attribute }"), OpenMS (MascotGenericFile) and MSToolkit (MSReader) all skip them silently.
            // Without these a Calibrate then Search chain through mgf loses roughly a quarter of its
            // identifications.
            //
            // Written only when the source actually had the value -- an empty header would be worse than
            // an absent one, since a reader cannot tell "unknown" from "zero".

            // TIC belongs to every scan, not just the fragmented ones. It is not the sum of the peaks
            // written below: the recorded value predates centroiding and thresholding. MetaMorpheus
            // scores a match partly as a fraction of it, so recomputing it moves scores.
            if (!double.IsNaN(scan.TotalIonCurrent) && !double.IsInfinity(scan.TotalIonCurrent))
            {
                output.WriteLine($"TIC={scan.TotalIonCurrent.ToString(CultureInfo.InvariantCulture)}");
            }

            // The rest describe a precursor, so MS2 and above only.
            if (scan.MsnOrder > 1)
            {
                if (scan.DissociationType.HasValue && scan.DissociationType.Value != DissociationType.Unknown)
                {
                    output.WriteLine($"ACTIVATIONMETHOD={scan.DissociationType.Value}");
                }

                if (scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    output.WriteLine($"PRECURSORSCAN={scan.OneBasedPrecursorScanNumber.Value.ToString(CultureInfo.InvariantCulture)}");
                }

                if (scan.IsolationWidth.HasValue)
                {
                    output.WriteLine($"ISOLATIONWIDTH={scan.IsolationWidth.Value.ToString(CultureInfo.InvariantCulture)}");
                }

                if (scan.IsolationMz.HasValue)
                {
                    output.WriteLine($"ISOLATIONMZ={scan.IsolationMz.Value.ToString(CultureInfo.InvariantCulture)}");
                }
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
