using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;

namespace Readers
{
    /// <summary>
    /// Writes deconvolution results to the msAlign plain-text format.
    /// The format is compatible with TopFD, FlashDeconv, IsoDec, and the
    /// mzLib <see cref="Ms1Align"/> / <see cref="Ms2Align"/> readers.
    ///
    /// <para>
    /// Format rules enforced by this writer:
    /// <list type="bullet">
    ///   <item><description>
    ///     <c>RETENTION_TIME</c> is written in <b>seconds</b>.
    ///     <see cref="MsDataScan.RetentionTime"/> is in minutes; the writer
    ///     multiplies by 60.  The readers divide by 60, so values round-trip
    ///     correctly through any mzLib <c>MsAlign</c> subclass.
    ///   </description></item>
    ///   <item><description>
    ///     Peak lines are tab-separated:
    ///     <c>monoisotopicMass\tintensity\tcharge</c>.
    ///     Charges are always written as positive integers regardless of polarity.
    ///   </description></item>
    ///   <item><description>
    ///     Scans with no envelopes are silently skipped.
    ///   </description></item>
    /// </list>
    /// </para>
    /// </summary>
    public static class MsAlignWriter
    {
        // ── Public API ────────────────────────────────────────────────────────

        /// <summary>
        /// Writes MS1 deconvolution results to an ms1.msalign file.
        /// </summary>
        /// <param name="results">
        ///   Sequence of (scan, envelopes) pairs to write.
        ///   Pairs where <paramref name="results"/> envelopes are null or empty
        ///   are skipped; the scan is not written to the file.
        /// </param>
        /// <param name="outputPath">
        ///   Destination file path.  If the path does not already end with
        ///   <c>.msalign</c> the extension is appended automatically.
        /// </param>
        public static void WriteMs1Align(
            IEnumerable<(MsDataScan Scan, IEnumerable<IsotopicEnvelope> Envelopes)> results,
            string outputPath)
        {
            outputPath = EnsureExtension(outputPath);
            using var writer = new StreamWriter(outputPath);
            foreach (var (scan, envelopes) in results)
                WriteEntry(writer, scan, envelopes, msnOrder: 1);
        }

        /// <summary>
        /// Writes MS2 deconvolution results to an ms2.msalign file.
        /// </summary>
        /// <param name="results">
        ///   Sequence of (scan, envelopes) pairs to write.
        ///   Pairs where envelopes are null or empty are skipped.
        /// </param>
        /// <param name="outputPath">
        ///   Destination file path.  The <c>.msalign</c> extension is appended
        ///   if absent.
        /// </param>
        public static void WriteMs2Align(
            IEnumerable<(MsDataScan Scan, IEnumerable<IsotopicEnvelope> Envelopes)> results,
            string outputPath)
        {
            outputPath = EnsureExtension(outputPath);
            using var writer = new StreamWriter(outputPath);
            foreach (var (scan, envelopes) in results)
                WriteEntry(writer, scan, envelopes, msnOrder: 2);
        }

        // ── Private helpers ───────────────────────────────────────────────────

        private static void WriteEntry(
            StreamWriter writer,
            MsDataScan scan,
            IEnumerable<IsotopicEnvelope> envelopes,
            int msnOrder)
        {
            if (envelopes is null) return;

            var envList = envelopes
                .OrderBy(e => e.MonoisotopicMass)
                .ToList();

            if (envList.Count == 0) return;

            // ── header ────────────────────────────────────────────────────────
            writer.WriteLine("BEGIN IONS");
            writer.WriteLine($"ID={scan.OneBasedScanNumber}");
            writer.WriteLine($"FRACTION_ID=0");
            writer.WriteLine($"SCANS={scan.OneBasedScanNumber}");

            // RetentionTime is stored in minutes; msAlign uses seconds
            writer.WriteLine($"RETENTION_TIME={scan.RetentionTime * 60.0:F2}");
            writer.WriteLine($"LEVEL={msnOrder}");

            if (msnOrder == 2)
                WriteMs2Headers(writer, scan);

            // ── peak lines ────────────────────────────────────────────────────
            foreach (var env in envList)
            {
                // Charges are always positive in the msAlign format
                int positiveCharge = Math.Abs(env.Charge);
                writer.WriteLine(
                    $"{env.MonoisotopicMass:F6}\t{env.TotalIntensity:F2}\t{positiveCharge}");
            }

            writer.WriteLine("END IONS");
            writer.WriteLine(); // blank line between entries
        }

        private static void WriteMs2Headers(StreamWriter writer, MsDataScan scan)
        {
            if (scan.DissociationType.HasValue)
                writer.WriteLine($"ACTIVATION={scan.DissociationType.Value}");

            if (scan.OneBasedPrecursorScanNumber.HasValue)
            {
                writer.WriteLine($"MS_ONE_ID={scan.OneBasedPrecursorScanNumber.Value}");
                writer.WriteLine($"MS_ONE_SCAN={scan.OneBasedPrecursorScanNumber.Value}");
            }

            if (scan.SelectedIonMZ.HasValue)
                writer.WriteLine($"PRECURSOR_MZ={scan.SelectedIonMZ.Value:F6}");

            if (scan.SelectedIonChargeStateGuess.HasValue)
                writer.WriteLine(
                    $"PRECURSOR_CHARGE={Math.Abs(scan.SelectedIonChargeStateGuess.Value)}");

            if (scan.SelectedIonMonoisotopicGuessIntensity.HasValue)
                writer.WriteLine(
                    $"PRECURSOR_INTENSITY={scan.SelectedIonMonoisotopicGuessIntensity.Value:F2}");
        }

        private static string EnsureExtension(string path)
        {
            if (!path.EndsWith(".msalign", StringComparison.OrdinalIgnoreCase))
                return path + ".msalign";
            return path;
        }
    }
}