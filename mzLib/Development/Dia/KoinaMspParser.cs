// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry.Dia;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace Development.Dia
{
    /// <summary>
    /// Lightweight parser for Koina/Prosit .msp spectral library files.
    /// 
    /// Koina .msp format (from real output):
    ///   Name: PEPTIDEK/2
    ///   MW: 1184.608784670464              (NEUTRAL monoisotopic mass, NOT m/z)
    ///   Comment: Parent=1184.60878... RT=67.30357...
    ///   Num peaks: 64
    ///   147.076  0.05209  "y1^1/0ppm"
    ///   ...
    ///   Name: NEXTPEPTIDE/3                (next entry starts immediately, NO blank line)
    /// 
    /// Two critical bugs fixed from v1:
    ///   1. Entries are delimited by the next "Name:" line, NOT blank lines.
    ///   2. MW/Parent= is the neutral mass. m/z = (MW + z * 1.00728) / z.
    /// </summary>
    public static class KoinaMspParser
    {
        public static List<LibraryPrecursorInput> Parse(
            string mspPath,
            Dictionary<string, double> retentionTimeLookup = null,
            float minIntensity = 0.05f)
        {
            if (!File.Exists(mspPath))
                throw new FileNotFoundException($"MSP file not found: {mspPath}");

            var results = new List<LibraryPrecursorInput>();
            int entriesParsed = 0;
            int entriesSkipped = 0;
            int totalFragments = 0;

            // Current entry state
            string currentSequence = null;
            int currentCharge = 0;
            double currentMW = 0;
            double currentMspRt = -1;
            var mzList = new List<float>(128);
            var intList = new List<float>(128);
            bool inPeaks = false;

            using var reader = new StreamReader(mspPath);
            string line;

            while ((line = reader.ReadLine()) != null)
            {
                string trimmed = line.Trim();
                if (string.IsNullOrEmpty(trimmed)) continue;

                // A new "Name:" line means: finalize previous entry, start new one
                if (trimmed.StartsWith("Name:", StringComparison.OrdinalIgnoreCase))
                {
                    // Finalize previous entry
                    if (currentSequence != null && mzList.Count > 0)
                    {
                        var input = BuildPrecursorInput(
                            currentSequence, currentCharge, currentMW, currentMspRt,
                            mzList, intList, minIntensity, retentionTimeLookup);
                        if (input.HasValue)
                        {
                            results.Add(input.Value);
                            totalFragments += input.Value.FragmentCount;
                            entriesParsed++;
                        }
                        else entriesSkipped++;
                    }

                    // Reset for new entry
                    currentMW = 0;
                    currentMspRt = -1;
                    mzList.Clear();
                    intList.Clear();
                    inPeaks = false;
                    ParseNameLine(trimmed, out currentSequence, out currentCharge);
                    continue;
                }

                if (trimmed.StartsWith("MW:", StringComparison.OrdinalIgnoreCase))
                {
                    string val = trimmed.Substring(trimmed.IndexOf(':') + 1).Trim();
                    double.TryParse(val, NumberStyles.Float, CultureInfo.InvariantCulture, out currentMW);
                    continue;
                }

                if (trimmed.StartsWith("Comment:", StringComparison.OrdinalIgnoreCase))
                {
                    currentMspRt = ParseFieldFromComment(trimmed, "RT=");
                    if (currentMW <= 0)
                        currentMW = ParseFieldFromComment(trimmed, "Parent=");
                    continue;
                }

                if (trimmed.StartsWith("Num peaks:", StringComparison.OrdinalIgnoreCase) ||
                    trimmed.StartsWith("Num Peaks:", StringComparison.OrdinalIgnoreCase))
                {
                    inPeaks = true;
                    continue;
                }

                // Peak lines start with digit or minus
                if (inPeaks && trimmed.Length > 0 && (char.IsDigit(trimmed[0]) || trimmed[0] == '-'))
                {
                    ParsePeakLine(trimmed, mzList, intList);
                }
            }

            // Finalize last entry
            if (currentSequence != null && mzList.Count > 0)
            {
                var input = BuildPrecursorInput(
                    currentSequence, currentCharge, currentMW, currentMspRt,
                    mzList, intList, minIntensity, retentionTimeLookup);
                if (input.HasValue)
                {
                    results.Add(input.Value);
                    totalFragments += input.Value.FragmentCount;
                    entriesParsed++;
                }
                else entriesSkipped++;
            }

            Console.WriteLine($"  MSP parsed: {entriesParsed:N0} entries, {totalFragments:N0} total fragments" +
                              (entriesSkipped > 0 ? $", {entriesSkipped} skipped" : ""));

            return results;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Line parsers
        // ═══════════════════════════════════════════════════════════════

        private static void ParseNameLine(string line, out string sequence, out int charge)
        {
            sequence = null;
            charge = 0;
            string value = line.Substring(line.IndexOf(':') + 1).Trim();
            int slashIdx = value.LastIndexOf('/');
            if (slashIdx > 0 && slashIdx < value.Length - 1 &&
                int.TryParse(value.Substring(slashIdx + 1), out charge))
            {
                string seqPart = value.Substring(0, slashIdx);
                sequence = System.Text.RegularExpressions.Regex.Replace(
                    seqPart, @"\[UNIMOD:(\d+)\]", "(UniMod:$1)",
                    System.Text.RegularExpressions.RegexOptions.IgnoreCase);
            }
        }

        /// <summary>
        /// Extracts a numeric field (e.g. "RT=" or "Parent=") from a Comment line.
        /// Returns -1 if not found.
        /// </summary>
        private static double ParseFieldFromComment(string line, string fieldName)
        {
            string comment = line.Substring(line.IndexOf(':') + 1);
            int idx = comment.IndexOf(fieldName, StringComparison.OrdinalIgnoreCase);
            if (idx < 0) return -1;

            int start = idx + fieldName.Length;
            int end = start;
            while (end < comment.Length && (char.IsDigit(comment[end]) || comment[end] == '.' || comment[end] == '-'))
                end++;

            if (end > start &&
                double.TryParse(comment.Substring(start, end - start),
                    NumberStyles.Float, CultureInfo.InvariantCulture, out double val))
                return val;
            return -1;
        }

        private static void ParsePeakLine(string line, List<float> mzList, List<float> intList)
        {
            // Format: "147.076416015625\t0.05209549888968468\t\"y1^1/0ppm\""
            var parts = line.Split(new[] { '\t', ' ' }, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length >= 2 &&
                float.TryParse(parts[0], NumberStyles.Float, CultureInfo.InvariantCulture, out float mz) &&
                float.TryParse(parts[1], NumberStyles.Float, CultureInfo.InvariantCulture, out float intensity))
            {
                mzList.Add(mz);
                intList.Add(intensity);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Input builder
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// In standard NIST .msp format (including Koina output), the MW and Parent=
        /// fields contain the PRECURSOR M/Z, not the neutral peptide mass.
        /// 
        /// Evidence: for AAAAAAAAAPAAAATAPTTAATTAATAAQ/2 (27 residues, neutral mass ~2367 Da),
        /// MW = 1184.6 ≈ (2367 + 2*1.007)/2, which is the [M+2H]²⁺ m/z.
        /// </summary>
        private static LibraryPrecursorInput? BuildPrecursorInput(
            string sequence, int charge, double mwOrMz, double mspRt,
            List<float> mzList, List<float> intList,
            float minIntensity,
            Dictionary<string, double> rtLookup)
        {
            if (sequence == null || charge <= 0 || mzList.Count == 0) return null;

            // MW/Parent in MSP is already the precursor m/z
            double precursorMz = mwOrMz;
            if (precursorMz <= 0) return null;

            // Filter by relative intensity
            float maxInt = 0;
            for (int i = 0; i < intList.Count; i++)
                if (intList[i] > maxInt) maxInt = intList[i];

            float threshold = maxInt * minIntensity;
            var filtMz = new List<float>();
            var filtInt = new List<float>();
            for (int i = 0; i < mzList.Count; i++)
            {
                if (intList[i] > threshold && intList[i] > 0)
                {
                    filtMz.Add(mzList[i]);
                    filtInt.Add(intList[i]);
                }
            }
            if (filtMz.Count == 0) return null;

            var mzArr = filtMz.ToArray();
            var intArr = filtInt.ToArray();
            Array.Sort(mzArr, intArr);

            // RT: prefer DIA-NN ground truth, fall back to MSP Comment RT
            string key = sequence + "/" + charge;
            double? rt = null;
            if (rtLookup != null && rtLookup.TryGetValue(key, out double rtVal))
                rt = rtVal;
            else if (mspRt >= 0)
                rt = mspRt;

            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: precursorMz,
                chargeState: charge,
                retentionTime: rt,
                isDecoy: false,
                fragmentMzs: mzArr,
                fragmentIntensities: intArr);
        }

        // ═══════════════════════════════════════════════════════════════
        //  RT lookup from DIA-NN TSV
        // ═══════════════════════════════════════════════════════════════

        public static Dictionary<string, double> BuildRtLookupFromDiannTsv(string tsvPath)
        {
            var lookup = new Dictionary<string, double>();
            var lines = File.ReadAllLines(tsvPath);
            if (lines.Length < 2) return lookup;

            var header = lines[0].Split('\t');
            var cols = new Dictionary<string, int>();
            for (int i = 0; i < header.Length; i++) cols[header[i].Trim()] = i;

            int seqCol = FindCol(cols, "Modified.Sequence", "Stripped.Sequence", "Sequence");
            int chargeCol = FindCol(cols, "Precursor.Charge", "Charge");
            int rtCol = FindCol(cols, "RT", "RetentionTime", "iRT");

            if (seqCol < 0 || chargeCol < 0 || rtCol < 0)
            {
                Console.WriteLine($"  WARNING: RT lookup columns not found. Available: {string.Join(", ", cols.Keys)}");
                return lookup;
            }

            for (int row = 1; row < lines.Length; row++)
            {
                var f = lines[row].Split('\t');
                if (f.Length <= Math.Max(seqCol, Math.Max(chargeCol, rtCol))) continue;

                string seq = f[seqCol].Trim();
                if (!int.TryParse(f[chargeCol], out int charge)) continue;
                if (!double.TryParse(f[rtCol], NumberStyles.Float,
                    CultureInfo.InvariantCulture, out double rt)) continue;

                string key = seq + "/" + charge;
                if (!lookup.ContainsKey(key))
                    lookup[key] = rt;
            }

            Console.WriteLine($"  RT lookup: {lookup.Count:N0} entries from {Path.GetFileName(tsvPath)}");
            return lookup;
        }

        private static int FindCol(Dictionary<string, int> cols, params string[] names)
        {
            foreach (var name in names)
                if (cols.TryGetValue(name, out int idx)) return idx;
            return -1;
        }
    }
}