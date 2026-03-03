// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace Development.Dia
{
    /// <summary>
    /// Generates a Koina/Prosit input table from a DIA-NN ground truth TSV.
    /// 
    /// Usage (add to Program.cs):
    ///     KoinaTableGenerator.Generate(
    ///         @"F:\DiaBenchmark\PXD005573\diann_output.tsv",
    ///         @"F:\DiaBenchmark\PXD005573\koina_input.tsv",
    ///         collisionEnergy: 27);
    /// 
    /// Then upload koina_input.tsv to https://koina.wilhelmlab.org/
    ///   - Model: Prosit_2020_intensity_HCD
    ///   - Input columns: modified_sequence, collision_energy, precursor_charge
    ///   - Download the predicted fragment library
    /// </summary>
    public static class KoinaTableGenerator
    {
        public static void Generate(
            string diannTsvPath,
            string outputPath,
            int collisionEnergy = 27)
        {
            Console.WriteLine("=== Koina/Prosit Input Table Generator ===");
            Console.WriteLine($"  Input:  {diannTsvPath}");
            Console.WriteLine($"  Output: {outputPath}");
            Console.WriteLine($"  NCE:    {collisionEnergy}");
            Console.WriteLine();

            if (!File.Exists(diannTsvPath))
            {
                Console.WriteLine($"ERROR: File not found: {diannTsvPath}");
                return;
            }

            var lines = File.ReadAllLines(diannTsvPath);
            if (lines.Length < 2)
            {
                Console.WriteLine("ERROR: File has no data rows.");
                return;
            }

            // Parse header to find column indices
            var header = lines[0].Split('\t');
            var colMap = new Dictionary<string, int>();
            for (int i = 0; i < header.Length; i++)
                colMap[header[i].Trim()] = i;

            int seqCol = FindColumn(colMap, "Modified.Sequence", "Stripped.Sequence", "Sequence");
            int chargeCol = FindColumn(colMap, "Precursor.Charge", "Charge");

            if (seqCol < 0 || chargeCol < 0)
            {
                Console.WriteLine("ERROR: Could not find required columns.");
                Console.WriteLine($"  Available columns: {string.Join(", ", colMap.Keys)}");
                return;
            }

            Console.WriteLine($"  Sequence column: '{header[seqCol]}' (index {seqCol})");
            Console.WriteLine($"  Charge column:   '{header[chargeCol]}' (index {chargeCol})");

            // Collect unique peptide+charge pairs
            var seen = new HashSet<string>();
            var entries = new List<(string koinaSeq, int charge, string originalSeq, string stripped)>();
            int skippedDup = 0, skippedParse = 0, skippedLong = 0;

            for (int row = 1; row < lines.Length; row++)
            {
                var fields = lines[row].Split('\t');
                if (fields.Length <= Math.Max(seqCol, chargeCol))
                {
                    skippedParse++;
                    continue;
                }

                string seq = fields[seqCol].Trim();
                if (!int.TryParse(fields[chargeCol].Trim(), out int charge))
                {
                    skippedParse++;
                    continue;
                }

                string key = seq + "/" + charge;
                if (seen.Contains(key))
                {
                    skippedDup++;
                    continue;
                }
                seen.Add(key);

                // Convert modification format for Koina/Prosit:
                //   C(UniMod:4)  -> C[UNIMOD:4]
                //   M(UniMod:35) -> M[UNIMOD:35]
                string koinaSeq = ConvertToKoinaFormat(seq);
                string stripped = StripModifications(seq);

                // Prosit 2020 accepts sequences up to 30 amino acids
                if (stripped.Length > 30)
                {
                    skippedLong++;
                    continue;
                }

                entries.Add((koinaSeq, charge, seq, stripped));
            }

            // Sort for reproducibility
            entries.Sort((a, b) =>
            {
                int cmp = string.Compare(a.stripped, b.stripped, StringComparison.Ordinal);
                return cmp != 0 ? cmp : a.charge.CompareTo(b.charge);
            });

            // Write output
            using (var writer = new StreamWriter(outputPath))
            {
                writer.WriteLine("modified_sequence\tcollision_energy\tprecursor_charge");

                foreach (var (koinaSeq, charge, originalSeq, stripped) in entries)
                {
                    writer.WriteLine($"{koinaSeq}\t{collisionEnergy}\t{charge}");
                }
            }

            // Summary
            Console.WriteLine();
            Console.WriteLine($"  Results:");
            Console.WriteLine($"    Unique precursors written: {entries.Count:N0}");
            Console.WriteLine($"    Skipped (duplicate):       {skippedDup:N0}");
            Console.WriteLine($"    Skipped (parse error):     {skippedParse:N0}");
            Console.WriteLine($"    Skipped (>30 AA):          {skippedLong:N0}");

            // Charge distribution
            var chargeCounts = new Dictionary<int, int>();
            int hasCAM = 0, hasOx = 0, unmod = 0;
            foreach (var (koinaSeq, charge, _, _) in entries)
            {
                if (!chargeCounts.ContainsKey(charge)) chargeCounts[charge] = 0;
                chargeCounts[charge]++;

                if (koinaSeq.Contains("[UNIMOD:4]")) hasCAM++;
                if (koinaSeq.Contains("[UNIMOD:35]")) hasOx++;
                if (!koinaSeq.Contains("[UNIMOD")) unmod++;
            }

            Console.WriteLine();
            foreach (var kvp in chargeCounts)
                Console.WriteLine($"    Charge {kvp.Key}: {kvp.Value:N0}");

            Console.WriteLine();
            Console.WriteLine($"    Unmodified:          {unmod:N0}");
            Console.WriteLine($"    Carbamidomethyl (C): {hasCAM:N0}");
            Console.WriteLine($"    Oxidation (M):       {hasOx:N0}");
            Console.WriteLine();
            Console.WriteLine($"  Output written to: {outputPath}");
            Console.WriteLine();
            Console.WriteLine("  Next steps:");
            Console.WriteLine("    1. Go to https://koina.wilhelmlab.org/");
            Console.WriteLine("    2. Select model: Prosit_2020_intensity_HCD");
            Console.WriteLine("    3. Upload this TSV file");
            Console.WriteLine("    4. Download the predicted .msp library");
            Console.WriteLine("    5. Use the .msp file as input for the Phase 10 benchmark");
        }

        /// <summary>
        /// Converts DIA-NN modification format to Koina/Prosit bracket notation.
        /// C(UniMod:4) -> C[UNIMOD:4], M(UniMod:35) -> M[UNIMOD:35], etc.
        /// </summary>
        private static string ConvertToKoinaFormat(string seq)
        {
            return Regex.Replace(seq, @"\(UniMod:(\d+)\)", "[UNIMOD:$1]");
        }

        /// <summary>
        /// Removes all UniMod modifications to get the bare amino acid sequence.
        /// </summary>
        private static string StripModifications(string seq)
        {
            return Regex.Replace(seq, @"\(UniMod:\d+\)", "");
        }

        /// <summary>
        /// Finds the first matching column name from a list of candidates.
        /// Returns -1 if none found.
        /// </summary>
        private static int FindColumn(Dictionary<string, int> colMap, params string[] candidates)
        {
            foreach (var name in candidates)
            {
                if (colMap.TryGetValue(name, out int idx))
                    return idx;
            }
            return -1;
        }
    }
}