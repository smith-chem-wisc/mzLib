// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// DecoyKoinaTableGenerator.cs
// Location: Development/Dia/ (or run standalone)
//
// Reads a Koina input TSV (modified_sequence, collision_energy, precursor_charge)
// and produces a decoy version where each sequence is pseudo-reversed:
//   - Last amino acid stays in place (preserves tryptic C-terminus)
//   - Remaining residues are reversed
//   - Modifications stay attached to their original amino acid
//
// Output: same format TSV, ready to submit to Koina for RT + fragmentation prediction.
//
// Usage:
//   DecoyKoinaTableGenerator.Generate(
//       @"F:\DiaBenchmark\PXD005573\DiannOut\koina_input.tsv",
//       @"F:\DiaBenchmark\PXD005573\DiannOut\koina_decoy_input.tsv");

using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace Development.Dia
{
    public static class DecoyKoinaTableGenerator
    {
        /// <summary>
        /// Reads a Koina TSV and writes a decoy version with reversed sequences.
        /// </summary>
        public static void Generate(string inputTsvPath, string outputTsvPath)
        {
            if (!File.Exists(inputTsvPath))
                throw new FileNotFoundException($"Input file not found: {inputTsvPath}");

            var lines = File.ReadAllLines(inputTsvPath);
            if (lines.Length < 2)
                throw new InvalidOperationException("Input file has no data rows.");

            // Parse header
            string header = lines[0];
            string[] headerCols = header.Split('\t');

            int seqCol = -1, energyCol = -1, chargeCol = -1;
            for (int i = 0; i < headerCols.Length; i++)
            {
                string col = headerCols[i].Trim().ToLowerInvariant();
                if (col == "modified_sequence") seqCol = i;
                else if (col == "collision_energy") energyCol = i;
                else if (col == "precursor_charge") chargeCol = i;
            }

            if (seqCol < 0 || energyCol < 0 || chargeCol < 0)
                throw new InvalidOperationException(
                    $"Expected columns: modified_sequence, collision_energy, precursor_charge. Found: {header}");

            int total = 0;
            int duplicatesAvoided = 0;
            var seenDecoys = new HashSet<string>();

            using var writer = new StreamWriter(outputTsvPath);
            writer.WriteLine(header); // same header

            for (int lineIdx = 1; lineIdx < lines.Length; lineIdx++)
            {
                string line = lines[lineIdx];
                if (string.IsNullOrWhiteSpace(line)) continue;

                string[] cols = line.Split('\t');
                if (cols.Length <= Math.Max(seqCol, Math.Max(energyCol, chargeCol)))
                    continue;

                string targetSeq = cols[seqCol].Trim();
                string energy = cols[energyCol].Trim();
                string charge = cols[chargeCol].Trim();

                string decoySeq = ReverseSequenceKeepLastAA(targetSeq);

                // Ensure decoy != target (rare but possible for palindromic sequences)
                if (decoySeq == targetSeq)
                {
                    // Swap first two residues as fallback
                    decoySeq = SwapFirstTwoResidues(decoySeq);
                    if (decoySeq == targetSeq)
                    {
                        // Give up — skip this entry (should be extremely rare)
                        duplicatesAvoided++;
                        continue;
                    }
                }

                // Deduplicate (same sequence + charge = redundant decoy)
                string key = $"{decoySeq}\t{charge}";
                if (!seenDecoys.Add(key))
                {
                    duplicatesAvoided++;
                    continue;
                }

                writer.WriteLine($"{decoySeq}\t{energy}\t{charge}");
                total++;
            }

            Console.WriteLine($"Decoy table generated: {outputTsvPath}");
            Console.WriteLine($"  Target entries read: {lines.Length - 1}");
            Console.WriteLine($"  Decoy entries written: {total}");
            if (duplicatesAvoided > 0)
                Console.WriteLine($"  Skipped (palindrome/duplicate): {duplicatesAvoided}");
        }

        /// <summary>
        /// Reverses a modified peptide sequence, keeping the last amino acid in place
        /// and keeping modifications attached to their original amino acid.
        /// 
        /// Input format examples:
        ///   PEPTIDER                          → EDITPEPR
        ///   AAAGEFADDPC[UNIMOD:4]SSVK         → VSSC[UNIMOD:4]PDDAFEGAAAK
        ///   AAAIGIDLGTTYSC[UNIMOD:4]VGVFQHGK  → GHQFVGVC[UNIMOD:4]SYTTGLDLGIAAAK
        ///   AAAPAPVSEAVC[UNIMOD:4]R            → C[UNIMOD:4]VAESVSAPAAPAR  (note: mod stays with C)
        /// 
        /// Algorithm:
        ///   1. Parse sequence into tokens: each token is an amino acid + optional modification
        ///   2. Separate last token (C-terminal AA)
        ///   3. Reverse the remaining tokens
        ///   4. Re-append the last token
        /// </summary>
        public static string ReverseSequenceKeepLastAA(string modifiedSequence)
        {
            if (string.IsNullOrEmpty(modifiedSequence))
                return modifiedSequence;

            // Parse into residue tokens: "AA" or "AA[UNIMOD:N]"
            var tokens = ParseResidueTokens(modifiedSequence);

            if (tokens.Count <= 2)
                return modifiedSequence; // Nothing meaningful to reverse

            // Last token stays in place
            string lastToken = tokens[tokens.Count - 1];

            // Reverse everything except the last
            var reversed = new StringBuilder();
            for (int i = tokens.Count - 2; i >= 0; i--)
                reversed.Append(tokens[i]);

            reversed.Append(lastToken);
            return reversed.ToString();
        }

        /// <summary>
        /// Parses a modified sequence into tokens where each token is one amino acid
        /// plus any immediately following modification bracket.
        /// 
        /// Example: "AAPC[UNIMOD:4]SVK" → ["A", "A", "P", "C[UNIMOD:4]", "S", "V", "K"]
        /// </summary>
        private static List<string> ParseResidueTokens(string sequence)
        {
            var tokens = new List<string>();
            int i = 0;

            while (i < sequence.Length)
            {
                if (!char.IsLetter(sequence[i]))
                {
                    // Shouldn't happen in well-formed input, but be safe
                    i++;
                    continue;
                }

                // Start of a residue: one uppercase letter
                int start = i;
                i++; // consume the amino acid letter

                // Check if followed by a modification in brackets
                if (i < sequence.Length && sequence[i] == '[')
                {
                    // Find matching closing bracket
                    int bracketDepth = 1;
                    i++; // skip opening '['
                    while (i < sequence.Length && bracketDepth > 0)
                    {
                        if (sequence[i] == '[') bracketDepth++;
                        else if (sequence[i] == ']') bracketDepth--;
                        i++;
                    }
                }

                tokens.Add(sequence.Substring(start, i - start));
            }

            return tokens;
        }

        /// <summary>
        /// Fallback for palindromic sequences: swap the first two residue tokens.
        /// </summary>
        private static string SwapFirstTwoResidues(string modifiedSequence)
        {
            var tokens = ParseResidueTokens(modifiedSequence);
            if (tokens.Count < 3) return modifiedSequence;

            // Swap tokens[0] and tokens[1]
            (tokens[0], tokens[1]) = (tokens[1], tokens[0]);

            var sb = new StringBuilder();
            foreach (var t in tokens) sb.Append(t);
            return sb.ToString();
        }
    }
}