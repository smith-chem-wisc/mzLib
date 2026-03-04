// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Text;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Parses modification-annotated peptide sequences used in spectral libraries.
    ///
    /// Handles bracket notation: <c>"PEPTM[Oxidation]IDEK"</c>
    ///   → base sequence <c>"PEPTMIDEK"</c>
    ///   → modifications <c>[(5, "Oxidation")]</c>  (1-based residue position)
    ///
    /// This class lives in mzLib because library sequences are a mzLib concern.
    /// MetaMorpheus calls it via <see cref="DiaPsmAdapter"/> to populate
    /// <c>BaseSequence</c> and <c>ModificationSummary</c> columns in .psmtsv output.
    /// </summary>
    public static class DiaSequenceParser
    {
        /// <summary>
        /// Returns the bare amino-acid sequence with all bracket modifications stripped.
        /// </summary>
        /// <example>
        /// <c>GetBaseSequence("PEPTM[Oxidation]IDEK")</c> → <c>"PEPTMIDEK"</c>
        /// </example>
        public static string GetBaseSequence(string annotatedSequence)
        {
            if (string.IsNullOrEmpty(annotatedSequence))
                return annotatedSequence ?? string.Empty;

            var sb = new StringBuilder(annotatedSequence.Length);
            bool inBracket = false;

            foreach (char c in annotatedSequence)
            {
                if (c == '[') { inBracket = true; continue; }
                if (c == ']') { inBracket = false; continue; }
                if (!inBracket) sb.Append(c);
            }

            return sb.ToString();
        }

        /// <summary>
        /// Parses an annotated sequence into a base sequence and a list of
        /// (1-based residue position, modification name) tuples.
        /// </summary>
        /// <param name="annotatedSequence">e.g. <c>"PEPTM[Oxidation]IDEK"</c></param>
        /// <param name="baseSequence">Output bare sequence, e.g. <c>"PEPTMIDEK"</c></param>
        /// <param name="modifications">
        /// Output list of <c>(int Position, string Name)</c> value tuples.
        /// Position is 1-based and refers to the residue immediately before the bracket.
        /// </param>
        public static void Parse(
            string annotatedSequence,
            out string baseSequence,
            out List<(int Position, string Name)> modifications)
        {
            modifications = new List<(int, string)>();

            if (string.IsNullOrEmpty(annotatedSequence))
            {
                baseSequence = annotatedSequence ?? string.Empty;
                return;
            }

            var sb = new StringBuilder(annotatedSequence.Length);
            var modName = new StringBuilder();
            bool inBracket = false;
            int residuePos = 0; // 1-based position in base sequence

            foreach (char c in annotatedSequence)
            {
                if (c == '[')
                {
                    inBracket = true;
                    modName.Clear();
                    continue;
                }

                if (c == ']')
                {
                    inBracket = false;
                    modifications.Add((residuePos, modName.ToString()));
                    continue;
                }

                if (inBracket)
                {
                    modName.Append(c);
                }
                else
                {
                    sb.Append(c);
                    residuePos++;
                }
            }

            baseSequence = sb.ToString();
        }

        /// <summary>
        /// Returns a compact modification summary string, e.g. <c>"Oxidation@5"</c>
        /// or <c>"Phospho@3; Oxidation@7"</c> for multiple modifications.
        /// Returns <see cref="string.Empty"/> when no modifications are present.
        /// </summary>
        public static string FormatModificationSummary(string annotatedSequence)
        {
            Parse(annotatedSequence, out _, out var mods);

            if (mods.Count == 0)
                return string.Empty;

            var parts = new string[mods.Count];
            for (int i = 0; i < mods.Count; i++)
                parts[i] = $"{mods[i].Name}@{mods[i].Position}";

            return string.Join("; ", parts);
        }
    }
}