// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
// Location: MassSpectrometry/Dia/DiaSequenceParser.cs

using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Parses modification-annotated peptide sequences stored in library spectra.
    ///
    /// Library sequences use bracket notation for modifications:
    ///   "PEPTM[Oxidation]IDEK"         → base: "PEPTMIDEK", mods: [(5, "Oxidation")]
    ///   "C[Carbamidomethyl]ARRIER"     → base: "CARRIER",   mods: [(1, "Carbamidomethyl")]
    ///   "M[Oxidation]PEPTM[Oxidation]" → base: "MPEPTM",    mods: [(1,"Oxidation"),(6,"Oxidation")]
    ///
    /// Also handles N-terminal modifications:
    ///   "[Acetyl]PEPTIDE"              → base: "PEPTIDE",   mods: [(0, "Acetyl")]
    ///
    /// Thread-safe: all methods are static with no shared state.
    /// Lives in MassSpectrometry.Dia (mzLib) to avoid circular dependency with MetaMorpheus.
    /// </summary>
    public static class DiaSequenceParser
    {
        // Matches any [...] block: modification name inside brackets
        private static readonly Regex ModPattern = new Regex(
            @"\[([^\[\]]+)\]",
            RegexOptions.Compiled);

        /// <summary>
        /// Strips all modification annotations from a sequence string.
        /// Returns the bare amino acid sequence (one-letter codes only).
        /// </summary>
        /// <param name="annotatedSequence">e.g. "PEPTM[Oxidation]IDEK"</param>
        /// <returns>e.g. "PEPTMIDEK"</returns>
        public static string GetBaseSequence(string annotatedSequence)
        {
            if (string.IsNullOrEmpty(annotatedSequence))
                return annotatedSequence ?? string.Empty;

            // Fast path: no brackets → already a base sequence
            int firstBracket = annotatedSequence.IndexOf('[');
            if (firstBracket < 0)
                return annotatedSequence;

            // Build base sequence by copying characters that are not inside brackets
            var sb = new StringBuilder(annotatedSequence.Length);
            int depth = 0;
            foreach (char c in annotatedSequence)
            {
                if (c == '[') { depth++; continue; }
                if (c == ']') { depth--; continue; }
                if (depth == 0)
                    sb.Append(c);
            }
            return sb.ToString();
        }

        /// <summary>
        /// Parses the full annotated sequence into base sequence + located modifications.
        /// </summary>
        /// <param name="annotatedSequence">e.g. "PEPTM[Oxidation]IDEK"</param>
        /// <param name="baseSequence">Output: the bare amino acid sequence</param>
        /// <param name="modifications">
        /// Output: list of (1-based residue index, modification name).
        /// Index 0 indicates an N-terminal modification (before residue 1).
        /// </param>
        public static void Parse(
            string annotatedSequence,
            out string baseSequence,
            out List<(int OneBasedResidueIndex, string ModificationName)> modifications)
        {
            modifications = new List<(int, string)>();

            if (string.IsNullOrEmpty(annotatedSequence))
            {
                baseSequence = annotatedSequence ?? string.Empty;
                return;
            }

            // Fast path
            if (annotatedSequence.IndexOf('[') < 0)
            {
                baseSequence = annotatedSequence;
                return;
            }

            var baseSb = new StringBuilder(annotatedSequence.Length);
            var modNameSb = new StringBuilder();
            int residueIndex = 0; // tracks position in base sequence (1-based after increment)
            bool inMod = false;

            for (int i = 0; i < annotatedSequence.Length; i++)
            {
                char c = annotatedSequence[i];

                if (c == '[')
                {
                    inMod = true;
                    modNameSb.Clear();
                    continue;
                }

                if (c == ']')
                {
                    inMod = false;
                    // residueIndex is the 1-based position of the preceding amino acid,
                    // or 0 if the mod comes before any residue (N-terminal)
                    modifications.Add((residueIndex, modNameSb.ToString()));
                    continue;
                }

                if (inMod)
                {
                    modNameSb.Append(c);
                }
                else
                {
                    baseSb.Append(c);
                    residueIndex++;
                }
            }

            baseSequence = baseSb.ToString();
        }

        /// <summary>
        /// Returns true if the sequence string contains any modification annotations.
        /// </summary>
        public static bool HasModifications(string annotatedSequence)
        {
            return !string.IsNullOrEmpty(annotatedSequence)
                   && annotatedSequence.IndexOf('[') >= 0;
        }

        /// <summary>
        /// Formats a human-readable modification summary string.
        /// e.g. "Oxidation@5; Carbamidomethyl@1"
        /// Returns empty string if no modifications.
        /// </summary>
        public static string FormatModificationSummary(string annotatedSequence)
        {
            if (!HasModifications(annotatedSequence))
                return string.Empty;

            Parse(annotatedSequence, out _, out var mods);
            if (mods.Count == 0)
                return string.Empty;

            var sb = new StringBuilder();
            for (int i = 0; i < mods.Count; i++)
            {
                if (i > 0) sb.Append("; ");
                sb.Append(mods[i].ModificationName);
                sb.Append('@');
                sb.Append(mods[i].OneBasedResidueIndex);
            }
            return sb.ToString();
        }
    }
}
