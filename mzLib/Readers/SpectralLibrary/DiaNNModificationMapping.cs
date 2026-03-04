using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Provides bidirectional mapping between DIA-NN's UniMod modification notation
    /// and mzLib's descriptive modification notation.
    /// 
    /// DIA-NN format:  _PEPTM[UniMod:35]IDE_   (flanking underscores, UniMod IDs)
    /// mzLib format:   PEPTM[Common Variable:Oxidation on M]IDE   (no underscores, descriptive names)
    /// 
    /// Also handles neutral loss type mapping between DIA-NN string names and numeric masses.
    /// </summary>
    public static class DiaNNModificationMapping
    {
        #region UniMod ID ↔ mzLib Name Mapping

        /// <summary>
        /// Maps UniMod accession IDs to mzLib modification names.
        /// This is the core mapping table for converting DIA-NN sequences to mzLib format.
        /// 
        /// When adding new modifications:
        /// 1. Find the UniMod accession at https://www.unimod.org/
        /// 2. Find the corresponding mzLib name in the MetaMorpheus modification database
        /// 3. Add both the forward and reverse mapping
        /// </summary>
        public static readonly Dictionary<int, string> UniModIdToMzLibName = new()
        {
            { 4,   "[Common Fixed:Carbamidomethyl on C]" },
            { 35,  "[Common Variable:Oxidation on M]" },
            { 1,   "[Common Biological:Acetylation on X]" },
            { 21,  "[Common Biological:Phosphorylation on S]" },   // Also T, Y — residue context determines
            { 7,   "[Common Biological:Deamidation on N]" },       // Also Q
            { 27,  "[Common Biological:Glu to PyroGlu on E]" },    // Glu->pyro-Glu (N-term E)
            { 28,  "[Common Biological:Gln to PyroGlu on Q]" },    // Gln->pyro-Glu (N-term Q)
            { 526, "[Common Biological:Sodium on D]" },            // Na adduct — also E
            { 34,  "[Common Variable:Methylation on K]" },         // Also R
            { 36,  "[Common Variable:Dimethylation on K]" },       // Also R
            { 37,  "[Common Variable:Trimethylation on K]" },
            { 121, "[Common Biological:Ubiquitination (GlyGly) on K]" },
        };

        /// <summary>
        /// Reverse mapping: mzLib modification name → UniMod accession ID.
        /// Built automatically from UniModIdToMzLibName at static initialization.
        /// </summary>
        public static readonly Dictionary<string, int> MzLibNameToUniModId;

        /// <summary>
        /// Maps UniMod accession IDs to their monoisotopic mass shifts (Da).
        /// Used as a fallback when a UniMod ID is not in the name mapping table.
        /// Masses sourced from https://www.unimod.org/
        /// </summary>
        public static readonly Dictionary<int, double> UniModIdToMonoisotopicMass = new()
        {
            { 4,   57.021464 },   // Carbamidomethyl
            { 35,  15.994915 },   // Oxidation
            { 1,   42.010565 },   // Acetyl
            { 21,  79.966331 },   // Phospho
            { 7,   0.984016 },    // Deamidated
            { 27,  -18.010565 },  // Glu->pyro-Glu
            { 28,  -17.026549 },  // Gln->pyro-Glu
            { 526, 21.981943 },   // Sodium (Na - H)
            { 34,  14.015650 },   // Methyl
            { 36,  28.031300 },   // Dimethyl
            { 37,  42.046950 },   // Trimethyl
            { 121, 114.042927 },  // GlyGly (ubiquitin)
        };

        #endregion

        #region Neutral Loss Mapping

        /// <summary>
        /// Maps DIA-NN neutral loss type strings to their monoisotopic mass values.
        /// DIA-NN uses string identifiers in TSV; mzLib uses numeric masses in Product.NeutralLoss.
        /// </summary>
        public static readonly Dictionary<string, double> NeutralLossNameToMass = new(StringComparer.OrdinalIgnoreCase)
        {
            { "noloss", 0.0 },
            { "",       0.0 },
            { "H2O",    18.010565 },
            { "NH3",    17.026549 },
            { "H3PO4",  97.976896 },
            { "HPO3",   79.966331 },
            { "CO",     27.994915 },
        };

        /// <summary>
        /// Reverse mapping: neutral loss mass → DIA-NN string name.
        /// Uses a tolerance of 0.001 Da for matching.
        /// </summary>
        public static string MassToNeutralLossName(double mass)
        {
            if (Math.Abs(mass) < 0.001) return "noloss";
            foreach (var kvp in NeutralLossNameToMass)
            {
                if (Math.Abs(kvp.Value - mass) < 0.001)
                    return kvp.Key;
            }
            // Unknown loss — return mass as string for round-trip safety
            return mass.ToString("F4");
        }

        #endregion

        #region Regex Patterns

        /// <summary>
        /// Matches UniMod modification tags in DIA-NN format: [UniMod:35]
        /// Group 1 captures the UniMod accession ID as a string.
        /// </summary>
        private static readonly Regex UniModTagRegex = new(@"\[UniMod:(\d+)\]", RegexOptions.Compiled);

        /// <summary>
        /// Matches mzLib modification tags: [Common Variable:Oxidation on M], [Common Fixed:Carbamidomethyl on C], etc.
        /// Captures the entire bracketed tag including brackets.
        /// </summary>
        private static readonly Regex MzLibModTagRegex = new(@"\[[^\]]+\]", RegexOptions.Compiled);

        #endregion

        #region Static Constructor

        static DiaNNModificationMapping()
        {
            // Build reverse mapping from the forward mapping
            MzLibNameToUniModId = new Dictionary<string, int>();
            foreach (var kvp in UniModIdToMzLibName)
            {
                // Only add if not already present (first mapping wins for duplicates)
                MzLibNameToUniModId.TryAdd(kvp.Value, kvp.Key);
            }
        }

        #endregion

        #region Sequence Conversion Methods

        /// <summary>
        /// Converts a DIA-NN modified peptide sequence to mzLib format.
        /// 
        /// Input:  _PEPTM[UniMod:35]IDE_       (DIA-NN format)
        /// Output: PEPTM[Common Variable:Oxidation on M]IDE   (mzLib format)
        /// 
        /// Steps:
        /// 1. Strip flanking underscores
        /// 2. Replace each [UniMod:N] tag with the corresponding mzLib name
        /// 3. For unknown UniMod IDs, fall back to mass notation: [+15.9949]
        /// </summary>
        /// <param name="diannSequence">Modified peptide in DIA-NN format (with or without flanking underscores)</param>
        /// <returns>Modified peptide in mzLib format</returns>
        public static string DiaNNToMzLib(string diannSequence)
        {
            if (string.IsNullOrEmpty(diannSequence))
                return diannSequence;

            // Step 1: Strip flanking underscores
            string sequence = StripFlankingUnderscores(diannSequence);

            // Step 2: Replace UniMod tags with mzLib names
            sequence = UniModTagRegex.Replace(sequence, match =>
            {
                int unimodId = int.Parse(match.Groups[1].Value);

                if (UniModIdToMzLibName.TryGetValue(unimodId, out var mzLibName))
                {
                    return mzLibName;
                }

                // Fallback: use mass notation if we know the mass
                if (UniModIdToMonoisotopicMass.TryGetValue(unimodId, out var mass))
                {
                    string sign = mass >= 0 ? "+" : "";
                    return $"[{sign}{mass:F6}]";
                }

                // Last resort: preserve the original tag with a warning-friendly format
                return $"[UniMod:{unimodId}]";
            });

            return sequence;
        }

        /// <summary>
        /// Converts an mzLib modified peptide sequence to DIA-NN format.
        /// 
        /// Input:  PEPTM[Common Variable:Oxidation on M]IDE   (mzLib format)
        /// Output: _PEPTM[UniMod:35]IDE_                      (DIA-NN format)
        /// 
        /// Steps:
        /// 1. Replace each mzLib modification tag with [UniMod:N]
        /// 2. Add flanking underscores
        /// </summary>
        /// <param name="mzLibSequence">Modified peptide in mzLib format</param>
        /// <returns>Modified peptide in DIA-NN format with flanking underscores</returns>
        public static string MzLibToDiaNN(string mzLibSequence)
        {
            if (string.IsNullOrEmpty(mzLibSequence))
                return mzLibSequence;

            string sequence = MzLibModTagRegex.Replace(mzLibSequence, match =>
            {
                string modTag = match.Value;

                if (MzLibNameToUniModId.TryGetValue(modTag, out var unimodId))
                {
                    return $"[UniMod:{unimodId}]";
                }

                // Unknown mod — preserve as-is for debugging visibility
                return modTag;
            });

            // Add flanking underscores
            return $"_{sequence}_";
        }

        /// <summary>
        /// Strips the flanking underscores that DIA-NN uses to delimit peptide sequences.
        /// Handles sequences with or without underscores gracefully.
        /// </summary>
        public static string StripFlankingUnderscores(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return sequence;

            if (sequence.StartsWith('_'))
                sequence = sequence.Substring(1);
            if (sequence.EndsWith('_'))
                sequence = sequence.Substring(0, sequence.Length - 1);

            return sequence;
        }

        /// <summary>
        /// Extracts the stripped (unmodified) sequence from a modified peptide string.
        /// Removes all bracketed modifications and flanking underscores.
        /// 
        /// Input:  _PEPTM[UniMod:35]IDE_  → PEPTMIDE
        /// Input:  PEPTM[Common Variable:Oxidation on M]IDE → PEPTMIDE
        /// </summary>
        public static string GetStrippedSequence(string modifiedSequence)
        {
            if (string.IsNullOrEmpty(modifiedSequence))
                return modifiedSequence;

            string stripped = StripFlankingUnderscores(modifiedSequence);
            stripped = MzLibModTagRegex.Replace(stripped, string.Empty);
            return stripped;
        }

        #endregion
    }
}
