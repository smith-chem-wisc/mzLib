using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Test.Omics.SequenceConversion
{
    /// <summary>
    /// Ground truth test cases with known-good sequence representations in all formats.
    /// Each test case represents the SAME logical sequence in different format representations.
    /// 
    /// FORMAT RULES:
    /// 
    /// mzLib: Modification annotations go AFTER the residue they modify
    /// - Format: `RESIDUE[ModificationName on Residue]`
    /// - Example: `M[Common Variable:Oxidation on M]` or `S[Common Biological:Phosphorylation on S]`
    /// - N-term: `[ModificationName on X]SEQUENCE`
    /// - C-term: `SEQUENCE-[ModificationName on X]`
    /// 
    /// MassShift: Mass shifts in brackets go AFTER the residue (same positioning as mzLib)
    /// - Format: `RESIDUE[+massShift]` or `RESIDUE[-massShift]`
    /// - Example: `M[+15.9949]` or `S[+79.9663]`
    /// - N-term: `[+mass]SEQUENCE`
    /// - C-term: `SEQUENCE[-mass]` or `SEQUENCE[+mass]`
    /// 
    /// Chronologer: Single-character encoding (lowercase for modified residues)
    /// - Format: `{NTerm}{SEQUENCE}{CTerm}`
    /// - NTerm codes: `-` (free), `^` (acetyl), `)` (pyroGlu), `(` (cyclized CAM-Cys), `&` (GlyGly), `*` (heavy GlyGly)
    /// - CTerm: always `_`
    /// - Modified residues: `m` (oxidized M), `c` (carbamidomethyl C), `s/t/y` (phospho S/T/Y), `a` (acetyl K), etc.
    /// - Example: `-PEPmTIDE_` (oxidized methionine), `^PEPTIDE_` (N-term acetyl)
    /// </summary>
    public static class GroundTruthTestData
    {
        /// <summary>
        /// Represents a single test case with the same sequence in multiple formats
        /// </summary>
        public class SequenceConversionTestCase
        {
            public string Description { get; set; }
            public string MzLibFormat { get; set; }
            public string MassShiftFormat { get; set; }
            public string ChronologerFormat { get; set; }
            public string UnimodUpperCaseFormat { get; set; }
            public string UnimodCamelCaseFormat => UnimodUpperCaseFormat.Replace("UNIMOD", "Unimod");
            public string UnimodLowerCaseFormat => UnimodUpperCaseFormat.Replace("UNIMOD", "unimod");
            public string UnimodNoLabelFormat => UnimodUpperCaseFormat.Replace("UNIMOD:", "");
            
            /// <summary>
            /// FDRBench format uses parentheses and CamelCase: "PEPC(UniMod:4)IDE"
            /// </summary>
            public string FdrBenchFormat => UnimodCamelCaseFormat.Replace("[", "(").Replace("]", ")");
            
            public string? UniProtFormat { get; set; }
            public string? UniProtFormatWithModType => UniProtFormat is null ? null : $"UniProt:{UniProtFormat}";
            public string ExpectedBaseSequence { get; set; }
            public bool HasNTermMod { get; set; }
            public bool HasCTermMod { get; set; }
            public int ExpectedResidueMods { get; set; }
            
            public SequenceConversionTestCase(string description, string mzLib, string massShift, string chronologer,
                string baseSeq, bool hasNTerm = false, bool hasCTerm = false, int residueMods = 0,
                string? unimodUpperCase = null, string? uniProtFormat = null)
            {
                Description = description;
                MzLibFormat = mzLib;
                MassShiftFormat = massShift;
                ChronologerFormat = chronologer;
                UnimodUpperCaseFormat = unimodUpperCase;
                ExpectedBaseSequence = baseSeq;
                HasNTermMod = hasNTerm;
                HasCTermMod = hasCTerm;
                ExpectedResidueMods = residueMods;
                UniProtFormat = uniProtFormat;
            }

            public override string ToString() => Description ?? MzLibFormat;
        }

        /// <summary>
        /// Core test cases covering essential functionality with Chronologer-compatible modifications only
        /// </summary>
        public static readonly SequenceConversionTestCase[] CoreTestCases = new[]
        {
            // 1. Unmodified sequence (simplest case)
            new SequenceConversionTestCase(
                description: "Unmodified peptide",
                mzLib: "PEPTIDE",
                massShift: "PEPTIDE",
                chronologer: "-PEPTIDE_",
                baseSeq: "PEPTIDE",
                unimodUpperCase: "PEPTIDE",
                uniProtFormat: "PEPTIDE"
            ),

            // 2. Common Variable:Oxidation on M - UNIMOD:35, mass +15.9949
            new SequenceConversionTestCase(
                description: "Single Common Variable:Oxidation on Methionine",
                mzLib: "PEPTM[Common Variable:Oxidation on M]IDE",
                massShift: "PEPTM[+15.9949]IDE",
                chronologer: "-PEPTmIDE_",
                baseSeq: "PEPTMIDE",
                residueMods: 1,
                unimodUpperCase: "PEPTM[UNIMOD:35]IDE",
                uniProtFormat: "PEPTM[Methionine sulfoxide]IDE"
            ),

            // 3. Common Fixed:Carbamidomethyl on C - UNIMOD:4, mass +57.0215
            new SequenceConversionTestCase(
                description: "Common Fixed:Carbamidomethyl on Cysteine",
                mzLib: "PEPC[Common Fixed:Carbamidomethyl on C]IDE",
                massShift: "PEPC[+57.0215]IDE",
                chronologer: "-PEPcIDE_",
                baseSeq: "PEPCIDE",
                residueMods: 1,
                unimodUpperCase: "PEPC[UNIMOD:4]IDE",
                uniProtFormat: "PEPCIDE"
            ),

            // 4. Phosphorylation on S - UNIMOD:21, mass +79.9663
            new SequenceConversionTestCase(
                description: "Phosphorylation on serine",
                mzLib: "PEPS[Common Biological:Phosphorylation on S]TIDE",
                massShift: "PEPS[+79.9663]TIDE",
                chronologer: "-PEPsTIDE_",
                baseSeq: "PEPSTIDE",
                residueMods: 1,
                unimodUpperCase: "PEPS[UNIMOD:21]TIDE",
                uniProtFormat: "PEPS[Phosphoserine]TIDE"
            ),

            // 5. Phosphorylation on T - UNIMOD:21, mass +79.9663
            new SequenceConversionTestCase(
                description: "Phosphorylation on threonine",
                mzLib: "PEPT[Common Biological:Phosphorylation on T]IDE",
                massShift: "PEPT[+79.9663]IDE",
                chronologer: "-PEPtIDE_",
                baseSeq: "PEPTIDE",
                residueMods: 1,
                unimodUpperCase: "PEPT[UNIMOD:21]IDE",
                uniProtFormat: "PEPT[Phosphothreonine]IDE"
            ),

            // 6. N-terminal acetylation - UNIMOD:1, mass +42.0106
            new SequenceConversionTestCase(
                description: "N-terminal acetylation",
                mzLib: "[Common Biological:Acetylation on X]PEPTIDE",
                massShift: "[+42.0106]PEPTIDE",
                chronologer: "^PEPTIDE_",
                baseSeq: "PEPTIDE",
                hasNTerm: true,
                unimodUpperCase: "[UNIMOD:1]PEPTIDE",
                uniProtFormat: "[N-acetylproline]PEPTIDE"
            ),

            // 7. Multiple modifications: oxidation and phosphorylation
            new SequenceConversionTestCase(
                description: "Common Variable:Oxidation on M and phosphorylation on S",
                mzLib: "PEPS[Common Biological:Phosphorylation on S]TM[Common Variable:Oxidation on M]IDE",
                massShift: "PEPS[+79.9663]TM[+15.9949]IDE",
                chronologer: "-PEPsTmIDE_",
                baseSeq: "PEPSTMIDE",
                residueMods: 2,
                unimodUpperCase: "PEPS[UNIMOD:21]TM[UNIMOD:35]IDE",
                uniProtFormat: "PEPS[Phosphoserine]TM[Methionine sulfoxide]IDE"
            ),

            // 8. N-term acetylation + internal oxidation
            new SequenceConversionTestCase(
                description: "N-terminal acetylation with internal oxidation",
                mzLib: "[Common Biological:Acetylation on X]PEPTM[Common Variable:Oxidation on M]IDE",
                massShift: "[+42.0106]PEPTM[+15.9949]IDE",
                chronologer: "^PEPTmIDE_",
                baseSeq: "PEPTMIDE",
                hasNTerm: true,
                residueMods: 1,
                unimodUpperCase: "[UNIMOD:1]PEPTM[UNIMOD:35]IDE",
                uniProtFormat: "[N-acetylproline]PEPTM[Methionine sulfoxide]IDE"
            ),

            // 9. Two oxidations on different methionines
            new SequenceConversionTestCase(
                description: "Two oxidations on different methionines",
                mzLib: "M[Common Variable:Oxidation on M]EPTM[Common Variable:Oxidation on M]IDE",
                massShift: "M[+15.9949]EPTM[+15.9949]IDE",
                chronologer: "-mEPTmIDE_",
                baseSeq: "MEPTMIDE",
                residueMods: 2,
                unimodUpperCase: "M[UNIMOD:35]EPTM[UNIMOD:35]IDE",
                uniProtFormat: "M[Methionine sulfoxide]EPTM[Methionine sulfoxide]IDE"
            ),

            // 10. Water loss (negative mass shift)
            new SequenceConversionTestCase(
                description: "dehydrobutyrine on T (negative mass shift)",
                mzLib: "PEPT[Less Common:Dehydrobutyrine on T]TIDE",
                massShift: "PEPT[-18.0106]TIDE",
                chronologer: "-PEPTTIDE_",
                baseSeq: "PEPTTIDE",
                residueMods: 1,
                unimodUpperCase: "PEPT[UNIMOD:23]TIDE",
                uniProtFormat: "PEPT[2,3-didehydrobutyrine]TIDE"
                ),

            // 11. N-term pyroglutamate formation (loss of water, -18.0106 Da)
            new SequenceConversionTestCase(
                description: "N-term pyroglutamate formation (loss of water, -18.0106 Da)",
                mzLib: "[Common Artifact:Water Loss on E]EPTIDE",
                massShift: "[-18.0106]EPTIDE",
                chronologer: ")EPTIDE_",
                baseSeq: "EPTIDE",
                hasNTerm: true,
                residueMods: 0,
                unimodUpperCase: "[UNIMOD:27]EPTIDE",
                uniProtFormat: "[Pyrrolidone carboxylic acid (Glu)]EPTIDE"
                ),

            // 12. C-term modification
            new SequenceConversionTestCase(
                description: "C-terminal amidation",
                mzLib: "PEPTIDE-[Common Biological:Amidation on X]",
                massShift: "PEPTIDE-[-0.984]",
                chronologer: "-PEPTIDE_",
                baseSeq: "PEPTIDE",
                hasCTerm: true,
                residueMods: 0,
                unimodUpperCase: "PEPTIDE-[UNIMOD:2]",
                uniProtFormat: "PEPTIDE-[Glutamic acid 1-amide]"
                )
        };

        /// <summary>
        /// Edge cases and special scenarios
        /// </summary>
        public static readonly SequenceConversionTestCase[] EdgeCases = new[]
        {
            // Single amino acid
            new SequenceConversionTestCase(
                description: "Single unmodified residue",
                mzLib: "K",
                massShift: "K",
                chronologer: "-K_",
                baseSeq: "K",
                unimodUpperCase: "K"
            ),

            // Single modified residue
            new SequenceConversionTestCase(
                description: "Single oxidized methionine",
                mzLib: "M[Common Variable:Oxidation on M]",
                massShift: "M[+15.9949]",
                chronologer: "-m_",
                baseSeq: "M",
                residueMods: 1,
                unimodUpperCase: "M[UNIMOD:35]",
                uniProtFormat: "M[Methionine sulfoxide]"
            ),

            // Adjacent modifications
            new SequenceConversionTestCase(
                description: "Two adjacent oxidized methionines",
                mzLib: "PEPTM[Common Variable:Oxidation on M]M[Common Variable:Oxidation on M]IDE",
                massShift: "PEPTM[+15.9949]M[+15.9949]IDE",
                chronologer: "-PEPTmmIDE_",
                baseSeq: "PEPTMMIDE",
                residueMods: 2,
                unimodUpperCase: "PEPTM[UNIMOD:35]M[UNIMOD:35]IDE",
                uniProtFormat: "PEPTM[Methionine sulfoxide]M[Methionine sulfoxide]IDE"
            ),
        };
    }
}
