using System.Linq;

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
        public class TestCase
        {
            public string Description { get; set; }
            public string MzLibFormat { get; set; }
            public string MassShiftFormat { get; set; }
            public string ChronologerFormat { get; set; }
            
            // Expected canonical representation
            public string ExpectedBaseSequence { get; set; }
            public bool HasNTermMod { get; set; }
            public bool HasCTermMod { get; set; }
            public int ExpectedResidueMods { get; set; }
            
            public TestCase(string description, string mzLib, string massShift, string chronologer,
                string baseSeq, bool hasNTerm = false, bool hasCTerm = false, int residueMods = 0)
            {
                Description = description;
                MzLibFormat = mzLib;
                MassShiftFormat = massShift;
                ChronologerFormat = chronologer;
                ExpectedBaseSequence = baseSeq;
                HasNTermMod = hasNTerm;
                HasCTermMod = hasCTerm;
                ExpectedResidueMods = residueMods;
            }
        }

        /// <summary>
        /// Core test cases covering essential functionality with Chronologer-compatible modifications only
        /// </summary>
        public static readonly TestCase[] CoreTestCases = new[]
        {
            // 1. Unmodified sequence (simplest case)
            new TestCase(
                description: "Unmodified peptide",
                mzLib: "PEPTIDE",
                massShift: "PEPTIDE",
                chronologer: "-PEPTIDE_",
                baseSeq: "PEPTIDE"
            ),

            // 2. Common Variable:Oxidation on M - UNIMOD:35, mass +15.9949
            new TestCase(
                description: "Single Common Variable:Oxidation on Methionine",
                mzLib: "PEPTM[Common Variable:Oxidation on M]IDE",
                massShift: "PEPTM[+15.9949]IDE",
                chronologer: "-PEPTmIDE_",
                baseSeq: "PEPTMIDE",
                residueMods: 1
            ),

            // 3. Common Fixed:Carbamidomethyl on C - UNIMOD:4, mass +57.0215
            new TestCase(
                description: "Common Fixed:Carbamidomethyl on Cysteine",
                mzLib: "PEPC[Common Fixed:Carbamidomethyl on C]IDE",
                massShift: "PEPC[+57.0215]IDE",
                chronologer: "-PEPcIDE_",
                baseSeq: "PEPCIDE",
                residueMods: 1
            ),

            // 4. Phosphorylation on S - UNIMOD:21, mass +79.9663
            new TestCase(
                description: "Phosphorylation on serine",
                mzLib: "PEPS[Common Biological:Phosphorylation on S]TIDE",
                massShift: "PEPS[+79.9663]TIDE",
                chronologer: "-PEPsTIDE_",
                baseSeq: "PEPSTIDE",
                residueMods: 1
            ),

            // 5. Phosphorylation on T - UNIMOD:21, mass +79.9663
            new TestCase(
                description: "Phosphorylation on threonine",
                mzLib: "PEPT[Common Biological:Phosphorylation on T]IDE",
                massShift: "PEPT[+79.9663]IDE",
                chronologer: "-PEPtIDE_",
                baseSeq: "PEPTIDE",
                residueMods: 1
            ),

            // 6. N-terminal acetylation - UNIMOD:1, mass +42.0106
            new TestCase(
                description: "N-terminal acetylation",
                mzLib: "[Common Biological:Acetylation on X]PEPTIDE",
                massShift: "[+42.0106]PEPTIDE",
                chronologer: "^PEPTIDE_",
                baseSeq: "PEPTIDE",
                hasNTerm: true
            ),

            // 7. Multiple modifications: oxidation and phosphorylation
            new TestCase(
                description: "Common Variable:Oxidation on M and phosphorylation on S",
                mzLib: "PEPS[Common Biological:Phosphorylation on S]TM[Common Variable:Oxidation on M]IDE",
                massShift: "PEPS[+79.9663]TM[+15.9949]IDE",
                chronologer: "-PEPsTmIDE_",
                baseSeq: "PEPSTMIDE",
                residueMods: 2
            ),

            // 8. N-term acetylation + internal oxidation
            new TestCase(
                description: "N-terminal acetylation with internal oxidation",
                mzLib: "[Common Biological:Acetylation on X]PEPTM[Common Variable:Oxidation on M]IDE",
                massShift: "[+42.0106]PEPTM[+15.9949]IDE",
                chronologer: "^PEPTmIDE_",
                baseSeq: "PEPTMIDE",
                hasNTerm: true,
                residueMods: 1
            ),

            // 9. Two oxidations on different methionines
            new TestCase(
                description: "Two oxidations on different methionines",
                mzLib: "M[Common Variable:Oxidation on M]EPTM[Common Variable:Oxidation on M]IDE",
                massShift: "M[+15.9949]EPTM[+15.9949]IDE",
                chronologer: "-mEPTmIDE_",
                baseSeq: "MEPTMIDE",
                residueMods: 2
            ),
        };

        /// <summary>
        /// Edge cases and special scenarios
        /// </summary>
        public static readonly TestCase[] EdgeCases = new[]
        {
            // Single amino acid
            new TestCase(
                description: "Single unmodified residue",
                mzLib: "K",
                massShift: "K",
                chronologer: "-K_",
                baseSeq: "K"
            ),

            // Single modified residue
            new TestCase(
                description: "Single oxidized methionine",
                mzLib: "M[Common Variable:Oxidation on M]",
                massShift: "M[+15.9949]",
                chronologer: "-m_",
                baseSeq: "M",
                residueMods: 1
            ),

            // Adjacent modifications
            new TestCase(
                description: "Two adjacent oxidized methionines",
                mzLib: "PEPTM[Common Variable:Oxidation on M]M[Common Variable:Oxidation on M]IDE",
                massShift: "PEPTM[+15.9949]M[+15.9949]IDE",
                chronologer: "-PEPTmmIDE_",
                baseSeq: "PEPTMMIDE",
                residueMods: 2
            ),
        };

        /// <summary>
        /// All test cases combined
        /// </summary>
        public static TestCase[] AllTestCases => CoreTestCases.Concat(EdgeCases).ToArray();
    }
}
