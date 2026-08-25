using Quantification.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Quantification.Strategies;

namespace Quantification
{
    /// <summary>
    /// QuantificationParameters contain all the strategies and settings needed to perform quantification. These are:
    /// 3 normalization strategies (for spectral matches, peptides, and proteins),
    /// 2 roll-up strategies (spectral matches to peptides, peptides to proteins),
    /// 1 collapse strategy (to collapse samples, e.g., fractions and technical replicates)
    /// Quantification is performed in the following order:
    /// 1) Spectral Match Normalization ->
    /// 2) Spectral Match to Peptide Roll-Up ->
    /// 3) Peptide Normalization ->
    /// 4) Collapse Samples ->
    /// 5) Peptide to Protein Roll-Up ->
    /// 6) Protein Normalization
    /// The strategies are listed in the order they are applied.
    /// The parameters also contain settings for output directory and whether to use shared peptides for protein quantification.
    /// </summary>
    public class QuantificationParameters
    {
        public INormalizationStrategy SpectralMatchNormalizationStrategy { get; set; }
        public IRollUpStrategy SpectralMatchToPeptideRollUpStrategy { get; set; }
        public INormalizationStrategy PeptideNormalizationStrategy { get; set; }
        public ICollapseStrategy CollapseStrategy { get; set; }
        public IRollUpStrategy PeptideToProteinRollUpStrategy { get; set; }
        public INormalizationStrategy ProteinNormalizationStrategy { get; set; }
        /// <summary>
        /// Where the output files are written. Optional: when this is left empty and any of the
        /// <c>Write*Information</c> flags is on, the engine writes beside the source data — the directory
        /// holding the spectral matches' files, or their nearest common ancestor when they span several.
        /// Setting it always wins over that default, and a run that writes nothing ignores it entirely.
        /// The directory actually used is reported as <see cref="QuantificationResults.OutputDirectory"/>.
        /// </summary>
        public string OutputDirectory { get; set; }
        public bool UseSharedPeptidesForProteinQuant { get; set; } = false;

        /// <summary>
        /// If true, the engine writes the PSM-level intensities to <see cref="OutputDirectory"/>
        /// exactly as they arrived, before any normalization or roll-up.
        ///
        /// On by default, because this file is what makes re-processing possible: with it, the same
        /// search can be re-quantified under different normalization, roll-up and collapse strategies
        /// without re-searching. Turning it off discards that option for the run.
        ///
        /// <see cref="OutputDirectory"/> need not be set: the engine falls back to the directory the
        /// source files came from. It rejects the run with a descriptive
        /// <see cref="QuantificationResults.Summary"/> only when neither is available, rather than
        /// writing into the working directory.
        /// </summary>
        public bool WriteRawInformation { get; set; } = true;

        /// <summary>
        /// If true, the final peptide matrix is written to <see cref="OutputDirectory"/>. On by default.
        /// </summary>
        public bool WritePeptideInformation { get; set; } = true;

        /// <summary>
        /// If true, the final protein-group matrix is written to <see cref="OutputDirectory"/>. On by default.
        /// </summary>
        public bool WriteProteinInformation { get; set; } = true;

        /// <summary>
        /// A no-op configuration — no normalization, no collapsing, sum roll-ups — for the test and
        /// development projects only. It is deliberately <c>internal</c>: production callers should
        /// construct a <see cref="QuantificationParameters"/> and choose strategies for their data.
        ///
        /// Writing is switched off so that a test run leaves nothing behind. That matters more than it
        /// looks: with no <see cref="OutputDirectory"/>, a write-enabled run would fall back to writing
        /// beside whatever data files the test pointed at.
        /// </summary>
        internal static QuantificationParameters GetSimpleParameters()
        {
            return new QuantificationParameters
            {
                SpectralMatchNormalizationStrategy = new NoNormalization(),
                SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
                PeptideNormalizationStrategy = new NoNormalization(),
                CollapseStrategy = new NoCollapse(),
                PeptideToProteinRollUpStrategy = new SumRollUp(),
                ProteinNormalizationStrategy = new NoNormalization(),
                OutputDirectory = string.Empty,
                UseSharedPeptidesForProteinQuant = false,
                WriteRawInformation = false,
                WritePeptideInformation = false,
                WriteProteinInformation = false
            };
        }
    }
}
