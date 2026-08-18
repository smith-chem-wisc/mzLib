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
        public string OutputDirectory { get; set; }
        public bool UseSharedPeptidesForProteinQuant { get; set; } = false;

        /// <summary>
        /// If true, the quantification engine will write raw quantification information to disk.
        /// This enables re-analysis later using different normalization or roll-up strategies without re-processing the raw data.
        /// </summary>
        public bool WriteRawInformation { get; set; } = true;

        public bool WritePeptideInformation { get; set; } = true;

        public bool WriteProteinInformation { get; set; } = true;

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
                WriteRawInformation = true,
                WritePeptideInformation = true,
                WriteProteinInformation = true
            };
        }
    }
}
