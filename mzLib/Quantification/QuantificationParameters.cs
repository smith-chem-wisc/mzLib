using Quantification.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Quantification
{
    public class QuantificationParameters
    {
        public INormalizationStrategy SpectralMatchNormalizationStrategy { get; set; }
        public INormalizationStrategy PeptideNormalizationStrategy { get; set; }
        public INormalizationStrategy ProteinNormalizationStrategy { get; set; }
        public IRollUpStrategy SpectralMatchToPeptideRollUpStrategy { get; set; }
        public IRollUpStrategy PeptideToProteinRollUpStrategy { get; set; }
        public string OutputDirectory { get; set; }
        public bool UseSharedPeptidesForProteinQuant { get; set; } = false;

        /// <summary>
        /// If true, the quantification engine will write raw quantification information to disk.
        /// This enables re-analysis later using different normalization or roll-up strategies without re-processing the raw data.
        /// </summary>
        public bool WriteRawInformation { get; set; } = true;

        public bool WritePeptideInformation { get; set; } = true;

        public bool WriteProteinInformation { get; set; } = true;
    }
}
