using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MassSpectrometry;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using Proteomics;
using static System.Net.Mime.MediaTypeNames;
using ThermoFisher.CommonCore.Data.Interfaces;
using Easy.Common.Extensions;

namespace Readers.ExternalResults.IndividualResultRecords
{
    public class MaxQuantEvidence
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        #region MaxQuant Fields

        [Name("Sequence")]
        public string BaseSequence { get; set; }

        [Name("Modified sequence")]
        public string ModifiedSequence { get; set; }

        [Name("Proteins")]
        public string Proteins { get; set; }

        [Name("Leading proteins")]
        public char LeadingProteins { get; set; }

        [Name("Next AA")]
        public char NextAminoAcid { get; set; }

        [Name("Peptide Length")]
        public int PeptideLength { get; set; }

        [Name("Charge")]
        public int Charge { get; set; }

        [Name("Retention")]
        public double RetentionTime { get; set; }

        [Name("Observed Mass")]
        public double ObservedMass { get; set; }

        [Name("Calibrated Observed Mass")]
        public double CalibratedObservedMass { get; set; }

        [Name("Observed M/Z")]
        public double ObservedMz { get; set; }

        [Name("Calibrated Observed M/Z")]
        public double CalibratedObservedMz { get; set; }

        [Name("Calculated Peptide Mass")]
        public double CalculatedPeptideMass { get; set; }

        [Name("Calculated M/Z")]
        public double CalculatedMz { get; set; }

        [Name("Delta Mass")]
        public double DeltaMass { get; set; }

        [Name("Expectation")]
        public double Expectation { get; set; }

        [Name("Hyperscore")]
        public double HyperScore { get; set; }

        [Name("Nextscore")]
        public double NextScore { get; set; }

        [Name("PeptideProphet Probability")]
        public double PeptideProphetProbability { get; set; }

        [Name("Number of Enzymatic Termini")]
        public int NumberOfEnzymaticTermini { get; set; }

        [Name("Number of Missed Cleavages")]
        public int NumberOfMissedCleavages { get; set; }
        #endregion

    }
}
