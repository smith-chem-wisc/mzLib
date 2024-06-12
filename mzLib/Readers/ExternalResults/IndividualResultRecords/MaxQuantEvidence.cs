using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
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

namespace Readers
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
            MissingFieldFound = null,
            HeaderValidated = null,
        };

        #region MaxQuant Fields

        [Name("Sequence")]
        public string BaseSequence { get; set; }
        [Name("Length")]
        public int PeptideLength { get; set; }
        [Name("Modifications")]
        public string Modifications { get; set; }
        [Name("Modified sequence")]
        public string ModifiedSequence { get; set; }
        [Name("Oxidation (M) Probabilities")]
        public string OxidationProbabilities { get; set; }
        [Name("Oxidation (M) Score Diffs")]
        public string OxidationScoreDiffs { get; set; }

        [Name("Acetyl (Protein N-term)")]
        public int AcetylProteinNterm { get; set; }
        
        [Name("Oxidation (M)")]
        public int Oxidation { get; set; }

        [Name("Missed cleavages")]
        public int MissedCleavages { get; set; }

        [Name("Proteins")]
        public string Proteins { get; set; }

        [Name("Leading proteins")]
        public string LeadingProteins { get; set; }

        [Name("Leading razor protein")]
        public string LeadingRazorProtein { get; set; }

        [Name("Gene names")]
        public string GeneNames { get; set; }

        [Name("Protein names")]
        public string ProteinNames { get; set; }

        [Name("Type")]
        public string Type { get; set; }

        [Name("Raw file")]
        public string FileNameWithoutExtension { get; set; }

        [Name("Experiment")]
        public string Experiment { get; set; }

        [Name("MS/MS m/z")]
        public double IsolatedMz { get; set; }

        [Name("Charge")]
        public int Charge { get; set; }

        [Name("m/z")]
        public double PrecursorMz { get; set; }

        [Name("Mass")]
        public double PrecursorMass { get; set; }

        [Name("Uncalibrated - Calibrated m/z [ppm]")]
        public double MassCalibrationDiffPpm { get; set; }

        [Name("Uncalibrated - Calibrated m/z [Da]")]
        public double MassCalibrationDiffDa { get; set; }

        [Name("Mass Error [ppm]")]
        public double MassDiffPpm { get; set; }

        [Name("Mass Error [Da]")]
        public double MassDiffDa { get; set; }

        [Name("Uncalibrated mass error [ppm]")]
        public double? UncalibratedMassDiffPpm { get; set; }

        [Name("Uncalibrated mass error [Da]")]
        public double? UncalibratedMassDiffDa { get; set; }

        [Name("Max intensity m/z 0")]
        public double LfqIntensity { get; set; }

        [Name("Retention time")]
        public double RetentionTime { get; set; } // in minutes

        [Name("Retention length")]
        public double RetentionLength { get; set; } //in minutes

        [Name("Calibrated retention time")]
        public double CalibratedRetentionTime { get; set; } // in minutes

        [Name("Calibrated retention time start")]
        public double CalibratedRetentionTimeStart { get; set; } // in minutes

        [Name("Calibrated retention time finish")]
        public double CalibratedRetentionTimeFinish { get; set; } // in minutes

        [Name("Match time difference")] 
        [Optional]
        public double? MatchTimeDifference { get; set; } // in minutes

        [Name("Match m/z difference")]
        [Optional]
        public double? MatchMzDifference { get; set; }

        [Name("Match q-value")]
        [Optional]
        public double? MatchQvalue { get; set; }

        [Name("Match score")]
        [Optional]
        public double? MatchScore { get; set; }

        [Name("Number of data points")]
        public int NumberOfDataPoints { get; set; }

        [Name("Number of scans")]
        public int NumberOfScans { get; set; }

        [Name("Number of isotopic peaks")]
        public int NumberOfIsotopicPeaks { get; set; }

        [Name("PIF")]
        public double ParentIonFraction { get; set; }

        [Name("Fraction of total spectrum")]
        public double FractionOfTotalSpectrum { get; set; }

        [Name("Base peak fraction")]
        public double BasePeakFraction { get; set; }

        [Name("PEP")]
        public double PEP { get; set; }

        [Name("MS/MS count")]
        public int PsmCount { get; set; }

        [Name("MS/MS scan number")]
        public int PrecursorScanNum { get; set; }

        [Name("MS/MS scan numbers")]
        public string PrecursorScanNumbers { get; set; }

        [Name("MS3 scan numbers")]
        public int? Ms3ScanNumbers { get; set; }

        [Name("Score")]
        public double Score { get; set; }

        [Name("Delta score")]
        public double DeltaScore { get; set; }

        [Name("Combinatorics")]
        public int Combinatorics { get; set; }

        [Name("Intensity")]
        public double SummedIntensity { get; set; }

        [Name("Reverse")]
        public char Reverse { get; set; }

        [Name("Potential contaminant")]
        public char PotentialContaminant { get; set; }

        [Name("id")]
        public int Id { get; set; }

        [Name("Protein group IDs")]
        public string ProteinGroupIds { get; set; }

        [Name("Peptide ID")]
        public int PeptideId { get; set; }

        [Name("Mod. peptide ID")]
        public int ModPeptideId { get; set; }

        [Name("Best MS/MS")]
        public int BestMsMs { get; set; } // This is an ID identifier that links to a different MaxQuant output file. It doesn't contain the scan number

        [Name("Oxidation (M) site IDs")]
        public string OxidationSiteIds { get; set; }

        [Name("Taxonomy IDs")]
        public string TaxonomyIds { get; set; }

        [Name("Taxonomy names")]
        public string TaxonomyNames { get; set; }

        [Name("Mass deficit")]
        public double MassDeficit { get; set; }

        #endregion

    }
}
