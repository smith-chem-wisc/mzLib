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

namespace Readers
{
    public class MsFraggerPsm
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        #region MsFragger Fields

        [Name("Spectrum")]
        public string Spectrum { get; set; }

        [Name("Spectrum File")]
        public string SpectrumFilePath { get; set; }

        [Name("Peptide")]
        public string BaseSequence { get; set; }

        [Name("Modified Peptide")]
        public string ModifiedSequence { get; set; }

        [Name("Extended Peptide")]
        public string ExtendedSequence { get; set; }

        [Name("Prev AA")]
        public char PreviousAminoAcid { get; set; }

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

        [Name("Protein Start")]
        public int ProteinStart { get; set; }

        [Name("Protein End")]
        public int ProteinEnd { get; set; }

        [Name("Intensity")]
        public double Intensity { get; set; }

        [Name("Assigned Modifications")]
        public string AssignedModifications { get; set; }

        [Name("Observed Modifications")]
        public string ObservedModifications { get; set; }

        [Name("Purity")]
        public double Purity { get; set; }

        [Name("Is Unique")]
        public bool IsUnique { get; set; }

        [Name("Protein")]
        public string Protein { get; set; }

        [Name("Protein ID")]
        public string ProteinAccession { get; set; }

        [Name("Entry Name")]
        public string EntryName { get; set; }

        [Name("Gene")]
        public string Gene { get; set; }

        [Name("Protein Description")]
        public string ProteinDescription { get; set; }

        [Name("Mapped Genes")]
        public string MappedGenes { get; set; }

        [Name("Mapped Proteins")]
        public string MappedProteins { get; set; }

        #endregion

        #region Interpreted Fields

        [Ignore] private string _fileNameWithoutExtension;
        [Ignore] public string FileNameWithoutExtension =>
            _fileNameWithoutExtension ??= Spectrum.Split('.')[0];

        [Ignore] public bool IsDecoy =>
            _isDecoy ??= ProteinAccession.Contains("rev_");

        [Ignore] private bool? _isDecoy;

        [Ignore] private int? _oneBasedScanNumber;

        [Ignore]
        public int OneBasedScanNumber => _oneBasedScanNumber ??= int.Parse(Spectrum.Split('.')[1]);

        [Ignore]
        public string FullSequence => ModifiedSequence.IsNullOrEmpty() ? BaseSequence : ModifiedSequence;

        #endregion
    }
}
