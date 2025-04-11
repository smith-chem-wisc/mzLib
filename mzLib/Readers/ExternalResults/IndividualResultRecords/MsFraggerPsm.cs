using System.Globalization;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using Easy.Common.Extensions;

namespace Readers
{
    public class MsFraggerPsm : IQuantifiableRecord
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
        public string ModifiedPeptide { get; set; }

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

        /// <summary>
        /// MsFragger v22.0 output renames the header "PeptideProphet Probability" as just "Probability".
        /// Headers are mutually exclusive, will not both occur in the same file. 
        /// </summary>
        [Name("PeptideProphet Probability", "Probability")]
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

        [Ignore] private int? _oneBasedScanNumber;

        [Ignore]
        public int OneBasedScanNumber => _oneBasedScanNumber ??= int.Parse(Spectrum.Split('.')[1]);

        #endregion

        #region IQuantifiableRecord Implementation

        [Ignore] public string FileName => SpectrumFilePath;

        [Ignore] public List<(string, string, string)> ProteinGroupInfos
        {
            get 
            {
                _proteinGroupInfos ??= AddProteinGroupInfos();
                return _proteinGroupInfos;
            }
        }

        [Ignore] public string FullSequence => ModifiedPeptide.IsNullOrEmpty() ? BaseSequence : ModifiedPeptide;

        /// <summary>
        /// Creates a list of tuples, each of which represents a protein.
        /// Each tuple contains the accession number, gene name, and organism.
        /// These parameters are used to create a ProteinGroup object, 
        /// which is needed to make an identification.
        /// </summary>
        /// <returns></returns>
        private List<(string, string, string)> AddProteinGroupInfos ()
        {
            _proteinGroupInfos = new List<(string, string, string)> ();
            string protein = Protein;

            char[] delimiterChars = { '|', '_'};
            string[] proteinInfo = protein.Split(delimiterChars);

            string proteinAccessions;
            string geneName;
            string organism;

            // Fasta header is parsed to separate the accession number, gene name, and organism.
            // If the protein does not have this information, it will be assigned an empty string.
            // Ideally, a future refactor would create a method for parsing fasta headers
            // that is shared by Readers and UsefulProteomicsDatabases.
            proteinAccessions = proteinInfo.Length >= 2 ? proteinInfo[1] : "";
            geneName = proteinInfo.Length >= 3 ? proteinInfo[2] : "";
            organism = proteinInfo.Length >= 4 ? proteinInfo[3] : ""; ;

            _proteinGroupInfos.Add((proteinAccessions, geneName, organism));

            if (MappedProteins.IsNullOrEmpty()) return _proteinGroupInfos;

            string mappedProteins = MappedProteins;
            string[] allMappedProteinInfo = mappedProteins.Split(',');
            foreach (var singleMappedProteinInfo in allMappedProteinInfo)
            {
                string[] mappedProteinInfo = singleMappedProteinInfo.Split(delimiterChars);

                proteinAccessions = mappedProteinInfo.Length >= 2 ? mappedProteinInfo[1] : "";
                geneName = mappedProteinInfo.Length >= 3 ? mappedProteinInfo[2] : "";
                organism = mappedProteinInfo.Length >= 4 ? mappedProteinInfo[3] : "";

                _proteinGroupInfos.Add((proteinAccessions, geneName, organism));
            }

            return _proteinGroupInfos;
        }

        [Ignore] private List<(string, string, string)> _proteinGroupInfos;

        

        [Ignore] public int ChargeState => Charge;

        // decoy reading isn't currently supported for MsFragger psms, this will be revisited later
        [Ignore] public bool IsDecoy => false;

        [Ignore] public double MonoisotopicMass => CalculatedPeptideMass;

        #endregion
    }
}