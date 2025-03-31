using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ExternalResults.IndividualResultRecords
{
    public class MsFraggerSpeclib
    {
        public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        #region MsFragger Fields

        [Name("FileName")]
        public string FileName { get; set; }

        [Name("PrecursorMz")]
        public double PrecursorMz { get; set; }

        [Name("ProductMz")]
        public double ProductMz { get; set; }

        [Name("Tr_recalibrated")]
        public double Tr_recalibrated { get; set; }

        [Name("IonMobility")]
        public double IonMobility { get; set; }

        [Name("transition_name")]
        public string transition_name { get; set; }

        [Name("LibraryIntensity")]
        public double LibraryIntensity { get; set; }

        [Name("transition_group_id")]
        public string transition_group_id { get; set; }

        [Name("decoy")]
        public int decoy { get; set; }

        [Name("PeptideSequence")]
        public string PeptideSequence { get; set; }

        [Name("Proteotypic")]
        public int Proteotypic { get; set; }

        [Name("QValue")]
        public double QValue { get; set; }

        [Name("PeptidoformQValue")]
        public double PeptidoformQValue { get; set; }

        [Name("PGQValue")]
        public double PGQValue { get; set; }

        [Name("Ms1ProfileCorr")]
        public double Ms1ProfileCorr { get; set; }

        [Name("ProteinGroup")]
        public string ProteinGroup { get; set; }

        [Name("ProteinName")]
        public string ProteinName { get; set; }

        [Name("Genes")]
        public string Genes { get; set; }

        [Name("FullUniModPeptideName")]
        public string FullUniModPeptideName { get; set; }

        [Name("ModifiedPeptide")]
        public string ModifiedPeptide { get; set; }

        [Name("PrecursorCharge")]
        public int PrecursorCharge { get; set; }

        [Name("PeptideGroupLabel")]
        public string PeptideGroupLabel { get; set; }

        [Name("UniprotID")]
        public string UniprotID { get; set; }

        [Name("NTerm")]
        public int NTerm { get; set; }

        [Name("CTerm")]
        public int CTerm { get; set; }

        [Name("FragmentType")]
        public string FragmentType { get; set; }

        [Name("FragmentCharge")]
        public int FragmentCharge { get; set; }

        [Name("FragmentSeriesNumber")]
        public int FragmentSeriesNumber { get; set; }

        [Name("FragmentLossType")]
        public string FragmentLossType { get; set; }

        [Name("ExcludeFromAssay")]
        public bool ExcludeFromAssay { get; set; }

        [Name("PTM.Informative")]
        public double PTM_Informative { get; set; }

        [Name("PTM.Specific")]
        public double PTM_Specific { get; set; }

        [Name("PTM.Localising")]
        public double PTM_Localising { get; set; }

        [Name("BestPTMQValue")]
        public double BestPTMQValue { get; set; }

        [Name("PTMSiteConfidence")]
        public double PTMSiteConfidence { get; set; }

        #endregion

        #region Interpreted Fields

        //[Ignore] private string _fileNameWithoutExtension;
        //[Ignore]
        //public string FileNameWithoutExtension =>
        //    _fileNameWithoutExtension ??= Spectrum.Split('.')[0];

        //[Ignore] private int? _oneBasedScanNumber;

        //[Ignore]
        //public int OneBasedScanNumber => _oneBasedScanNumber ??= int.Parse(Spectrum.Split('.')[1]);

        //#endregion

        //#region IQuantifiableRecord Implementation

        //[Ignore] public string FileName => SpectrumFilePath;

        //[Ignore]
        //public List<(string, string, string)> ProteinGroupInfos
        //{
        //    get
        //    {
        //        _proteinGroupInfos ??= AddProteinGroupInfos();
        //        return _proteinGroupInfos;
        //    }
        //}

        ///// <summary>
        ///// Creates a list of tuples, each of which represents a protein.
        ///// Each tuple contains the accession number, gene name, and organism.
        ///// These parameters are used to create a ProteinGroup object, 
        ///// which is needed to make an identification.
        ///// </summary>
        ///// <returns></returns>
        //private List<(string, string, string)> AddProteinGroupInfos()
        //{
        //    _proteinGroupInfos = new List<(string, string, string)>();
        //    string protein = Protein;

        //    char[] delimiterChars = { '|', '_' };
        //    string[] proteinInfo = protein.Split(delimiterChars);

        //    string proteinAccessions;
        //    string geneName;
        //    string organism;

        //    // Fasta header is parsed to separate the accession number, gene name, and organism.
        //    // If the protein does not have this information, it will be assigned an empty string.
        //    // Ideally, a future refactor would create a method for parsing fasta headers
        //    // that is shared by Readers and UsefulProteomicsDatabases.
        //    proteinAccessions = proteinInfo.Length >= 2 ? proteinInfo[1] : "";
        //    geneName = proteinInfo.Length >= 3 ? proteinInfo[2] : "";
        //    organism = proteinInfo.Length >= 4 ? proteinInfo[3] : ""; ;

        //    _proteinGroupInfos.Add((proteinAccessions, geneName, organism));

        //    if (MappedProteins.IsNullOrEmpty()) return _proteinGroupInfos;

        //    string mappedProteins = MappedProteins;
        //    string[] allMappedProteinInfo = mappedProteins.Split(',');
        //    foreach (var singleMappedProteinInfo in allMappedProteinInfo)
        //    {
        //        string[] mappedProteinInfo = singleMappedProteinInfo.Split(delimiterChars);

        //        proteinAccessions = mappedProteinInfo.Length >= 2 ? mappedProteinInfo[1] : "";
        //        geneName = mappedProteinInfo.Length >= 3 ? mappedProteinInfo[2] : "";
        //        organism = mappedProteinInfo.Length >= 4 ? mappedProteinInfo[3] : "";

        //        _proteinGroupInfos.Add((proteinAccessions, geneName, organism));
        //    }

        //    return _proteinGroupInfos;
        //}

        //[Ignore] private List<(string, string, string)> _proteinGroupInfos;

        //[Ignore] public string ModifiedSequence => FullSequence.IsNullOrEmpty() ? BaseSequence : FullSequence;

        //[Ignore] public int ChargeState => Charge;

        //// decoy reading isn't currently supported for MsFragger psms, this will be revisited later
        //[Ignore] public bool IsDecoy => false;

        //[Ignore] public double MonoisotopicMass => CalculatedPeptideMass;

        #endregion
    }
}
