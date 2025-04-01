using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.SpectrumLibraries
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

    }
}
