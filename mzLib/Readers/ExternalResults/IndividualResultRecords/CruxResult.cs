using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MzLibUtil;

namespace Readers
{
    public class CruxResult
    { 
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t",
        };

        [Name("file")]
        public string FilePath { get; set; }

        [Name("scan")]
        public int OneBasedScanNumber { get; set; }

        [Name("charge")]
        public int Charge { get; set; }

        [Name("retention time")]
        public double RetentionTime { get; set; }

        [Name("spectrum precursor m/z")]
        public double PrecursorMz { get; set; }

        [Name("spectrum neutral mass")]
        public double NeutralMass { get; set; }

        [Name("peptide mass")]
        public double PeptideMass { get; set; }

        [Name("delta_cn")]
        public double DeltaCn { get; set; }

        [Name("xcorr score")]
        public double XCorrScore { get; set; }

        [Name("xcorr rank")]
        public int XCorrRank { get; set; }

        [Name("tailor score")]
        public double TailorScore { get; set; }

        [Name("tdc q-value")]
        public double TdcQValue { get; set; }

        [Name("b/y ions matched")]
        public int BAndYIonsMatched { get; set; }

        [Name("b/y ions total")]
        public int BAndYIonsTotal { get; set; }

        [Name("b/y ions fraction")]
        public double BAndYIonsFraction { get; set; }

        [Name("b/y ion repeat match")]
        public int BAndYIonRepeatMatch { get; set; }

        [Name("distinct matches/spectrum")]
        public int DistinctMatchesPerSpectrum { get; set; }

        [Name("sequence")]
        public string FullSequence { get; set; }

        [Name("unmodified sequence")]
        public string BaseSequence { get; set; }

        [Name("protein id")]
        public string ProteinId { get; set; }

        [Name("flanking aa")]
        public string FlankingAa { get; set; }

        #region Interpreted properties

        [Ignore] private string? _fileNameWithoutExtension = null;
        [Ignore] public string FileNameWithoutExtension => _fileNameWithoutExtension ??= FilePath.GetPeriodTolerantFilenameWithoutExtension();

        [Ignore] private string? _accession = null;
        [Ignore] public string Accession => _accession ??= ProteinId.Split('|')[1].Trim();

        #endregion
    }
}
