using System.Globalization;
using System.Text;
using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;

namespace Readers
{
    /// <summary>
    /// A class representing a single entry in FlashDeconv's .tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvTsv
    {
        [Ignore]
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t",
        };

        [Name("FeatureIndex")]
        public int FeatureIndex { get; set; }

        [Name("FileName")]
        public string FilePath { get; set; }

        [Name("MonoisotopicMass")]
        public double MonoisotopicMass { get; set; }

        [Name("AverageMass")]
        public double AverageMass { get; set; }

        [Name("MassCount")]
        public double MassCount { get; set; }

        [Name("StartRetentionTime")]
        public double RetentionTimeBegin { get; set; }

        [Name("EndRetentionTime")]
        public double RetentionTimeEnd { get; set; }

        [Name("RetentionTimeDuration")]
        public double RetentionTimeDuration { get; set; }

        [Name("ApexRetentionTime")]
        public double RetentionTimeApex { get; set; }

        [Name("SumIntensity")]
        public double IntensityOfAllPeaks { get; set; }

        [Name("MaxIntensity")]
        public double IntensityOfMostAbundantPeak { get; set; }

        [Name("FeatureQuantity")]
        public double FeatureQuantity { get; set; }

        [Name("MinCharge")]
        public int ChargeStateMin { get; set; }

        [Name("MaxCharge")]
        public int ChargeStateMax { get; set; }

        [Name("ChargeCount")]
        public int ChargeCount { get; set; }

        [Name("IsotopeCosineScore")]
        public double IsotopeCosineSimilarity { get; set; }

        [Name("MaxQScore")]
        public double MaxQScore { get; set; }

        [Name("PerChargeIntensity")]
        [TypeConverter(typeof(SemicolonDelimitedToDoubleListConverter))]
        public List<double> IntensityPerChargeState { get; set; }

        [Name("PerIsotopeIntensity")]
        [TypeConverter(typeof(SemicolonDelimitedToDoubleListConverter))]
        public List<double> IntensityPerIsotope { get; set; }
    }
}
