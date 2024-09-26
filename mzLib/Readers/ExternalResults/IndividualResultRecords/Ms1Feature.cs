using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// A class representing a single entry in a ms1.feature file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class Ms1Feature
    {
        [Ignore]
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        [Name("Sample_ID")]
        public int SampleId { get; set; }

        [Name("ID")]
        public int Id { get; set; }

        [Name("Mass")]
        public double Mass { get; set; }

        [Name("Intensity")]
        public double Intensity { get; set; }

        [Name("Time_begin")]
        public double RetentionTimeBegin { get; set; }

        [Name("Time_end")]
        public double RetentionTimeEnd { get; set; }

        [Name("Time_apex", "Apex_time")]
        public double RetentionTimeApex { get; set; }

        [Name("Apex_intensity", "Intensity_Apex")]
        [Optional]
        public double? IntensityApex { get; set; }

        [Name("Minimum_charge_state")]
        public int ChargeStateMin { get; set; }

        [Name("Maximum_charge_state")]
        public int ChargeStateMax { get; set; }

        [Name("Minimum_fraction_id")]
        public int FractionIdMin { get; set; }

        [Name("Maximum_fraction_id")]
        public int FractionIdMax { get; set; }
    }
}
