using System.Globalization;
using System.Text;
using Chemistry;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MassSpectrometry;

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

        // Sample_ID is an int in the original FlashDeconv / TopFD-v1.6.2 _ms1.feature schema.
        // Newer TopFD output uses "File_name" (a string path) in its place and omits Sample_ID
        // entirely. Sample_ID isn't read by GetSingleChargeFeatures, so mark it Optional rather
        // than try to coerce a path into an int.
        [Name("Sample_ID")]
        [Optional]
        public int SampleId { get; set; }

        [Name("ID", "Feature_ID")]
        [Optional]
        public int Id { get; set; }

        [Name("Mass")]
        public double Mass { get; set; }

        [Name("Intensity")]
        public double Intensity { get; set; }

        [Name("Time_begin", "Min_time")]
        public double RetentionTimeBegin { get; set; }

        [Name("Time_end", "Max_time")]
        public double RetentionTimeEnd { get; set; }

        [Name("Time_apex", "Apex_time")]
        public double RetentionTimeApex { get; set; }

        [Name("Apex_intensity", "Intensity_Apex")]
        [Optional]
        public double? IntensityApex { get; set; }

        [Name("Minimum_charge_state", "Min_charge")]
        public int ChargeStateMin { get; set; }

        [Name("Maximum_charge_state", "Max_charge")]
        public int ChargeStateMax { get; set; }

        // Newer TopFD output has a single "Fraction_ID" column rather than separate min/max
        // fraction-id bounds. These fields aren't consumed downstream by GetSingleChargeFeatures,
        // so making them Optional keeps both schemas readable without forcing an alias collision.
        [Name("Minimum_fraction_id")]
        [Optional]
        public int FractionIdMin { get; set; }

        [Name("Maximum_fraction_id")]
        [Optional]
        public int FractionIdMax { get; set; }

        /// <summary>
        /// Expands this row into one <see cref="ISingleChargeMs1Feature"/> per charge in
        /// [<see cref="ChargeStateMin"/>, <see cref="ChargeStateMax"/>]. The per-charge
        /// <c>Intensity</c> is taken from <see cref="IntensityApex"/> (the apex intensity),
        /// not the summed <see cref="Intensity"/> column: the summed value is preserved on
        /// the record for fidelity, but downstream MS2 pairing weights by the charge-state
        /// apex. Changing this to use the summed column would also change how external
        /// TopFD / FLASHDeconv files are interpreted.
        /// </summary>
        public IEnumerable<ISingleChargeMs1Feature> GetSingleChargeFeatures()
        {
            for (int z = ChargeStateMin; z <= ChargeStateMax ; z++)
            {
                yield return new SingleChargeMs1Feature(Mass.ToMz(z), z, RetentionTimeBegin, RetentionTimeEnd, IntensityApex ?? 0);
            }
        }
    }

    public record SingleChargeMs1Feature(double Mz, int Charge, double RetentionTimeStart, double RetentionTimeEnd, double Intensity,
        int? NumberOfIsotopes = null) : ISingleChargeMs1Feature;
}
