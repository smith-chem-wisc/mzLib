using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// A class representing a single entry in a .feature.tsv file, the output from Dinosaur
    /// </summary>
    public class DinosaurFeature : ISingleChargeMs1Feature
    {
        #region IMs1FeatureProperties

        public double RetentionTimeStart => RtStart;
        public double RetentionTimeEnd => RtEnd;
        public double Intensity => IntensityApex;
        public int? NumberOfIsotopes => NIsotopes;

        #endregion

        [Ignore]
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t"
        };

        [Name("mz")]
        public double Mz { get; set; }

        [Name("mostAbundantMz")]
        public double MostAbundantMz { get; set; }

        [Name("charge")]
        public int Charge { get; set; }

        [Name("rtStart")]
        public double RtStart { get; set; }

        [Name("rtApex")]
        public double RtApex { get; set; }

        [Name("rtEnd")]
        public double RtEnd { get; set; }

        [Name("fwhm")]
        public double Fwhm { get; set; }

        [Name("nIsotopes")]
        public int NIsotopes { get; set; }

        [Name("nScans")]
        public int NScans { get; set; }

        [Name("averagineCorr")]
        public double AveragineCorr { get; set; }

        [Name("mass")]
        public double Mass { get; set; }

        [Name("massCalib")]
        public double MassCalib { get; set; }

        [Name("intensityApex")]
        public double IntensityApex { get; set; }

        [Name("intensitySum")]
        public double IntensitySum { get; set; }
    }
}
