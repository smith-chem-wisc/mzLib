using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    public class TopFdMzrt : IResult
    {

        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = ",",
        };

        [Name("ID")]
        public int Id { get; set; }

        [Name("Fraction_ID")]
        public int FractionId { get; set; }

        [Name("Envelope_num")]
        public int EnvelopeNumber { get; set; }

        [Name("Mass")]
        public double Mass { get; set; }

        [Name("MonoMz")]
        public double MonoisotopicMz { get; set; }

        [Name("Charge")]
        public int Charge { get; set; }

        [Name("Intensity")]
        public double Intensity { get; set; }

        [Name("mzLo")]
        public double MinimumMz { get; set; }

        [Name("mzHi")]
        public double MaximumMz { get; set; }

        [Name("rtLo")]
        public double RetentionTimeMinimum { get; set; }

        [Name("rtHi")]
        public double RetentionTimeMaximum { get; set; }

        [Name("color")]
        public string ColorHexcode { get; set; }

        [Name("opacity")]
        public double Opacity { get; set; }

        [Name("promex_score")]
        public double PromexScore { get; set; }
    }
}
