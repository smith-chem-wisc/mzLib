using System.Globalization;
using System.Text;
using Chemistry;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    public class TopFdMzrt : IHasMass
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
        public double MonoisotopicMass { get; set; }

        [Name("Charge")]
        public int Charge { get; set; }

        [Name("Intensity")]
        public double Intensity { get; set; }

        [Name("mzLo")]
        public double MzMin { get; set; }

        [Name("mzHi")]
        public double MzMax { get; set; }

        [Name("rtLo")]
        public double RetentionTimeBegin { get; set; }

        [Name("rtHi")]
        public double RetentionTimeEnd { get; set; }

        [Name("color")]
        public string ColorHexcode { get; set; }

        [Name("opacity")]
        public double Opacity { get; set; }

        [Name("promex_score")]
        public double PromexScore { get; set; }
    }
}
