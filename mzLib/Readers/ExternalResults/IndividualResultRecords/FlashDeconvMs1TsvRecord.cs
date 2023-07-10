using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System.Globalization;
using System.Text;
using Chemistry;

namespace Readers
{
    /// <summary>
    /// A class representing a single entry in FlashDeconv's ms1.tsv output file
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class FlashDeconvMs1Tsv : IHasMass
    {
        [Ignore]
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t",
        };

        [Name("Index")]
        public int Index { get; set; }

        [Name("FileName")]
        public string FileName { get; set; }

        [Name("ScanNum")]
        public int ZeroBasedScanNumber { get; set; }

        [Ignore]
        public int OneBasedScanNumber => ZeroBasedScanNumber + 1;

        /// <summary>
        /// Target Decoy of FlashDeconv
        ///  0 -> Target
        ///  1 -> Isotope Decoy
        ///  2 -> Noise Decoy
        ///  3 -> Charge Decoy
        /// </summary>
        [Name("Decoy")]
        public int Decoy { get; set; }

        [Name("RetentionTime")]
        public double RetentionTime { get; set; }

        /// <summary>
        /// Number of m/z in the spectrum after deconvolution
        /// </summary>
        [Name("MassCountInSpec")]
        public int MassCountInSpec { get; set; }

        [Name("AverageMass")]
        public double AverageMass { get; set; }

        [Name("MonoisotopicMass")]
        public double MonoisotopicMass { get; set; }

        [Name("SumIntensity")]
        public double SumIntensity { get; set; }

        /// <summary>
        /// Min charge of all charge states identified
        /// </summary>
        [Name("MinCharge")]
        public int MinCharge { get; set; }

        /// <summary>
        /// Max charge of all charge states identified
        /// </summary>
        [Name("MaxCharge")]
        public int MaxCharge { get; set; }

        /// <summary>
        /// The number of isotopic peaks from all detected charge states 
        /// </summary>
        [Name("PeakCount")]
        public int PeakCount { get; set; }

        [Name("IsotopeCosine")]
        public double IsotopeCosine { get; set; }

        [Name("ChargeScore")]
        public double ChargeScore { get; set; }

        /// <summary>
        /// Calculation can be found here
        /// https://www.nature.com/articles/s41467-022-31922-z
        /// </summary>
        [Name("MassSNR")]
        public double MassSNR { get; set; }

        /// <summary>
        /// Calculation can be found here
        /// https://www.nature.com/articles/s41467-022-31922-z
        /// </summary>
        [Name("ChargeSNR")]
        public double ChargeSNR { get; set; }

        /// <summary>
        /// Charge of the entries representative peak
        /// </summary>
        [Name("RepresentativeCharge")]
        public int RepresentativeCharge { get; set; }

        /// <summary>
        /// Start mz of the entries representative peak
        /// </summary>
        [Name("RepresentativeMzStart")]
        public double RepresentativeMzMin { get; set; }

        /// <summary>
        /// End mz of the entries representative peak
        /// </summary>
        [Name("RepresentativeMzEnd")]
        public double RepresentativeMzMax { get; set; }

        [Name("QScore")]
        public double QScore { get; set; }

        /// <summary>
        /// IMPORTANT: Decoys will always have their Q value set to 1
        /// </summary>
        [Name("Qvalue")]
        public double Qvalue { get; set; }

        [Name("QvalueWithIsotopeDecoyOnly")]
        public double QvalueWithIsotopeDecoyOnly { get; set; }

        [Name("QvalueWithNoiseDecoyOnly")]
        public double QvalueWithNoiseDecoyOnly { get; set; }

        [Name("QvalueWithChargeDecoyOnly")]
        public double QvalueWithChargeDecoyOnly { get; set; }

        /// <summary>
        /// Intensity of each charge state in the charge state envelope
        /// If a charge state is not found, value will be 0
        /// </summary>
        [Name("PerChargeIntensity")]
        [TypeConverter(typeof(SemicolonDelimitedToDoubleListConverter))]
        public List<double> PerChargeIntensity { get; set; }

        /// <summary>
        /// Intensity of each peak in the representative charge state
        /// </summary>
        [Name("PerIsotopeIntensity")]
        [TypeConverter(typeof(SemicolonDelimitedToDoubleListConverter))]
        public List<double> PerIsotopeIntensity { get; set; }
    }
}
