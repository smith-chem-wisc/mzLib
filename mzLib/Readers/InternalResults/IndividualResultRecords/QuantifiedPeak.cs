using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    public class QuantifiedPeak
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = CsvHelper.Configuration.TrimOptions.Trim,
        };

        [Name("File Name")]
        public string FileName { get; set; }

        [Name("Base Sequence")]
        public string BaseSequence { get; set; }

        [Name("Full Sequence")]
        public string FullSequence { get; set; }

        [Name("Protein Group")]
        public string ProteinGroup { get; set; }

        [Name("Peptide Monoisotopic Mass")]
        public double PeptideMonoisotopicMass { get; set; }

        [Name("MS2 Retention Time")]
        public double? MS2RetentionTime { get; set; }

        [Name("Precursor Charge")]
        public int PrecursorCharge { get; set; }

        [Name("Theoretical MZ")]
        public double TheoreticalMZ { get; set; }

        [Name("Peak intensity")]
        public double PeakIntensity { get; set;}

        [Name("Peak RT Start")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double? PeakRTStart { get; set; }

        [Name("Peak RT Apex")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double? PeakRTApex { get; set; }

        [Name("Peak RT End")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double? PeakRTEnd { get; set; }

        [Name("Peak MZ")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double? PeakMz { get; set; }

        [Name("Peak Charge")]
        [TypeConverter(typeof(DashToNullOrIntegerConverter))]
        public int? PeakCharge { get; set; }

        [Name("Num Charge States Observed")]
        public int NumChargeStatesObserved { get; set; }

        [Name("Peak Detection Type")]
        public string PeakDetectionType { get; set; }

        [Name("MBR Score")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double MBRScore { get; set; }

        [Name("PSMs Mapped")]
        public int PSMsMapped { get; set; }

        [Name("Base Sequences Mapped")]
        public int BaseSequencesMapped { get; set; }

        [Name("Full Sequences Mapped")]
        public int FullSequencesMapped { get; set; }

        [Name("Peak Split Valley RT")]
        public double PeakSplitValleyRT { get; set; }

        [Name("Peak Apex Mass Error (ppm)")]
        [TypeConverter(typeof(DashToNullOrDoubleConverter))]
        public double? PeakApexMassError { get; set; }
    }
}
