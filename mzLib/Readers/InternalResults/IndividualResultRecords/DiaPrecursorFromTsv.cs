using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One row of MetaMorpheus's DIA precursor table (AllDiaPrecursors.tsv): one library precursor identified in one run.
    /// The schema was agreed with dataRepo (dia 003). Column names reuse MetaMorpheus's psmtsv headers where the concept
    /// exists (<see cref="SpectrumMatchFromTsvHeader"/>), and the rest are <see cref="DiaTsvHeader"/>; retention columns
    /// carry their unit, and iRT and minutes are never mixed in one column. Columns are read by name, so their order in a
    /// file does not matter; the writer writes them in the order the properties are declared here.
    /// </summary>
    public record DiaPrecursorFromTsv
    {
        [Ignore]
        public static CsvConfiguration CsvConfiguration => new(CultureInfo.InvariantCulture)
        {
            Encoding = Encoding.UTF8,
            HasHeaderRecord = true,
            Delimiter = "\t",
            AllowComments = true,
            Comment = '#',
        };

        /// <summary>The run key: the spectra file's name exactly as deposited, extension included, never shortened.</summary>
        [Name(SpectrumMatchFromTsvHeader.FileName)]
        public string FileName { get; set; } = "";

        /// <summary>The precursor's sequence with modifications, in MetaMorpheus notation.</summary>
        [Name(SpectrumMatchFromTsvHeader.FullSequence)]
        public string FullSequence { get; set; } = "";

        [Name(SpectrumMatchFromTsvHeader.BaseSequence)]
        public string BaseSequence { get; set; } = "";

        [Name(SpectrumMatchFromTsvHeader.PrecursorCharge)]
        public int PrecursorCharge { get; set; }

        /// <summary>The library precursor's m/z.</summary>
        [Name(SpectrumMatchFromTsvHeader.PrecursorMz)]
        public double PrecursorMz { get; set; }

        /// <summary>The library's protein accessions, `|`-joined; the order means nothing.</summary>
        [Name(SpectrumMatchFromTsvHeader.ProteinAccession)]
        public string ProteinAccession { get; set; } = "";

        /// <summary>T target, D decoy, ET entrapment target, ED entrapment decoy, C contaminant.</summary>
        [Name(SpectrumMatchFromTsvHeader.DecoyContaminantTarget)]
        public string Label { get; set; } = "";

        /// <summary>The classifier score the q-values are computed on; higher is better.</summary>
        [Name(SpectrumMatchFromTsvHeader.Score)]
        public double Score { get; set; }

        /// <summary>Precursor q-value from target-decoy competition within this run.</summary>
        [Name(DiaTsvHeader.QValuePrecursorRun)]
        public double QValuePrecursorRun { get; set; }

        /// <summary>Precursor q-value from target-decoy competition across all runs searched together.</summary>
        [Name(DiaTsvHeader.QValuePrecursorGlobal)]
        public double QValuePrecursorGlobal { get; set; }

        /// <summary>The library's iRT for the precursor.</summary>
        [Name(DiaTsvHeader.LibraryIrt)]
        public double LibraryIrt { get; set; }

        /// <summary>The apex's retention time on the library's iRT scale, under this run's calibration.</summary>
        [Name(DiaTsvHeader.ApexIrt)]
        public double ApexIrt { get; set; }

        /// <summary>The apex's retention time in this run, minutes.</summary>
        [Name(DiaTsvHeader.ApexRtMin)]
        public double ApexRtMin { get; set; }

        /// <summary>The one-based scan number of the apex MS2 scan, so a row can point at a spectrum through a USI.</summary>
        [Name(DiaTsvHeader.ApexScanNumber)]
        public int ApexScanNumber { get; set; }

        /// <summary>
        /// The precursor's quantity: its fragments' summed peak areas. Empty when not reported; a 0 is an integrated area
        /// of exactly 0 over a detected peak.
        /// </summary>
        [Name(DiaTsvHeader.PrecursorQuantity)]
        public double? PrecursorQuantity { get; set; }
    }
}
