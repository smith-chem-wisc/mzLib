using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One row of MetaMorpheus's DIA precursor table (AllDiaPrecursors.tsv): one library precursor identified in one run.
    /// The schema was agreed with dataRepo (dia 003). Column names reuse MetaMorpheus's psmtsv headers where the concept
    /// exists; retention columns carry their unit, and iRT and minutes are never mixed in one column.
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
        [Name("File Name"), Index(0)]
        public string FileName { get; set; } = "";

        /// <summary>The precursor's sequence with modifications, in MetaMorpheus notation.</summary>
        [Name("Full Sequence"), Index(1)]
        public string FullSequence { get; set; } = "";

        [Name("Base Sequence"), Index(2)]
        public string BaseSequence { get; set; } = "";

        [Name("Precursor Charge"), Index(3)]
        public int PrecursorCharge { get; set; }

        /// <summary>The library precursor's m/z.</summary>
        [Name("Precursor MZ"), Index(4)]
        public double PrecursorMz { get; set; }

        /// <summary>The library's protein accessions, `|`-joined; the order means nothing.</summary>
        [Name("Protein Accession"), Index(5)]
        public string ProteinAccession { get; set; } = "";

        /// <summary>T target, D decoy, ET entrapment target, ED entrapment decoy, C contaminant.</summary>
        [Name("Decoy/Contaminant/Target"), Index(6)]
        public string Label { get; set; } = "";

        /// <summary>The classifier score the q-values are computed on; higher is better.</summary>
        [Name("Score"), Index(7)]
        public double Score { get; set; }

        /// <summary>Precursor q-value from target-decoy competition within this run.</summary>
        [Name("QValue_Precursor_Run"), Index(8)]
        public double QValuePrecursorRun { get; set; }

        /// <summary>Precursor q-value from target-decoy competition across all runs searched together.</summary>
        [Name("QValue_Precursor_Global"), Index(9)]
        public double QValuePrecursorGlobal { get; set; }

        /// <summary>The library's iRT for the precursor.</summary>
        [Name("LibraryIrt"), Index(10)]
        public double LibraryIrt { get; set; }

        /// <summary>The apex's retention time on the library's iRT scale, under this run's calibration.</summary>
        [Name("ApexIrt"), Index(11)]
        public double ApexIrt { get; set; }

        /// <summary>The apex's retention time in this run, minutes.</summary>
        [Name("ApexRtMin"), Index(12)]
        public double ApexRtMin { get; set; }

        /// <summary>The one-based scan number of the apex MS2 scan, so a row can point at a spectrum through a USI.</summary>
        [Name("ApexScanNumber"), Index(13)]
        public int ApexScanNumber { get; set; }

        /// <summary>
        /// The precursor's quantity: its fragments' summed peak areas. Empty when not reported; a 0 is an integrated area
        /// of exactly 0 over a detected peak.
        /// </summary>
        [Name("Precursor Quantity"), Index(14)]
        public double? PrecursorQuantity { get; set; }
    }
}
