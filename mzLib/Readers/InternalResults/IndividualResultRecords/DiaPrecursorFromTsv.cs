using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MassSpectrometry;

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

        private double? _precursorQuantity;
        private double? _ms1Quantity;

        /// <summary>
        /// The precursor's quantity (QuantProject's DEF-DIA-PREC): its selected fragments' summed peak areas, MS2 only
        /// (<see cref="FragmentQuantity"/>). Empty for every kind of "nothing", never 0, with the kind in
        /// <see cref="QuantityState"/>. A 0, set here or read from a table an older search wrote (where it meant that nothing
        /// was measured), is no value.
        /// </summary>
        [Name(DiaTsvHeader.PrecursorQuantity)]
        public double? PrecursorQuantity
        {
            get => _precursorQuantity;
            set => _precursorQuantity = value == 0 ? null : value;
        }

        /// <summary>
        /// What <see cref="PrecursorQuantity"/> is (DEF-DIA-KIND): Value, BelowFloor or NotLooked, written by name. Empty in
        /// a table written before the column existed.
        /// </summary>
        [Name(DiaTsvHeader.QuantityState)]
        [Optional]
        public PrecursorQuantityKind? QuantityState { get; set; }

        /// <summary>The floor <see cref="PrecursorQuantity"/>'s fragments were judged against, in its units.</summary>
        [Name(DiaTsvHeader.QuantityFloor)]
        [Optional]
        public double? QuantityFloor { get; set; }

        /// <summary>
        /// The precursor's MS1 quantity (DEF-DIA-MS1). It is never added to, averaged with or substituted for
        /// <see cref="PrecursorQuantity"/>. Empty for every kind of "nothing", never 0, with the kind in
        /// <see cref="Ms1QuantityState"/>.
        /// </summary>
        [Name(DiaTsvHeader.Ms1Quantity)]
        [Optional]
        public double? Ms1Quantity
        {
            get => _ms1Quantity;
            set => _ms1Quantity = value == 0 ? null : value;
        }

        /// <summary>What <see cref="Ms1Quantity"/> is (DEF-DIA-KIND), written by name.</summary>
        [Name(DiaTsvHeader.Ms1QuantityState)]
        [Optional]
        public PrecursorQuantityKind? Ms1QuantityState { get; set; }

        /// <summary>The floor <see cref="Ms1Quantity"/> was judged against, in its units.</summary>
        [Name(DiaTsvHeader.Ms1QuantityFloor)]
        [Optional]
        public double? Ms1QuantityFloor { get; set; }
    }
}
