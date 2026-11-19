using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;

namespace Readers
{
    /// <summary>
    /// One row of MetaMorpheus's DIA peptide table (AllDiaPeptides.tsv): one peptide across the runs searched together.
    /// A peptide is a modified sequence; it carries the score of its best precursor in any run, and a peptide-level q-value
    /// from a fresh target-decoy competition among peptides. Column names are the same constants
    /// <see cref="DiaPrecursorFromTsv"/> uses, so a row joins the precursor table on <see cref="FileName"/>,
    /// <see cref="FullSequence"/> and <see cref="PrecursorCharge"/>. Columns are read by name, so their order in a file does
    /// not matter; the writer writes them in the order the properties are declared here.
    /// </summary>
    public record DiaPeptideFromTsv
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

        /// <summary>The peptide's sequence with modifications, in MetaMorpheus notation.</summary>
        [Name(SpectrumMatchFromTsvHeader.FullSequence)]
        public string FullSequence { get; set; } = "";

        [Name(SpectrumMatchFromTsvHeader.BaseSequence)]
        public string BaseSequence { get; set; } = "";

        /// <summary>The library's protein accessions, `|`-joined; the order means nothing.</summary>
        [Name(SpectrumMatchFromTsvHeader.ProteinAccession)]
        public string ProteinAccession { get; set; } = "";

        /// <summary>T target, D decoy, ET entrapment target, ED entrapment decoy, C contaminant.</summary>
        [Name(SpectrumMatchFromTsvHeader.DecoyContaminantTarget)]
        public string Label { get; set; } = "";

        /// <summary>The best score of the peptide's precursors over every run; higher is better.</summary>
        [Name(SpectrumMatchFromTsvHeader.Score)]
        public double Score { get; set; }

        /// <summary>Peptide q-value from target-decoy competition among peptides, across all runs searched together.</summary>
        [Name(DiaTsvHeader.QValuePeptideGlobal)]
        public double QValuePeptideGlobal { get; set; }

        /// <summary>The run of the best-scoring precursor: the spectra file's name exactly as deposited, extension included.</summary>
        [Name(SpectrumMatchFromTsvHeader.FileName)]
        public string FileName { get; set; } = "";

        /// <summary>The charge of the best-scoring precursor.</summary>
        [Name(SpectrumMatchFromTsvHeader.PrecursorCharge)]
        public int PrecursorCharge { get; set; }
    }
}
