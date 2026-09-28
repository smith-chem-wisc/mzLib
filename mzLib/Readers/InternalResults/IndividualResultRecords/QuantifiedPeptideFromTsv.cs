using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System.Globalization;

namespace Readers;

/// <summary>
/// One row of a FlashLFQ peptide table: <c>QuantifiedPeptides.tsv</c> from FlashLFQ, or
/// <c>AllQuantifiedPeptides.tsv</c> from MetaMorpheus. The writer is <c>FlashLFQ.Peptide.ToString</c>
/// under <c>Peptide.TabSeparatedHeader</c>.
/// </summary>
/// <remarks>
/// Unlike the protein-group table, a peptide that was not quantified in a sample is written as a
/// literal <c>0</c>, not a blank. The reader keeps what was written: <c>0</c> reads as 0, and only a
/// blank reads as null. Use <see cref="QuantifiedPeptideSample.DetectionType"/> to tell "not
/// detected" from a measured value.
/// </remarks>
public class QuantifiedPeptideFromTsv
{
    public static CsvConfiguration CsvConfiguration => new(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        IgnoreBlankLines = true,
        Mode = CsvHelper.CsvMode.NoEscape,
        BadDataFound = null,
    };

    /// <summary>Full sequence, with modifications.</summary>
    [Name("Sequence")] public string Sequence { get; set; } = "";
    [Name("Base Sequence")] public string BaseSequence { get; set; } = "";

    /// <summary>IsoTracker output only.</summary>
    [Name("Peak Order")] [Optional] public int? PeakOrder { get; set; }

    /// <summary>Protein groups the peptide maps to, <c>;</c>-joined.</summary>
    [Name("Protein Groups")] [Optional] public string? ProteinGroups { get; set; }
    [Name("Gene Names")] [Optional] public string? GeneNames { get; set; }
    [Name("Organism")] [Optional] public string? Organism { get; set; }

    /// <summary>Per-sample values keyed by the label after the column prefix, verbatim.</summary>
    [Ignore]
    public IReadOnlyDictionary<string, QuantifiedPeptideSample> Samples { get; internal set; }
        = new Dictionary<string, QuantifiedPeptideSample>();
}

/// <summary>One sample's columns in a peptide row.</summary>
public sealed class QuantifiedPeptideSample
{
    internal QuantifiedPeptideSample(string label) => Label = label;

    /// <summary>The label, verbatim from the header.</summary>
    public string Label { get; }

    /// <summary><c>Intensity_</c>. FlashLFQ writes 0 when the peptide was not quantified; null only for a blank cell.</summary>
    public double? Intensity { get; internal set; }

    /// <summary><c>Detection Type_</c>: MSMS, MBR, NotDetected, and FlashLFQ's other values, verbatim.</summary>
    public string? DetectionType { get; internal set; }

    /// <summary><c>RetentionTime (min)_</c>: IsoTracker output only; null otherwise.</summary>
    public double? RetentionTime { get; internal set; }
}
