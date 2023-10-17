using System.Globalization;
using System.Text;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MassSpectrometry;

namespace Readers.ExternalResults.IndividualResultRecords;

public class ToppicPrsm
{
    [Ignore]
    public static CsvConfiguration CsvConfiguration => new CsvConfiguration(CultureInfo.InvariantCulture)
    {
        Encoding = Encoding.UTF8,
        HasHeaderRecord = true,
        Delimiter = "\t",
    };

    public ToppicPrsm()
    {
        AlternativeIdentifications =
            new List<(int ScanNum, string Accession, string ProteinDescription, int FirstResidue, int LastResidue
                )>();
    }

    private string _fileNameWithoutExtension;
    [Ignore]
    public string FileNameWithoutExtension => _fileNameWithoutExtension ??= Path.GetFileNameWithoutExtension(FilePath);

    [Name("Data file name")]
    public string FilePath { get; set; }

    [Name("Prsm ID")]
    public int PrsmID { get; set; }

    [Name("Spectrum ID")]
    public int ScanNum { get; set; }

    [Name("Fragmentation")]
    public DissociationType DissociationType { get; set; }

    [Name("Scan(s)")]
    public int Scans { get; set; }

    [Name("Retention time")]
    public double RetentionTime { get; set; }

    [Name("#peaks")]
    public int PeakCount { get; set; }

    [Name("Charge")]
    public int PrecursorCharge { get; set; }

    [Name("Precursor mass")]
    public double PrecursorMass { get; set; }

    [Name("Adjusted precursor mass")]
    public double AdjustedPrecursorMass { get; set; }

    [Name("Proteoform ID")]
    public double ProteoformId { get; set; }

    [Name("Feature intensity")]
    public double FeatureIntensity { get; set; }

    [Name("Feature score")]
    public double FeatureScore { get; set; }

    [Name("Feature apex time")]
    public double FeatureApexTime { get; set; }

    [Name("#Protein hits")]
    public int ProteinHitsCount { get; set; }

    [Name("Protein accession")]
    public string ProteinAccession { get; set; }

    [Name("Protein description")]
    public string ProteinDescription { get; set; }

    [Name("First residue")]
    public int FirstResidue { get; set; }

    [Name("Last residue")]
    public int LastResidue { get; set; }

    [Name("Special amino acids")]
    public string? SpecialAminoAcids { get; set; }

    [Name("Database protein sequence")]
    public string BaseSequence { get; set; }

    [Name("Proteoform")]
    public string FullSequence { get; set; }

    [Name("Proteoform mass")]
    public double FullSequenceMass { get; set; }

    [Name("Protein N-terminal form")]
    public string ProteinNTerminalForm { get; set; }

    [Name("Fixed PTMs")]
    public string? FixedPTMs { get; set; }

    [Name("#unexpected modifications")]
    public int UnexpectedModificationsCount { get; set; }

    /// <summary>
    /// The mass shift of the mod and its semi-localization
    /// -47:[10-14] means a mass shift of -47 Da, and the semi-localization is between the 10th and 14th amino acids
    /// </summary>
    [Name("unexpected modifications")]
    public string UnexpectedModifications { get; set; }

    [Name("#variable PTMs")]
    public int VariableModificationsCount { get; set; }

    [Name("variable PTMs")]
    public string VariableModifications { get; set; }

    [Name("MIScore")]
    [TypeConverter(typeof(DashToNullOrDoubleConverter))]
    public double? MIScore { get; set; }

    [Name("#matched peaks")]
    public int MatchedPeaksCount { get; set; }

    [Name("#matched fragment ions")]
    public int MatchedFragmentIonsCount { get; set; }

    [Name("E-value")]
    public double EValue { get; set; }

    [Name("Spectrum-level Q-value")]
    [TypeConverter(typeof(DashToNullOrDoubleConverter))]
    public double? QValueSpectrumLevel { get; set; }

    [Name("Proteoform-level Q-value")]
    [TypeConverter(typeof(DashToNullOrDoubleConverter))]
    public double? QValueProteoformLevel { get; set; }

    [Ignore]
    public List<(int PrsmId, string Accession, string ProteinDescription, int FirstResidue, int LastResidue)> AlternativeIdentifications { get; set; }
}