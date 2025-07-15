using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using MassSpectrometry;
using static System.Net.Mime.MediaTypeNames;

namespace Readers;

/// <summary>
/// Class Representing a TopPIC prsm or proteoform
/// For supported versions and software this file type can come from see
///     Readers.ExternalResources.SupportedVersions.txt
/// </summary>
/// <remarks>
/// Things that could be done to improve compatibility:
/// Convert Variable Modifications to a list of Modification objects
/// Convert NTerminalForm to a Modification object
/// </remarks>
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
        AlternativeIdentifications = new List<AlternativeToppicId>();
    }

    private string? _fileNameWithoutExtension;
    [Ignore]
    public string FileNameWithoutExtension => _fileNameWithoutExtension ??= Path.GetFileNameWithoutExtension(FilePath);

    [Name("Data file name")]
    public string FilePath { get; set; }

    [Name("Prsm ID")]
    public int PrsmID { get; set; }

    [Name("Spectrum ID")]
    public int SpectrumId { get; set; }

    [Name("Fragmentation")]
    public DissociationType DissociationType { get; set; }

    [Name("Scan(s)")]
    public int OneBasedScanNumber { get; set; }

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
    [Format("#.00#E+00")]
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

    [Ignore]
    private string? _baseSequence;

    [Optional]
    [Name("Database protein sequence")]
    public string BaseSequence
    {
        get => _baseSequence ??= GetBaseSequenceFromFullSequence();
        set => _baseSequence = value;
    }

    [Name("Proteoform")]
    public string FullSequence { get; set; }

    [Name("Proteoform mass")]
    public double FullSequenceMass { get; set; }

    [Name("Protein N-terminal form")]
    public string ProteinNTerminalForm { get; set; }

    [Optional]
    [Name("Fixed PTMs")]
    public string? FixedPTMs { get; set; }
    
    [Name("#unexpected modifications")]
    public int UnexpectedModificationsCount { get; set; }

    /// <summary>
    /// The mass shift of the mod and its semi-localization
    /// -47:[10-14] means a mass shift of -47 Da, and the semi-localization is between the 10th and 14th amino acids
    /// </summary>
    [Optional]
    [Name("unexpected modifications")]
    public string UnexpectedModifications { get; set; }
    
    [Name("#variable PTMs")]
    public int VariableModificationsCount { get; set; }

    [Optional]
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
    [Format("0.00E+00")]
    public double EValue { get; set; }

    [Name("Spectrum-level Q-value")]
    [TypeConverter(typeof(DashToNullOrDoubleConverter))]
    public double? QValueSpectrumLevel { get; set; }

    [Name("Proteoform-level Q-value")]
    [TypeConverter(typeof(DashToNullOrDoubleConverter))]
    public double? QValueProteoformLevel { get; set; }

    [Ignore]
    public List<AlternativeToppicId> AlternativeIdentifications { get; set; }

    public string GetBaseSequenceFromFullSequence()
    {
        // Remove text within square brackets
        var text = Regex.Replace(FullSequence, @"\[[^\]]*\]", "");

        // Remove parentheses
        text = Regex.Replace(text, @"[()]", "");

        // Remove periods
        text = Regex.Replace(text, @"(^[^.]+)|(\.[^.]+$)", "")
            .Replace(".","");
        return text;
    }
}

/// <summary>
/// Class representing an alternative Identification from the tsv file. 
/// </summary>
public class AlternativeToppicId
{
    public int PrsmId { get; set; }
    public string Accession { get; set; }
    public string ProteinDescription { get; set; }
    public int FirstResidue { get; set; }
    public int LastResidue { get; set; }

    public AlternativeToppicId(int prsmId, string accession, string proteinDescription, int firstResidue,
        int lastResidue)
    {
        PrsmId = prsmId;
        Accession = accession;
        ProteinDescription = proteinDescription;
        FirstResidue = firstResidue;
        LastResidue = lastResidue;
    }

    public override string ToString()
    {
        return $"{PrsmId}\t\t\t\t\t\t\t\t\t\t\t\t\t\t{Accession}\t{ProteinDescription}\t{FirstResidue}\t{LastResidue}\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
    }
}