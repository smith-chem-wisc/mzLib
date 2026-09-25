using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using Omics.BioPolymer;
using Omics.BioPolymerGroup;
using System.Globalization;

namespace Readers;

/// <summary>
/// One row of a MetaMorpheus protein-group table (<c>AllQuantifiedProteinGroups.tsv</c>), read back.
/// The writer is <see cref="BioPolymerGroupTsvSchema"/>. Both its column vocabulary
/// ("BioPolymer Accession", "Unique Sequences") and the one MetaMorpheus 1.1.x wrote
/// ("Protein Accession", "Unique Peptides") are accepted.
/// </summary>
/// <remarks>
/// <para><b>The file is written unfiltered.</b> Decoys, contaminants and groups above 1% FDR are all
/// rows. Filter on <see cref="QValue"/> and <see cref="DecoyContaminantTarget"/> before counting.</para>
/// <para><b>No member is a leading protein.</b> MetaMorpheus sorts a group's members by accession
/// (<c>BioPolymerGroup.ListOfBioPolymersOrderedByAccession</c>), so the first accession is simply the
/// one that sorts first.</para>
/// </remarks>
public class ProteinGroupFromTsv
{
    /// <summary>Tab-separated, invariant culture, and no quote handling: MetaMorpheus does not quote,
    /// and a quote inside a modification name must stay a character.</summary>
    public static CsvConfiguration CsvConfiguration => new(CultureInfo.InvariantCulture)
    {
        Delimiter = "\t",
        HasHeaderRecord = true,
        IgnoreBlankLines = true,
        Mode = CsvHelper.CsvMode.NoEscape,
        BadDataFound = null,
    };

    /// <summary>The group's accessions, <c>|</c>-joined, sorted by accession.</summary>
    [Name("BioPolymer Accession", "Protein Accession")]
    public string ProteinGroupName { get; set; } = "";

    [Name("Gene")] [Optional] public string? Gene { get; set; }
    [Name("Organism")] [Optional] public string? Organism { get; set; }
    [Name("BioPolymer Full Name", "Protein Full Name")] [Optional] public string? FullName { get; set; }

    /// <summary><c>|</c>-joined, one per member. Kept as text: a member with no computable mass is "NaN".</summary>
    [Name("BioPolymer Unmodified Mass", "Protein Unmodified Mass")] [Optional] public string? UnmodifiedMass { get; set; }

    [Name("Number of BioPolymers in Group", "Number of Proteins in Group")] [Optional] public int? NumberOfMembers { get; set; }
    [Name("Unique Sequences", "Unique Peptides")] [Optional] public string? UniqueSequences { get; set; }
    [Name("Shared Sequences", "Shared Peptides")] [Optional] public string? SharedSequences { get; set; }
    [Name("Number of Sequences", "Number of Peptides")] [Optional] public int? NumberOfSequences { get; set; }
    [Name("Number of Unique Sequences", "Number of Unique Peptides")] [Optional] public int? NumberOfUniqueSequences { get; set; }

    /// <summary><c>|</c>-joined, one per member.</summary>
    [Name("Sequence Coverage Fraction")] [Optional] public string? SequenceCoverageFraction { get; set; }
    [Name("Sequence Coverage")] [Optional] public string? SequenceCoverage { get; set; }
    [Name("Sequence Coverage with Mods")] [Optional] public string? SequenceCoverageWithMods { get; set; }
    [Name("Fragment Sequence Coverage")] [Optional] public string? FragmentSequenceCoverage { get; set; }

    /// <summary>PSMs at 1% FDR assigned to the group.</summary>
    [Name("Number of PSMs")] [Optional] public int? NumberOfPsms { get; set; }

    /// <summary><c>T</c>, <c>D</c>, <c>C</c>, or entrapment <c>ET</c>/<c>ED</c>.</summary>
    [Name("BioPolymer Decoy/Contaminant/Target", "Protein Decoy/Contaminant/Target")]
    public string DecoyContaminantTarget { get; set; } = "";

    [Name("BioPolymer Cumulative Target", "Protein Cumulative Target")] [Optional] public int? CumulativeTarget { get; set; }
    [Name("BioPolymer Cumulative Decoy", "Protein Cumulative Decoy")] [Optional] public int? CumulativeDecoy { get; set; }

    /// <summary>The group's q-value. The table is unfiltered; this is what filtering reads.</summary>
    [Name("BioPolymer QValue", "Protein QValue")]
    public double QValue { get; set; }

    [Name("Best Sequence Score", "Best Peptide Score")] [Optional] public double? BestScore { get; set; }
    [Name("Best Sequence Notch QValue", "Best Peptide Notch QValue")] [Optional] public double? BestNotchQValue { get; set; }

    /// <summary>Written by MetaMorpheus 1.1.x only.</summary>
    [Name("Best Peptide PEP")] [Optional] public double? BestPep { get; set; }

    private string? _accessionsSource;
    private string[] _accessions = [];

    /// <summary>The members, in the order written, which is accession order. Split once per
    /// <see cref="ProteinGroupName"/> value, not on every read.</summary>
    [Ignore]
    public IReadOnlyList<string> Accessions
    {
        get
        {
            if (!ReferenceEquals(_accessionsSource, ProteinGroupName))
            {
                _accessions = ProteinGroupName.Split('|');
                _accessionsSource = ProteinGroupName;
            }
            return _accessions;
        }
    }

    [Ignore] public bool IsDecoy => DecoyContaminantTargetLabel.IsDecoy(DecoyContaminantTarget);
    [Ignore] public bool IsContaminant => DecoyContaminantTarget == DecoyContaminantTargetLabel.Contaminant;
    [Ignore] public bool IsEntrapment => DecoyContaminantTargetLabel.IsEntrapment(DecoyContaminantTarget);

    /// <summary>
    /// The per-sample-group measurements, keyed by the label that follows the column prefix, verbatim
    /// (for example <c>QE-002106_GM1_a-calib</c>: a calibrated file keeps its <c>-calib</c> stem).
    /// Labels are opaque: they may contain underscores, and condition, replicate or channel cannot be
    /// recovered from them. An isobaric file has one counting label and one intensity label per channel,
    /// so an entry may carry only the counting half or only the intensity half.
    /// </summary>
    [Ignore]
    public IReadOnlyDictionary<string, SampleGroupMeasurement> SampleGroups { get; internal set; }
        = new Dictionary<string, SampleGroupMeasurement>();
}

/// <summary>
/// One sample group's columns in a protein-group row. Every value is null when its column is absent
/// from the file or its cell is blank; a blank is "not reported", which is not zero.
/// </summary>
public sealed class SampleGroupMeasurement
{
    internal SampleGroupMeasurement(string label) => Label = label;

    /// <summary>The label, verbatim from the header.</summary>
    public string Label { get; }

    /// <summary><c>SpectralCount_</c>: PSMs from this file assigned to the group.</summary>
    public int? SpectralCount { get; internal set; }

    /// <summary><c>Intensity_</c>: blank when the group was not quantified in this sample group.</summary>
    public double? Intensity { get; internal set; }

    private string? _countOccupancyText;
    private string? _intensityOccupancyText;
    private ModificationOccupancyCell? _countOccupancy;
    private ModificationOccupancyCell? _intensityOccupancy;

    /// <summary>The <c>CountOccupancy_</c> cell, verbatim.</summary>
    public string? CountOccupancyText
    {
        get => _countOccupancyText;
        internal set { _countOccupancyText = value; _countOccupancy = null; }
    }

    /// <summary>The <c>IntensityOccupancy_</c> cell, verbatim.</summary>
    public string? IntensityOccupancyText
    {
        get => _intensityOccupancyText;
        internal set { _intensityOccupancyText = value; _intensityOccupancy = null; }
    }

    /// <summary>The count cell, parsed on first read and kept. Trust its (modified/total) pair over its
    /// rounded fraction. Throws <see cref="FormatException"/>, on every read, if the cell is malformed.</summary>
    public ModificationOccupancyCell CountOccupancy =>
        _countOccupancy ??= ModificationOccupancyCell.Parse(CountOccupancyText);

    /// <summary>The intensity cell, parsed on first read and kept. Trust its fraction; the pair is rounded
    /// to four significant digits. Throws <see cref="FormatException"/>, on every read, if the cell is malformed.</summary>
    public ModificationOccupancyCell IntensityOccupancy =>
        _intensityOccupancy ??= ModificationOccupancyCell.Parse(IntensityOccupancyText);
}
