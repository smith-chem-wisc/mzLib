using CsvHelper.Configuration.Attributes;
using Omics;
using Omics.Modifications;

namespace Readers.ExternalResults.IndividualResultRecords;

public class CasanovoMzTabRecord : ISpectralMatch
{
    [Name("PSM_ID")]
    public int Id { get; set; }

    [Name("sequence")]
    public string ModifiedSequence { get; set; }

    [Name("accession")]
    [NullValues("null")]
    public string? Accession { get; set; }

    [Name("unique")]
    [NullValues("null")]
    public string? Unique { get; set; }
    
    [Name("database")]
    [NullValues("null")]
    public string? Database { get; set; }

    [Name("database_version")]
    [NullValues("null")]
    public string? DatabaseVersion { get; set; }

    [Name("search_engine")]
    public string SearchEngine { get; set; }

    [Name("search_engine_score[1]")]
    public double Score { get; set; }

    [Name("modifications")]
    [NullValues("null")]
    public string? Modifications { get; set; }

    [Name("retention_time")]
    [NullValues("null")]
    public double? RetentionTime { get; set; }

    [Name("charge")]
    public int Charge { get; set; }

    [Name("exp_mass_to_charge")]
    public double ExperimentalMassToCharge { get; set; }

    [Name("calc_mass_to_charge")]
    public double CalculatedMassToCharge { get; set; }

    [Name("spectra_ref")]
    public string SpectraRef { get; set; }

    [Name("pre")]
    [NullValues("null")]
    public string? Pre { get; set; }

    [Name("post")]
    [NullValues("null")]
    public string? Post { get; set; }

    [Name("start")]
    [NullValues("null")]
    public int? Start { get; set; }

    [Name("end")]
    [NullValues("null")]
    public int? End { get; set; }

    [Name("opt_ms_run[1]_aa_scores")]
    [TypeConverter(typeof(CommaDelimitedToDoubleArrayTypeConverter))]
    public required double[] AminoAcidScores { get; set; }


    #region Interpreted Fields - Set during reading

    /// <summary>
    /// One Based Scan Number. CAREFUL: With MGF files, casanovo uses the index of the spectrum in the file, not the actual scan number.
    /// </summary>
    [Ignore] public int OneBasedScanNumber { get; set; }
    [Ignore] public bool IsDecoy { get; set; } = false;
    [Ignore] public string BaseSequence { get; internal set; } = string.Empty;
    [Ignore] public Dictionary<int, Modification> AllModsOneIsNterminus { get; internal set; } = [];
    [Ignore] public string FullSequence { get; internal set; } = string.Empty;
    [Ignore] public string FileNameWithoutExtension { get; internal set; } = string.Empty;

    #endregion
}
