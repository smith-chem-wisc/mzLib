using System.Text;
using Chemistry;
using Chromatography.RetentionTimePrediction;
using Omics.Modifications;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers.DeNovoMapping;

/// <summary>
/// One Casanovo prediction prepared for mapping: the sequence used for comparison, the normalized score
/// and the observed and predicted retention times.
/// </summary>
public class DeNovoPrediction : IRetentionPredictable
{
    public DeNovoPrediction(CasanovoMzTabRecord record, string fileKey, bool replaceIWithL)
    {
        Record = record;
        FileKey = fileKey;

        // Casanovo writes N-terminal modifications as "[mod]-PEPTIDE"; the dash is not a residue.
        BaseSequence = record.BaseSequence.Replace("-", string.Empty);
        MappingSequence = NormalizeSequence(BaseSequence, replaceIWithL);

        // Casanovo reports a score in [-1, 0) when the precursor mass is outside its tolerance.
        // Shift those into [0, 1) so every score lies between 0 and 1.
        Score = record.Score < 0 ? record.Score + 1 : record.Score;
        ObservedRetentionTime = record.RetentionTime;
    }

    public CasanovoMzTabRecord Record { get; }

    /// <summary>
    /// Name of the mzTab file the prediction came from, without extension. Used for per-file counts.
    /// </summary>
    public string FileKey { get; }

    /// <summary>
    /// The predicted peptide with modifications and dashes removed and, when requested, I replaced by L.
    /// </summary>
    public string MappingSequence { get; }

    /// <summary>
    /// Casanovo score shifted into [0, 1].
    /// </summary>
    public double Score { get; }

    /// <summary>
    /// Retention time in minutes, from the mzTab or, when missing there, from the spectra file.
    /// </summary>
    public double? ObservedRetentionTime { get; set; }

    /// <summary>
    /// Chronologer prediction, in Chronologer's own units (% acetonitrile equivalent).
    /// </summary>
    public double? PredictedRetentionTime { get; set; }

    public string BaseSequence { get; }
    public string FullSequence => Record.FullSequence;
    public double MonoisotopicMass => Record.CalculatedMassToCharge.ToMass(Record.Charge);
    public string FullSequenceWithMassShifts => BuildFullSequenceWithMassShifts();

    public static string NormalizeSequence(string sequence, bool replaceIWithL) =>
        replaceIWithL ? sequence.Replace('I', 'L') : sequence;

    private string BuildFullSequenceWithMassShifts()
    {
        var mods = Record.AllModsOneIsNterminus;
        var builder = new StringBuilder();

        if (mods.TryGetValue(1, out Modification? mod))
            AppendMass(builder, mod);

        for (int r = 0; r < BaseSequence.Length; r++)
        {
            builder.Append(BaseSequence[r]);
            if (mods.TryGetValue(r + 2, out mod))
                AppendMass(builder, mod);
        }

        if (mods.TryGetValue(BaseSequence.Length + 2, out mod))
        {
            builder.Append('-');
            AppendMass(builder, mod);
        }

        return builder.ToString();
    }

    private static void AppendMass(StringBuilder builder, Modification mod)
    {
        double mass = mod.MonoisotopicMass ?? 0;
        builder.Append(mass > 0 ? $"[+{mass.RoundedDouble(6)}]" : $"[{mass.RoundedDouble(6)}]");
    }
}
