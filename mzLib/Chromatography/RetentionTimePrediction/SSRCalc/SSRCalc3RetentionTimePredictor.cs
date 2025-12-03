using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction.SSRCalc;

/// <summary>
/// SSRCalc3-based retention time predictor.
/// SSRCalc3 only uses base sequence, so no modification checking needed.
/// </summary>
public class SSRCalc3RetentionTimePredictor : RetentionTimePredictor
{
    private readonly SSRCalc3 _calculator;
    public override string PredictorName => "SSRCalc3";
    public override SeparationType SeparationType => SeparationType.HPLC;

    public SSRCalc3RetentionTimePredictor(SSRCalc3.Column column = SSRCalc3.Column.A300)
    {
        _calculator = new SSRCalc3("SSRCalc3", column);
    }

    public override string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;
        return peptide.BaseSequence;
    }

    protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        var baseSequence = peptide.BaseSequence;
        if (baseSequence.Any(aa => Array.IndexOf(CanonicalAminoAcids, aa) == -1))
        {
            failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
            return false;
        }

        return base.ValidateBasicConstraints(peptide, out failureReason);
    }

    protected override double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null)
    {
        return _calculator.ScoreSequence(peptide.BaseSequence);
    }
}