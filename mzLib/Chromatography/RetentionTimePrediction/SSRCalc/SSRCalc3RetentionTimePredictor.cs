using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction.SSRCalc
{
    /// <summary>
    /// SSRCalc3-based retention time predictor.
    /// SSRCalc3 only uses base sequence, so no modification checking needed.
    /// </summary>
    public class SSRCalc3RetentionTimePredictor : RetentionTimePredictorBase
    {
        private readonly SSRCalc3 _calculator;

        public override string PredictorName => "SSRCalc3";
        public override SeparationType SeparationType => SeparationType.HPLC;

        // SSRCalc3 doesn't use modifications, so no checking needed
        public override bool RequiresModificationChecking => false;

        public SSRCalc3RetentionTimePredictor(
            IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.UsePrimarySequence,
            SSRCalc3.Column column = SSRCalc3.Column.A300)
            : base(modHandlingMode)
        {
            _calculator = new SSRCalc3("SSRCalc3", column);
        }

        protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out string? failureReason)
        {
            if (peptide.BaseSequence.Contains('U'))
            {
                failureReason = "SSRCalc3 does not support selenocysteine (U)";
                return false;
            }

            if (peptide.BaseSequence.Length < 4)
            {
                failureReason = "SSRCalc3 requires peptides of at least 4 amino acids";
                return false;
            }

            failureReason = null;
            return true;
        }

        protected override double? PredictCore(IRetentionPredictable peptide)
        {
            return _calculator.ScoreSequence(peptide.BaseSequence);
        }

        protected override double? PredictWithBaseSequence(string baseSequence)
        {
            return _calculator.ScoreSequence(baseSequence);
        }
    }
}