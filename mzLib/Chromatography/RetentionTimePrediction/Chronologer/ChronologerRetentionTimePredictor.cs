using Chromatography.RetentionTimePrediction.Util;
using TorchSharp;

namespace Chromatography.RetentionTimePrediction.Chronologer
{
    /// <summary>
    /// Chronologer-based retention time predictor using deep learning.
    /// Predicts C18 retention times reported in % ACN.
    /// </summary>
    public class ChronologerRetentionTimePredictor : RetentionTimePredictorBase, IDisposable
    {
        private readonly Chronologer _model;

        public override string PredictorName => "Chronologer";
        public override bool RequiresModificationChecking => true;

        /// <summary>
        /// Initializes a new Chronologer predictor with default weights.
        /// </summary>
        public ChronologerRetentionTimePredictor(
            IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.RemoveIncompatibleMods)
            : this(modHandlingMode, weightsPath: null)
        {
        }

        /// <summary>
        /// Initializes a new Chronologer predictor with custom weights file.
        /// </summary>
        public ChronologerRetentionTimePredictor(
            IncompatibleModHandlingMode modHandlingMode,
            string? weightsPath)
            : base(modHandlingMode)
        {
            _model = weightsPath != null
                ? new Chronologer(weightsPath)
                : new Chronologer();
        }

        protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out string? failureReason)
        {
            return ChronologerSequenceEncoder.CanEncode(
                peptide.BaseSequence,
                peptide.GetSequenceWithMassShifts(),
                out failureReason);
        }

        protected override double? PredictCore(IRetentionPredictable peptide)
        {
            string? massShiftSequence = peptide.GetSequenceWithMassShifts();
            if (massShiftSequence == null)
            {
                // Fallback to base sequence
                return PredictFromSequence(peptide.BaseSequence);
            }

            return PredictFromSequence(massShiftSequence);
        }

        protected override double? PredictWithBaseSequence(string baseSequence)
        {
            return PredictFromSequence(baseSequence);
        }

        protected override double? PredictWithFilteredSequences(
            string baseSequence,
            string filteredFullSequence,
            string? filteredMassShiftSequence)
        {
            if (filteredMassShiftSequence != null)
            {
                return PredictFromSequence(filteredMassShiftSequence);
            }

            return PredictWithBaseSequence(baseSequence);
        }

        /// <summary>
        /// Predicts retention time from a sequence string (with or without mass shifts).
        /// </summary>
        private double? PredictFromSequence(string sequence)
        {
            var tensor = ChronologerSequenceEncoder.EncodeTensor(sequence);
            if (tensor is null)
                return null;

            try
            {
                var prediction = _model.Predict(tensor);
                return prediction[0].ToDouble();
            }
            finally
            {
                tensor?.Dispose();
            }
        }

        public void Dispose()
        {
            _model.Dispose();
        }
    }
}