using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction
{
    /// <summary>
    /// Base class providing common functionality for RT predictors.
    /// Handles modification checking and strategy application without depending on Modification class.
    /// </summary>
    public abstract class RetentionTimePredictorBase : IRetentionTimePredictor
    {
        protected IncompatibleModHandlingMode ModHandlingMode { get; }
        protected IModificationCompatibilityChecker? ModChecker { get; private set; }

        protected RetentionTimePredictorBase(IncompatibleModHandlingMode modHandlingMode)
        {
            ModHandlingMode = modHandlingMode;
        }
        /// <summary>
        /// Gets the separation type this predictor is designed for
        /// </summary>
        public virtual SeparationType SeparationType => SeparationType.HPLC;

        public abstract string PredictorName { get; }
        public abstract bool RequiresModificationChecking { get; }

        public void SetModificationChecker(IModificationCompatibilityChecker checker)
        {
            ModChecker = checker;
        }

        public double? PredictRetentionTime(IRetentionPredictable peptide)
        {
            // Quick validation of basic constraints
            if (!ValidateBasicConstraints(peptide, out string? basicFailure))
            {
                if (ModHandlingMode == IncompatibleModHandlingMode.ThrowException)
                    throw new InvalidOperationException(
                        $"Cannot predict RT for peptide '{peptide.FullSequence}' with {PredictorName}: {basicFailure}");
                return null;
            }

            // If this predictor doesn't require modification checking, use fast path
            if (!RequiresModificationChecking)
            {
                return PredictCore(peptide);
            }

            // Check modifications if checker is available
            if (ModChecker == null)
            {
                // No checker available - decide based on strategy
                if (ModHandlingMode == IncompatibleModHandlingMode.ThrowException)
                {
                    throw new InvalidOperationException(
                        $"Predictor '{PredictorName}' requires modification checking but no checker was provided.");
                }
                // For other modes, proceed without checking (assume compatible or use base sequence)
                return HandleNoModificationChecker(peptide);
            }

            // Check for incompatible modifications
            bool areCompatible = ModChecker.AreModificationsCompatible(
                peptide.FullSequence,
                out List<string>? incompatibleModIds);

            if (areCompatible)
            {
                // Fast path - all modifications compatible
                return PredictCore(peptide);
            }

            // Handle incompatible modifications
            return HandleIncompatibleModifications(peptide, incompatibleModIds!);
        }

        public bool CanPredict(IRetentionPredictable peptide, out string? failureReason)
        {
            if (!ValidateBasicConstraints(peptide, out failureReason))
                return false;

            if (RequiresModificationChecking && ModChecker != null)
            {
                bool compatible = ModChecker.AreModificationsCompatible(
                    peptide.FullSequence,
                    out List<string>? incompatibleModIds);

                if (!compatible)
                {
                    if (ModHandlingMode == IncompatibleModHandlingMode.ThrowException ||
                        ModHandlingMode == IncompatibleModHandlingMode.ReturnNull)
                    {
                        failureReason = $"{incompatibleModIds!.Count} incompatible modification(s)";
                        return false;
                    }
                }
            }

            failureReason = null;
            return true;
        }

        /// <summary>
        /// Validate basic constraints (length, amino acid content, etc.)
        /// </summary>
        protected abstract bool ValidateBasicConstraints(IRetentionPredictable peptide, out string? failureReason);

        /// <summary>
        /// Core prediction logic - called when peptide passes all validation
        /// </summary>
        protected abstract double? PredictCore(IRetentionPredictable peptide);

        /// <summary>
        /// Predict using only base sequence (no modifications)
        /// </summary>
        protected abstract double? PredictWithBaseSequence(string baseSequence);

        /// <summary>
        /// Predict with incompatible modifications filtered out.
        /// Default implementation falls back to base sequence.
        /// Override if predictor can use filtered modifications.
        /// </summary>
        protected virtual double? PredictWithFilteredModifications(IRetentionPredictable peptide)
        {
            // Get filtered sequences from the checker
            string filteredFullSequence = ModChecker!.FilterIncompatibleModifications(peptide.FullSequence);
            string? originalMassShifts = peptide.GetSequenceWithMassShifts();
            
            if (originalMassShifts != null)
            {
                string filteredMassShifts = ModChecker!.FilterIncompatibleModificationsFromMassShifts(
                    originalMassShifts,
                    peptide.FullSequence);
                
                // Create temporary wrapper for prediction
                return PredictWithFilteredSequences(peptide.BaseSequence, filteredFullSequence, filteredMassShifts);
            }

            return PredictWithFilteredSequences(peptide.BaseSequence, filteredFullSequence, null);
        }

        /// <summary>
        /// Predict using filtered sequences. Override if predictor needs custom handling.
        /// </summary>
        protected virtual double? PredictWithFilteredSequences(
            string baseSequence,
            string filteredFullSequence,
            string? filteredMassShiftSequence)
        {
            // Default: fall back to base sequence
            return PredictWithBaseSequence(baseSequence);
        }

        private double? HandleNoModificationChecker(IRetentionPredictable peptide)
        {
            // If no checker available, apply strategy
            switch (ModHandlingMode)
            {
                case IncompatibleModHandlingMode.UsePrimarySequence:
                case IncompatibleModHandlingMode.RemoveIncompatibleMods:
                    // Without checker, we can't filter - use base sequence
                    return PredictWithBaseSequence(peptide.BaseSequence);

                case IncompatibleModHandlingMode.ReturnNull:
                    return null;

                default:
                    return PredictCore(peptide);
            }
        }

        private double? HandleIncompatibleModifications(IRetentionPredictable peptide, List<string> incompatibleModIds)
        {
            switch (ModHandlingMode)
            {
                case IncompatibleModHandlingMode.ReturnNull:
                    return null;

                case IncompatibleModHandlingMode.ThrowException:
                    throw new IncompatibleModificationException(
                        peptide.FullSequence,
                        incompatibleModIds,
                        PredictorName);

                case IncompatibleModHandlingMode.UsePrimarySequence:
                    return PredictWithBaseSequence(peptide.BaseSequence);

                case IncompatibleModHandlingMode.RemoveIncompatibleMods:
                    return PredictWithFilteredModifications(peptide);

                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

    }
}