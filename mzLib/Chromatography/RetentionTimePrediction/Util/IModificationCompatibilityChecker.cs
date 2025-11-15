namespace Chromatography.RetentionTimePrediction.Util
{
    /// <summary>
    /// Checks whether modifications in a peptide are compatible with a specific predictor.
    /// This interface allows the Omics layer to provide modification-checking logic
    /// without the Chromatography layer depending on the Modification class.
    /// </summary>
    public interface IModificationCompatibilityChecker
    {
        /// <summary>
        /// Checks if all modifications in the sequence are compatible with the predictor
        /// </summary>
        /// <param name="fullSequence">The full sequence with modification identifiers</param>
        /// <param name="incompatibleModIds">List of incompatible modification identifiers (if any)</param>
        /// <returns>True if all modifications are compatible</returns>
        bool AreModificationsCompatible(string fullSequence, out List<string>? incompatibleModIds);

        /// <summary>
        /// Builds a full sequence string with only compatible modifications.
        /// Used for RemoveIncompatibleMods mode.
        /// </summary>
        /// <param name="fullSequence">Original full sequence</param>
        /// <returns>Full sequence with incompatible modifications removed</returns>
        string FilterIncompatibleModifications(string fullSequence);

        /// <summary>
        /// Builds a mass shift sequence with only compatible modifications.
        /// Used for RemoveIncompatibleMods mode with mass-based predictors.
        /// </summary>
        /// <param name="originalMassShiftSequence">Original sequence with mass shifts</param>
        /// <param name="fullSequence">Full sequence with modification identifiers (for reference)</param>
        /// <returns>Mass shift sequence with incompatible modifications removed</returns>
        string FilterIncompatibleModificationsFromMassShifts(string originalMassShiftSequence, string fullSequence);
    }
}