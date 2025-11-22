namespace Chromatography.RetentionTimePrediction
{
    /// <summary>
    /// Represents an entity for which retention time can be predicted.
    /// Designed for minimal allocation and no dependencies on higher layers (e.g., Omics).
    /// </summary>
    public interface IRetentionPredictable
    {
        /// <summary>
        /// Gets the base (unmodified) sequence
        /// </summary>
        string BaseSequence { get; }

        /// <summary>
        /// Gets the full sequence representation with modification identifiers
        /// e.g., "PEPTIDE[Variable:Oxidation on M]K[Variable:Acetylation on K]"
        /// </summary>
        string FullSequence { get; }

        /// <summary>
        /// Gets the monoisotopic mass of the peptide.
        /// Required for CZE electrophoretic mobility predictions.
        /// </summary>
        double MonoisotopicMass { get; }

        /// <summary>
        /// Builds a sequence string with mass shifts for modifications.
        /// Format: "PEPTIDE[+15.995]K[+42.011]" or "PEPTIDE[-17.026]K"
        /// This is used by predictors that work with mass-based representations (e.g., Chronologer).
        /// </summary>
        /// <returns>Sequence with mass shift annotations, or null if not applicable</returns>
        string FullSequenceWithMassShifts { get; }
    }
}