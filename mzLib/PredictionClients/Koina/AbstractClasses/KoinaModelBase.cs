using Chromatography.RetentionTimePrediction;
using System.Text.RegularExpressions;

namespace PredictionClients.Koina.AbstractClasses
{
    public abstract class KoinaModelBase: IDisposable
    {
        #region Model Metadata
        /// <summary>
        /// Gets the model name as registered in the Koina API.
        /// </summary>
        public abstract string ModelName { get; }
        
        /// <summary>
        /// Gets the maximum number of sequences allowed per API request batch.
        /// Can dig in the Koina github repo to find these values if needed.
        /// </summary>
        public abstract int MaxBatchSize { get; }

        public abstract int MaxNumberOfBatchesPerRequest { get; }
        #endregion

        #region Input Sequence Validation Constraints
        /// <summary>
        /// Gets the maximum allowed peptide base sequence length.
        /// </summary>
        public abstract int MaxPeptideLength { get; }
        
        /// <summary>
        /// Gets the minimum allowed peptide base sequence length.
        /// </summary>
        public abstract int MinPeptideLength { get; }

        public virtual IncompatibleModHandlingMode ModHandlingMode { get; } 

        /// <summary>
        /// Gets the mapping of valid modification annotations from mzLib format to UNIMOD format.
        /// Leaving this empty will only allow unmodified peptides. Models that allow modifications
        /// must override this property.
        /// </summary>
        public virtual Dictionary<string, string> ValidModificationUnimodMapping => new();
        
        /// <summary>
        /// Gets the regex pattern for validating amino acid sequences.
        /// </summary>
        public virtual string CanonicalAminoAcidPattern => @"^[ACDEFGHIKLMNPQRSTVWY]+$";
        
        /// <summary>
        /// Gets the regex pattern for identifying modifications.
        /// </summary>
        public virtual string ModificationPattern => @"\[[^\]]+\]";
        #endregion

        /// <summary>
        /// Peptides arte temporarily stored in the model instance to allow for batching and asynchronous prediction.
        /// </summary>
        protected virtual List<string> PeptideSequences { get; set; } = new();

        protected bool _disposed = false;

        #region Required Client Methods for Koina API Interaction
        /// <summary>
        /// Converts peptide sequences and associated data into batched request payloads for the Koina API.
        /// Implementations should group input sequences into batches respecting the MaxBatchSize constraint
        /// and format them according to the specific model's input requirements.
        /// </summary>
        /// <returns>List of request dictionaries, each containing a batch of sequences and parameters</returns>
        /// <remarks>
        /// Each dictionary in the returned list represents one API request batch and should contain:
        /// - Peptide sequences (formatted according to model requirements)
        /// - Model-specific parameters (e.g., charge states, collision energies, NCE values)
        /// - Any additional metadata required by the specific Koina model
        /// 
        /// The total number of sequences across all batches should equal PeptideSequences.Count.
        /// </remarks>
        protected abstract List<Dictionary<string, object>> ToBatchedRequests();
        #endregion

        #region Validation and Modification Handling
        /// <summary>
        /// Validates a peptide sequence against model constraints for modifications and basic sequence requirements.
        /// Handles incompatible modifications according to the specified ModHandlingMode.
        /// </summary>
        protected virtual bool ValidateBasicConstraints(string sequence, out string? failureReason)
        {
            switch (ModHandlingMode)
            {
                case IncompatibleModHandlingMode.RemoveIncompatibleMods:
                    var allMods = Regex.Matches(sequence, ModificationPattern).Select(m => m.Value).ToList();
                    foreach (var mod in allMods)
                    {
                        if (!ValidModificationUnimodMapping.ContainsKey(mod))
                        {
                            sequence = sequence.Replace(mod, ""); // Remove incompatible modification
                        }
                    }
                    break;

                case IncompatibleModHandlingMode.UsePrimarySequence:
                    sequence = Regex.Replace(sequence, ModificationPattern, ""); // Use primary sequence only
                    break;

                case IncompatibleModHandlingMode.ThrowException:
                    if (Regex.IsMatch(sequence, ModificationPattern))
                    {
                        throw new InvalidOperationException("Sequence contains unsupported modifications.");
                    }
                    break;

                case IncompatibleModHandlingMode.ReturnNull:
                    if (!HasValidModifications(sequence))
                    {
                        failureReason = "Sequence contains unsupported modifications.";
                        return false;
                    }
                    break;
            }
            failureReason = null;
            return true;
        }

        /// <summary>
        /// Validates that all modifications in a peptide sequence are supported by the model.
        /// Returns true for sequences without modifications or when all modifications are recognized.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications in mzLib format</param>
        /// <returns>True if all modifications are valid or no modifications present; otherwise false</returns>
        /// <remarks>
        /// Modification validation process:
        /// 1. Uses ModificationPattern regex to find all modification annotations
        /// 2. Checks each modification against ValidModificationUnimodMapping
        /// 3. Empty modification list (unmodified peptides) is considered valid
        /// </remarks>
        protected virtual bool HasValidModifications(string sequence)
        {
            var matches = Regex.Matches(sequence, ModificationPattern);
            if (matches.Count == 0)
            {
                return true; // No modifications found - valid for all models
            }

            // Check if any modifications are not in the valid mapping
            return matches.All(m => ValidModificationUnimodMapping.ContainsKey(m.Value));
        }

        /// <summary>
        /// Validates that a peptide sequence meets model constraints for length and amino acid composition.
        /// Removes modifications before validation to check only the base amino acid sequence.
        /// </summary>
        /// <param name="sequence">Peptide sequence with potential modifications</param>
        /// <returns>True if sequence meets model length and composition requirements; otherwise false</returns>
        /// <remarks>
        /// Validation criteria:
        /// - Strips modification annotations using ModificationPattern
        /// - Checks against CanonicalAminoAcidPattern for valid amino acids
        /// - Validates length is within [MinPeptideLength, MaxPeptideLength] range
        /// </remarks>
        protected virtual bool IsValidSequence(string sequence)
        {
            // Remove modification annotations to get base sequence
            var baseSequence = Regex.Replace(sequence, ModificationPattern, "");

            return Regex.IsMatch(baseSequence, CanonicalAminoAcidPattern) // Valid amino acids only
                && baseSequence.Length <= MaxPeptideLength                 // Within max length
                && baseSequence.Length >= MinPeptideLength;                // Above min length (implicit from abstract property)
        }
        #endregion

        #region Full Sequence Modification Conversion Methods
        /// <summary>
        /// Converts peptide sequence from mzLib modification format to UNIMOD format required by the model.
        /// Base implementation performs standard mzLib to UNIMOD conversion using the ValidModificationUnimodMapping.
        /// </summary>
        /// <param name="sequence">Peptide sequence in mzLib modification format</param>
        /// <returns>Sequence converted to UNIMOD format</returns>
        /// <remarks>
        /// Conversion process:
        /// 1. Replaces mzLib modification names with UNIMOD identifiers using ValidModificationUnimodMapping
        /// 
        /// Example transformations:
        /// - "PEPT[Common Variable:Oxidation on M]IDE" -> "PEPT[UNIMOD:35]IDE" (oxidation converted)
        /// - "C[Common Fixed:Carbamidomethyl on C]PEPTIDE" -> "C[UNIMOD:4]PEPTIDE" (carbamidomethyl converted)
        /// 
        /// Derived classes may override this method to implement model-specific modification handling,
        /// such as automatic carbamidomethylation of cysteines or other required modifications.
        /// </remarks>
        protected virtual string ConvertMzLibModificationsToUnimod(string sequence)
        {
            // Apply custom modification mappings first
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Key, mod.Value);
            }
            return sequence;
        }

        /// <summary>
        /// Converts peptide sequence from UNIMOD format back to mzLib modification format.
        /// Inverse operation of ConvertMzLibModificationsToUnimod.
        /// </summary>
        protected virtual string ConvertUnimodToMzLibModifications(string sequence)
        {
            foreach (var mod in ValidModificationUnimodMapping)
            {
                sequence = sequence.Replace(mod.Value, mod.Key);
            }
            return sequence;
        }
        #endregion
        public void Dispose()
        {
            if (!_disposed)
            {
                // Dispose of any resources if necessary.
                _disposed = true;
            }
        }
    }
}