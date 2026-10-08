#nullable enable
using System;
using System.Collections.Generic;
using Quantification.Differential;
using Readers;

namespace FlashLFQ
{
    /// <summary>
    /// Builds the differential-statistics input (<see cref="ObservationTable"/>) from FlashLFQ's peptide quantification,
    /// from either source (GR-18): the in-memory <see cref="FlashLfqResults"/>, or the stored peptide table FlashLFQ wrote
    /// from them (<c>AllQuantifiedPeptides.tsv</c>, read by <see cref="QuantifiedPeptideFile"/>). Both go through one
    /// mapping, so the same search gives the same table either way.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Detection types</b> (<c>DEF-PEP-DT</c>) map to states as follows. <c>MSMS</c> with an intensity above 0 is
    /// <see cref="ObservationState.Quantified"/>; <c>MBR</c> with an intensity above 0 is
    /// <see cref="ObservationState.MbrTransferred"/>. <c>MSMS</c> or <c>MBR</c> with intensity 0 is
    /// <see cref="ObservationState.AmbiguousPeak"/>: FlashLFQ zeroes every fraction of a sample whose strongest fraction is
    /// ambiguous, and leaves the type alone. <c>MSMSAmbiguousPeakfinding</c> is <see cref="ObservationState.AmbiguousPeak"/>
    /// whatever its intensity. <c>MSMSIdentifiedButNotQuantified</c> and <c>NotDetected</c> map to their namesakes.
    /// IsoTracker output is refused: its peptides are isobaric peaks, not sequences.
    /// </para>
    /// <para>
    /// <b>Protein groups</b> are taken from the producer's protein-group table or groups (GR-23), because FlashLFQ's own
    /// <see cref="ProteinGroup"/> carries no q-value and no decoy or contaminant flag.
    /// </para>
    /// </remarks>
    public static class FlashLfqObservations
    {
        /// <summary>The table from in-memory FlashLFQ results.</summary>
        /// <param name="results">FlashLFQ's results, after peptide quantification.</param>
        /// <param name="runs">The experimental design's runs; each matches one of <see cref="FlashLfqResults.SpectraFiles"/> by file name.</param>
        /// <param name="proteinGroups">The producer's protein groups (MetaMorpheus's, with q-values and flags).</param>
        /// <param name="basis">Which match-between-runs values enter.</param>
        /// <exception cref="ArgumentException">IsoTracker results, or input <see cref="ObservationTable.Build"/> refuses.</exception>
        public static ObservationTable FromResults(FlashLfqResults results, IEnumerable<ObservationRun> runs,
            IEnumerable<ProteinGroupInfo> proteinGroups, QuantBasis basis)
        {
            ArgumentNullException.ThrowIfNull(results);
            if (results.IsoTracker)
                throw new ArgumentException("IsoTracker results are not supported: their peptides are isobaric peaks.", nameof(results));

            var peptides = new List<PeptideRunValues>(results.PeptideModifiedSequences.Count);
            foreach (Peptide peptide in results.PeptideModifiedSequences.Values)
            {
                var byRun = new Dictionary<string, RunValue>(StringComparer.Ordinal);
                foreach (SpectraFileInfoLabel file in Labels(results))
                    byRun[file.Label] = ToRunValue(peptide.Sequence, file.Label, peptide.GetIntensity(file.File),
                        peptide.GetDetectionType(file.File));

                var groups = new List<string>();
                foreach (ProteinGroup group in peptide.ProteinGroups)
                    groups.Add(group.ProteinGroupName);
                peptides.Add(new PeptideRunValues(peptide.Sequence, peptide.BaseSequence, groups, byRun));
            }

            var labels = new List<string>();
            foreach (SpectraFileInfoLabel file in Labels(results))
                labels.Add(file.Label);
            return ObservationTable.Build(runs, labels, peptides, proteinGroups, basis);
        }

        /// <summary>The table from a stored peptide table (<c>AllQuantifiedPeptides.tsv</c>).</summary>
        /// <param name="peptideTable">The peptide table.</param>
        /// <param name="runs">The experimental design's runs; each matches one of the table's columns by file name.</param>
        /// <param name="proteinGroups">The producer's protein groups, e.g. from <see cref="ReadProteinGroups"/>.</param>
        /// <param name="basis">Which match-between-runs values enter.</param>
        /// <exception cref="ArgumentException">
        /// IsoTracker output, a blank intensity or an unreadable detection type (naming the peptide and column), or input
        /// <see cref="ObservationTable.Build"/> refuses.
        /// </exception>
        public static ObservationTable FromPeptideTable(QuantifiedPeptideFile peptideTable, IEnumerable<ObservationRun> runs,
            IEnumerable<ProteinGroupInfo> proteinGroups, QuantBasis basis)
        {
            ArgumentNullException.ThrowIfNull(peptideTable);
            ArgumentNullException.ThrowIfNull(runs);
            var runList = new List<ObservationRun>(runs);

            List<string>? labels = null;
            var peptides = new List<PeptideRunValues>(peptideTable.Results.Count);
            foreach (QuantifiedPeptideFromTsv row in peptideTable.Results)
            {
                if (row.PeakOrder != null || HasRetentionTimes(row))
                    throw new ArgumentException(
                        $"'{peptideTable.FilePath}' is IsoTracker output, which is not supported: its peptides are isobaric peaks.",
                        nameof(peptideTable));

                labels ??= new List<string>(row.Samples.Keys);

                var byRun = new Dictionary<string, RunValue>(StringComparer.Ordinal);
                foreach (var (label, sample) in row.Samples)
                {
                    if (sample.Intensity is not double intensity)
                        throw new ArgumentException(
                            $"Peptide '{row.Sequence}' has a blank Intensity_{label}; FlashLFQ writes 0 for no value.",
                            nameof(peptideTable));
                    byRun[label] = ToRunValue(row.Sequence, label, intensity, ParseDetectionType(row.Sequence, label, sample.DetectionType));
                }

                string[] groups = (row.ProteinGroups ?? "").Split(';', StringSplitOptions.RemoveEmptyEntries);
                peptides.Add(new PeptideRunValues(row.Sequence, row.BaseSequence, groups, byRun));
            }

            // An empty table names no columns; its runs are then the design's own.
            labels ??= runList.ConvertAll(r => r.FileName);
            return ObservationTable.Build(runList, labels, peptides, proteinGroups, basis);
        }

        /// <summary>Every row of a MetaMorpheus protein-group table (<c>AllQuantifiedProteinGroups.tsv</c>) as a <see cref="ProteinGroupInfo"/>.</summary>
        /// <remarks>
        /// The decoy, contaminant and entrapment flags are the reader's reading of the <c>Protein Decoy/Contaminant/Target</c>
        /// label (<see cref="ProteinGroupFromTsv.IsDecoy"/> and its siblings).
        /// </remarks>
        public static IReadOnlyList<ProteinGroupInfo> ReadProteinGroups(ProteinGroupFromTsvFile proteinGroupTable)
        {
            ArgumentNullException.ThrowIfNull(proteinGroupTable);
            return proteinGroupTable.Results.ConvertAll(row => new ProteinGroupInfo(row.ProteinGroupName, row.QValue,
                row.IsDecoy, row.IsContaminant, row.IsEntrapment, row.Gene));
        }

        /// <summary>
        /// One run's value from FlashLFQ's intensity and detection type (the mapping in the class remarks). Both sources
        /// call this, so they cannot map a type differently.
        /// </summary>
        private static RunValue ToRunValue(string sequence, string label, double intensity, DetectionType detectionType)
        {
            switch (detectionType)
            {
                case DetectionType.MSMS:
                    return intensity > 0
                        ? new RunValue(intensity, ObservationState.Quantified)
                        : new RunValue(0, ObservationState.AmbiguousPeak);
                case DetectionType.MBR:
                    return intensity > 0
                        ? new RunValue(intensity, ObservationState.MbrTransferred)
                        : new RunValue(0, ObservationState.AmbiguousPeak);
                case DetectionType.MSMSAmbiguousPeakfinding:
                    return new RunValue(0, ObservationState.AmbiguousPeak);
                case DetectionType.MSMSIdentifiedButNotQuantified:
                    return new RunValue(0, ObservationState.IdentifiedNotQuantified);
                case DetectionType.NotDetected:
                    return new RunValue(0, ObservationState.NotDetected);
                default:
                    throw new ArgumentException(
                        $"Peptide '{sequence}' in '{label}' has detection type {detectionType}, which is IsoTracker's and not supported.");
            }
        }

        private static DetectionType ParseDetectionType(string sequence, string label, string? text)
        {
            if (!string.IsNullOrEmpty(text) && char.IsLetter(text[0])
                && Enum.TryParse(text, ignoreCase: false, out DetectionType parsed) && Enum.IsDefined(parsed))
                return parsed;

            throw new ArgumentException(
                $"Peptide '{sequence}' has Detection Type_{label} '{text}', which is not a FlashLFQ detection type.");
        }

        private static bool HasRetentionTimes(QuantifiedPeptideFromTsv row)
        {
            foreach (QuantifiedPeptideSample sample in row.Samples.Values)
                if (sample.RetentionTime != null)
                    return true;
            return false;
        }

        private readonly record struct SpectraFileInfoLabel(MassSpectrometry.SpectraFileInfo File, string Label);

        /// <summary>Each file with the label FlashLFQ writes for it: its name without extension.</summary>
        private static IEnumerable<SpectraFileInfoLabel> Labels(FlashLfqResults results)
        {
            foreach (MassSpectrometry.SpectraFileInfo file in results.SpectraFiles)
                yield return new SpectraFileInfoLabel(file, file.FilenameWithoutExtension);
        }
    }
}
