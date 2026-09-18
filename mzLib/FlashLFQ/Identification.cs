using Chemistry;
using System.Collections.Generic;
using MassSpectrometry;

namespace FlashLFQ
{
    /// <summary>
    /// One identification handed to FlashLFQ as a starting point for quantification: a sequence observed
    /// in a particular file at a particular retention time, mass and charge. FlashLFQ does not identify
    /// anything itself - these come from a search engine - and it uses them as seeds for peak finding.
    /// </summary>
    public class Identification
    {
        public readonly string BaseSequence;
        public readonly string ModifiedSequence;
        public readonly double Ms2RetentionTimeInMinutes;
        public readonly double MonoisotopicMass;
        public readonly SpectraFileInfo FileInfo;
        public readonly int PrecursorChargeState;
        public readonly HashSet<ProteinGroup> ProteinGroups;
        public readonly ChemicalFormula OptionalChemicalFormula;
        public readonly bool UseForProteinQuant;
        public double PeakfindingMass;
        public double PsmScore { get; init; }
        public double QValue { get; init; }
        public bool IsDecoy { get; }

        /// <summary>
        /// Name of the digestion agent that produced this analyte, or null if unknown. Match-between-runs
        /// will not transfer an identification into a file digested with a different agent, because the
        /// peptide is not in that file's analyte pool. Null on either side means unrestricted.
        /// </summary>
        public string DigestionAgentName { get; init; }

        public Identification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups,
            ChemicalFormula optionalChemicalFormula = null, bool useForProteinQuant = true,
            double psmScore = 0, double qValue = 0, bool decoy = false, string digestionAgentName = null)
        {
            DigestionAgentName = digestionAgentName;
            this.FileInfo = fileInfo;
            this.BaseSequence = BaseSequence;
            this.ModifiedSequence = ModifiedSequence;
            this.MonoisotopicMass = monoisotopicMass;
            this.Ms2RetentionTimeInMinutes = ms2RetentionTimeInMinutes;
            this.PrecursorChargeState = chargeState;
            this.ProteinGroups = new HashSet<ProteinGroup>(proteinGroups);
            this.OptionalChemicalFormula = optionalChemicalFormula;
            UseForProteinQuant = !decoy && useForProteinQuant; // ensure that decoy peptides aren't used for protein quant
            QValue = qValue;
            PsmScore = psmScore;
            IsDecoy = decoy;
        }

        public override string ToString()
        {
            return ModifiedSequence;
        }
    }
}