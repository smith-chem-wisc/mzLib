using Chemistry;
using System.Collections.Generic;

namespace FlashLFQ
{
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

        public Identification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups,
            ChemicalFormula optionalChemicalFormula = null, bool useForProteinQuant = true,
            double psmScore = 0, double qValue = 0, bool decoy = false)
        {
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