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
        public double PosteriorErrorProbability;

        public Identification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups,
            ChemicalFormula optionalChemicalFormula = null, bool useForProteinQuant = true, double posteriorErrorProbability = 0)
        {
            this.FileInfo = fileInfo;
            this.BaseSequence = BaseSequence;
            this.ModifiedSequence = ModifiedSequence;
            this.MonoisotopicMass = monoisotopicMass;
            this.Ms2RetentionTimeInMinutes = ms2RetentionTimeInMinutes;
            this.PrecursorChargeState = chargeState;
            this.ProteinGroups = new HashSet<ProteinGroup>(proteinGroups);
            this.OptionalChemicalFormula = optionalChemicalFormula;
            UseForProteinQuant = useForProteinQuant;
            PosteriorErrorProbability = posteriorErrorProbability;
        }

        public override string ToString()
        {
            return ModifiedSequence;
        }
    }
}