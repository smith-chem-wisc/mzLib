using Chemistry;
using System.Collections.Generic;

namespace FlashLFQ
{
    public class Identification
    {
        public readonly string BaseSequence;
        public readonly string ModifiedSequence;
        public readonly double ms2RetentionTimeInMinutes;
        public readonly double monoisotopicMass;
        public readonly SpectraFileInfo fileInfo;
        public readonly int precursorChargeState;
        public readonly HashSet<ProteinGroup> proteinGroups;
        public readonly ChemicalFormula OptionalChemicalFormula;
        public readonly bool UseForProteinQuant;
        public double massToLookFor;

        public Identification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups,
            ChemicalFormula optionalChemicalFormula = null, bool useForProteinQuant = true)
        {
            this.fileInfo = fileInfo;
            this.BaseSequence = BaseSequence;
            this.ModifiedSequence = ModifiedSequence;
            this.monoisotopicMass = monoisotopicMass;
            this.ms2RetentionTimeInMinutes = ms2RetentionTimeInMinutes;
            this.precursorChargeState = chargeState;
            this.proteinGroups = new HashSet<ProteinGroup>(proteinGroups);
            this.OptionalChemicalFormula = optionalChemicalFormula;
            UseForProteinQuant = useForProteinQuant;
        }

        public override string ToString()
        {
            return ModifiedSequence;
        }
    }
}