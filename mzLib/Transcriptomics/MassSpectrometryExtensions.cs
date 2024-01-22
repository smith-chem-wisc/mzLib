using Omics.Modifications;
using Transcriptomics.Digestion;

namespace Transcriptomics
{
    public static class MassSpectrometryExtensions
    {
        public static T CreateNew<T>(this T target, string? sequence = null, IDictionary<int, List<Modification>>? modifications = null,
        bool? isDecoy = null) 
            where T : INucleicAcid
        {
            // set new object parameters where not null
            object? returnObj = null;
            string newSequence = sequence ?? target.BaseSequence;
            IDictionary<int, List<Modification>> newModifications = modifications ?? target.OneBasedPossibleLocalizedModifications;
            

            if (target is RNA rna)
            {
                bool newIsDecoy = isDecoy ?? rna.IsDecoy;
                returnObj = new RNA(newSequence, rna.Name, rna.Accession, rna.Organism, rna.DatabaseFilePath,
                    rna.FivePrimeTerminus, rna.ThreePrimeTerminus, newModifications, rna.IsContaminant, newIsDecoy, rna.AdditionalDatabaseFields);
            }
            else if (target is OligoWithSetMods oligo)
            {
                var oldParent = oligo.Parent as RNA ?? throw new NullReferenceException();
                var newParent = new RNA(
                    newSequence,
                    oldParent.Name, 
                    oldParent.Accession,
                    oldParent.Organism,
                    oldParent.DatabaseFilePath,
                    oldParent.FivePrimeTerminus, 
                    oldParent.ThreePrimeTerminus, 
                    newModifications,
                    oldParent.IsContaminant, 
                    oldParent.IsDecoy,
                    oldParent.AdditionalDatabaseFields);

                returnObj = new OligoWithSetMods(
                    newParent,
                    oligo.DigestionParams as RnaDigestionParams,
                    oligo.OneBasedStartResidue, 
                    oligo.OneBasedEndResidue, 
                    oligo.MissedCleavages,
                    oligo.CleavageSpecificityForFdrCategory,
                    newModifications.ToDictionary(p => p.Key, p => p.Value.First()),
                    oligo.NumFixedMods,
                    oligo.FivePrimeTerminus,
                    oligo.ThreePrimeTerminus);
            }
            else 
                throw new ArgumentException("Target must be RNA or OligoWithSetMods");

            return (T)returnObj ?? throw new NullReferenceException();
        }
    }
}
