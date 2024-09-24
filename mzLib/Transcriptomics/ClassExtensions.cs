using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics.Digestion;

namespace Transcriptomics
{
    public static class ClassExtensions
    {
        /// <summary>
        /// Creates a new instance of a nucleic acid or oligo with set modifications, optionally updating its sequence, modifications, and decoy status.
        /// </summary>
        /// <typeparam name="T">The type of the nucleic acid, which must implement <see cref="INucleicAcid"/>.</typeparam>
        /// <param name="target">The target nucleic acid or oligo with set modifications to base the new instance on.</param>
        /// <param name="sequence">The new sequence string, if any. If null, the original sequence is used.</param>
        /// <param name="modifications">A dictionary of modifications to apply, if any. If null, the original modifications are used.</param>
        /// <param name="isDecoy">A flag indicating whether the sequence is a decoy, if any. If null, the original decoy status is used.</param>
        /// <returns>A new instance of the specified nucleic acid type with the provided or existing properties.</returns>
        /// <remarks>
        /// This method facilitates the generation of new sequences for both nucleic acids and oligos with set modifications by allowing
        /// optional updates to the sequence string, modifications, and decoy status. It ensures that the new instances are properly
        /// initialized with the provided or existing properties, enabling further analysis of modified sequences and future generation of decoys on the fly.
        /// </remarks>
        public static T CreateNew<T>(this T target, string? sequence = null, IDictionary<int, List<Modification>>? modifications = null,
        bool? isDecoy = null)
            where T : INucleicAcid
        {
            // set new object parameters where not null
            object? returnObj = null;
            string newSequence = sequence ?? target.BaseSequence;
            IDictionary<int, List<Modification>> newModifications = modifications ?? target.OneBasedPossibleLocalizedModifications;

            switch (target)
            {
                case RNA rna:
                {
                    bool newIsDecoy = isDecoy ?? rna.IsDecoy;
                    returnObj = new RNA(newSequence, rna.Name, rna.Accession, rna.Organism, rna.DatabaseFilePath,
                        rna.FivePrimeTerminus, rna.ThreePrimeTerminus, newModifications, rna.IsContaminant, newIsDecoy, rna.AdditionalDatabaseFields);
                    break;
                }
                case OligoWithSetMods oligo:
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
                        (oligo.DigestionParams as RnaDigestionParams)!,
                        oligo.OneBasedStartResidue,
                        oligo.OneBasedEndResidue,
                        oligo.MissedCleavages,
                        oligo.CleavageSpecificityForFdrCategory,
                        newModifications.ToDictionary(p => p.Key, p => p.Value.First()),
                        oligo.NumFixedMods,
                        oligo.FivePrimeTerminus,
                        oligo.ThreePrimeTerminus);
                    break;
                }
                default:
                    throw new ArgumentException("INucleicAcid type not yet implemented");
            }

            return (T)returnObj ?? throw new NullReferenceException("Error creating new INucleicAcid");
        }
    }
}
