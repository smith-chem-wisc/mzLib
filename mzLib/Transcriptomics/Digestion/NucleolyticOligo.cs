using Chemistry;
using Omics.Digestion;
using Omics.Modifications;

namespace Transcriptomics.Digestion
{
    /// <summary>
    /// The most basic form of a digested oligo, this class does not care about mass or formula, just base sequence
    /// </summary>
    public class NucleolyticOligo : DigestionProduct
    {
        protected IHasChemicalFormula _fivePrimeTerminus;
        protected IHasChemicalFormula _threePrimeTerminus;

        internal NucleolyticOligo(NucleicAcid nucleicAcid, int oneBaseStartResidue,
            int oneBasedEndResidue, int missedCleavages, CleavageSpecificity cleavageSpecificity,
            IHasChemicalFormula? fivePrimeTerminus, IHasChemicalFormula? threePrimeTerminus, string? description = null)
        : base(nucleicAcid, oneBaseStartResidue, oneBasedEndResidue, missedCleavages, cleavageSpecificity, description)
        {
            _fivePrimeTerminus = fivePrimeTerminus ?? NucleicAcid.DefaultFivePrimeTerminus;
            _threePrimeTerminus = threePrimeTerminus ?? NucleicAcid.DefaultThreePrimeTerminus;
        }

        /// <summary>
        /// Nucleic acid this oligo was digested from
        /// </summary>
        public NucleicAcid NucleicAcid
        {
            get => Parent as NucleicAcid;
            protected set => Parent = value;
        }

        public override string ToString()
        {
            return BaseSequence;
        }

        /// <summary>
        /// Generates a collection of oligos with set modifications based on the provided fixed and variable modifications,
        /// digestion parameters, and the nucleic acid sequence.
        /// </summary>
        /// <param name="allKnownFixedMods">A collection of all known fixed modifications.</param>
        /// <param name="digestionParams">Parameters for RNA digestion.</param>
        /// <param name="variableModifications">A list of variable modifications to consider.</param>
        /// <returns>An enumerable collection of oligos with set modifications.</returns>
        /// <remarks>
        /// Code heavily borrowed from ProteolyticPeptide.GetModifiedPeptides
        /// </remarks>
        internal IEnumerable<OligoWithSetMods> GenerateModifiedOligos(List<Modification> allKnownFixedMods,
            RnaDigestionParams digestionParams, List<Modification> variableModifications)
        {
            int variableModificationIsoforms = 0;
            int oligoLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForOligo = digestionParams.MaxMods;
            var twoBasedPossibleVariableAndLocalizeableModifications = DictionaryPool.Get();
            var fixedModDictionary = FixedModDictionaryPool.Get();

            try
            {
                PopulateVariableModifications(variableModifications, in twoBasedPossibleVariableAndLocalizeableModifications);
                PopulateFixedModsOneIsNorFivePrimeTerminus(oligoLength, allKnownFixedMods, in fixedModDictionary);

                // Add the mods to the oligo by return numerous OligoWithSetMods
                foreach (Dictionary<int, Modification> variableModPattern in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForOligo, oligoLength))
                {
                    AppendFixedModificationsToVariable(in fixedModDictionary, in variableModPattern, out int numFixedMods);

                    yield return new OligoWithSetMods(NucleicAcid, digestionParams, OneBasedStartResidue, OneBasedEndResidue, MissedCleavages,
                        CleavageSpecificityForFdrCategory, variableModPattern, numFixedMods, _fivePrimeTerminus, _threePrimeTerminus);

                    variableModificationIsoforms++;
                    if (variableModificationIsoforms == maximumVariableModificationIsoforms)
                    {
                        yield break;
                    }
                }
            }
            finally
            {
                DictionaryPool.Return(twoBasedPossibleVariableAndLocalizeableModifications);
                FixedModDictionaryPool.Return(fixedModDictionary);
            }
        }
    }
}
