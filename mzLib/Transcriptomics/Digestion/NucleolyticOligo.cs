using Chemistry;
using Omics.Digestion;
using Omics.Modifications;

namespace Transcriptomics
{
    /// <summary>
    /// The most basic form of a digested oligo, this class does not care about mass or formula, just base sequence
    /// </summary>
    public class NucleolyticOligo : LysisProduct
    {
        protected IHasChemicalFormula _fivePrimeTerminus;
        protected IHasChemicalFormula _threePrimeTerminus;

        internal NucleolyticOligo(NucleicAcid nucleicAcid, int oneBaseStartResidue,
            int oneBasedEndResidue, int missedCleavages, CleavageSpecificity cleavageSpecificity,
            IHasChemicalFormula? fivePrimeTerminus, IHasChemicalFormula? threePrimeTerminus) 
        :base(nucleicAcid, oneBaseStartResidue, oneBasedEndResidue, missedCleavages, cleavageSpecificity)
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

        internal IEnumerable<OligoWithSetMods> GetModifiedOligos(IEnumerable<Modification> allKnownFixedMods,
            RnaDigestionParams digestionParams, List<Modification> variableModifications)
        {
            int oligoLength = OneBasedEndResidue - OneBasedStartResidue + 1;
            int maximumVariableModificationIsoforms = digestionParams.MaxModificationIsoforms;
            int maxModsForOligo = digestionParams.MaxMods;
            var twoBasedPossibleVariableAndLocalizeableModifications = new Dictionary<int, List<Modification>>(oligoLength + 4);

            var fivePrimeVariableMods = new List<Modification>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(1, fivePrimeVariableMods);

            var threePrimeVariableMods = new List<Modification>();
            twoBasedPossibleVariableAndLocalizeableModifications.Add(oligoLength + 2, threePrimeVariableMods);

            foreach (Modification variableModification in variableModifications)
            {
                // Check if can be a n-term mod
                if (CanBeFivePrime(variableModification, oligoLength)/* && !ModificationLocalization.UniprotModExists(NucleicAcid, 1, variableModification)*/)
                {
                    fivePrimeVariableMods.Add(variableModification);
                }

                for (int r = 0; r < oligoLength; r++)
                {
                    if (variableModification.LocationRestriction == "Anywhere." &&
                        ModificationLocalization.ModFits(variableModification, NucleicAcid.BaseSequence, r + 1, oligoLength, OneBasedStartResidue + r)
                         /*&& !ModificationLocalization.UniprotModExists(NucleicAcid, r + 1, variableModification)*/)
                    {
                        if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<Modification> residueVariableMods))
                        {
                            residueVariableMods = new List<Modification> { variableModification };
                            twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                        }
                        else
                        {
                            residueVariableMods.Add(variableModification);
                        }
                    }
                }
                // Check if can be a c-term mod
                if (CanBeThreePrime(variableModification, oligoLength) /*&& !ModificationLocalization.UniprotModExists(NucleicAcid, oligoLength, variableModification)*/)
                {
                    threePrimeVariableMods.Add(variableModification);
                }
            }

            // LOCALIZED MODS
            foreach (var kvp in NucleicAcid.OneBasedPossibleLocalizedModifications)
            {
                bool inBounds = kvp.Key >= OneBasedStartResidue && kvp.Key <= OneBasedEndResidue;
                if (!inBounds)
                {
                    continue;
                }

                int locInPeptide = kvp.Key - OneBasedStartResidue + 1;
                foreach (Modification modWithMass in kvp.Value)
                {
                    if (modWithMass is Modification variableModification)
                    {
                        // Check if can be a n-term mod
                        if (locInPeptide == 1 && CanBeFivePrime(variableModification, oligoLength) && !NucleicAcid.IsDecoy)
                        {
                            fivePrimeVariableMods.Add(variableModification);
                        }

                        int r = locInPeptide - 1;
                        if (r >= 0 && r < oligoLength
                            && (NucleicAcid.IsDecoy ||
                            (ModificationLocalization.ModFits(variableModification, NucleicAcid.BaseSequence, r + 1, oligoLength, OneBasedStartResidue + r)
                             && variableModification.LocationRestriction == "Anywhere.")))
                        {
                            if (!twoBasedPossibleVariableAndLocalizeableModifications.TryGetValue(r + 2, out List<Modification> residueVariableMods))
                            {
                                residueVariableMods = new List<Modification> { variableModification };
                                twoBasedPossibleVariableAndLocalizeableModifications.Add(r + 2, residueVariableMods);
                            }
                            else
                            {
                                residueVariableMods.Add(variableModification);
                            }
                        }

                        // Check if can be a c-term mod
                        if (locInPeptide == oligoLength && CanBeThreePrime(variableModification, oligoLength) && !NucleicAcid.IsDecoy)
                        {
                            threePrimeVariableMods.Add(variableModification);
                        }
                    }
                }
            }

            int variable_modification_isoforms = 0;

            foreach (Dictionary<int, Modification> kvp in GetVariableModificationPatterns(twoBasedPossibleVariableAndLocalizeableModifications, maxModsForOligo, oligoLength))
            {
                int numFixedMods = 0;
                foreach (var ok in GetFixedModsOneIsNterminusOrFivePrime(oligoLength, allKnownFixedMods))
                {
                    if (!kvp.ContainsKey(ok.Key))
                    {
                        numFixedMods++;
                        kvp.Add(ok.Key, ok.Value);
                    }
                }
                yield return new OligoWithSetMods(NucleicAcid, digestionParams, OneBasedStartResidue, OneBasedEndResidue, MissedCleavages,
                    CleavageSpecificityForFdrCategory, kvp, numFixedMods, _fivePrimeTerminus, _threePrimeTerminus);
                variable_modification_isoforms++;
                if (variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    yield break;
                }
            }

            //yield return new OligoWithSetMods(NucleicAcid, digestionParams, OneBasedStartResidue,
            //    OneBasedEndResidue, MissedCleavages, CleavageSpecificityForFdrCategory,
            //    new Dictionary<int, Modification>(), 0, _fivePrimeTerminus,
            //    _threePrimeTerminus);
        }

        private bool CanBeFivePrime(Modification variableModification, int peptideLength)
        {
            return (variableModification.LocationRestriction == "5'-terminal." || variableModification.LocationRestriction == "Oligo 5'-terminal.")
                && ModificationLocalization.ModFits(variableModification, NucleicAcid.BaseSequence, 1, peptideLength, OneBasedStartResidue);
        }

        private bool CanBeThreePrime(Modification variableModification, int peptideLength)
        {
            return (variableModification.LocationRestriction == "3'-terminal." || variableModification.LocationRestriction == "Oligo 3'-terminal.")
                && ModificationLocalization.ModFits(variableModification, NucleicAcid.BaseSequence, peptideLength, peptideLength, OneBasedStartResidue + peptideLength - 1);
        }
    }
}
