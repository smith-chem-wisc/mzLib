using Chemistry;
using MassSpectrometry;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.ProteolyticDigestion
{
    [Serializable]
    public class PeptideWithSetModifications : ProteolyticPeptide
    {
        public string FullSequence { get; private set; } //sequence with modifications
        public readonly int NumFixedMods;

        /// <summary>
        /// Dictionary of modifications on the peptide. The N terminus is index 1.
        /// The key indicates which residue modification is on (with 1 being N terminus).
        /// </summary>
        [NonSerialized] private Dictionary<int, Modification> _allModsOneIsNterminus; //we currently only allow one mod per position

        [NonSerialized] private bool? _hasChemicalFormulas;
        [NonSerialized] private string _sequenceWithChemicalFormulas;
        [NonSerialized] private double? _monoisotopicMass;
        [NonSerialized] private Dictionary<FragmentationTerminus, CompactPeptide> _compactPeptides;
        [NonSerialized] private DigestionParams _digestionParams;
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private readonly string DigestionParamString; // used to get digestion param object after deserialization
        private readonly string ProteinAccession; // used to get protein object after deserialization

        /// <summary>
        /// Creates a PeptideWithSetModifications object from a protein. Used when a Protein is digested.
        /// </summary>
        public PeptideWithSetModifications(Protein protein, DigestionParams digestionParams, int oneBasedStartResidueInProtein,
            int oneBasedEndResidueInProtein, string peptideDescription, int missedCleavages,
           Dictionary<int, Modification> allModsOneIsNterminus, int numFixedMods)
           : base(protein, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, peptideDescription)
        {
            _allModsOneIsNterminus = allModsOneIsNterminus;
            NumFixedMods = numFixedMods;
            _digestionParams = digestionParams;
            DetermineFullSequence();
            ProteinAccession = protein.Accession;
            DigestionParamString = digestionParams.ToString();
        }

        /// <summary>
        /// Creates a PeptideWithSetModifications object from a sequence string.
        /// Useful for reading in MetaMorpheus search engine output into mzLib objects.
        /// </summary>
        public PeptideWithSetModifications(string sequence, Dictionary<string, Modification> allKnownMods, int numFixedMods = 0,
            DigestionParams digestionParams = null, Protein p = null, int oneBasedStartResidueInProtein = int.MinValue,
            int oneBasedEndResidueInProtein = int.MinValue, int missedCleavages = int.MinValue, string peptideDescription = null)
            : base(p, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, peptideDescription)
        {
            if (sequence.Contains("|"))
            {
                throw new MzLibUtil.MzLibException("Ambiguous peptide cannot be parsed from string: " + sequence);
            }

            FullSequence = sequence;
            GetModsAfterDeserialization(allKnownMods, out _baseSequence);
            NumFixedMods = numFixedMods;
            _digestionParams = digestionParams;

            if (p != null)
            {
                ProteinAccession = p.Accession;
            }
            if (digestionParams != null)
            {
                DigestionParamString = digestionParams.ToString();
            }
        }

        public DigestionParams DigestionParams
        {
            get { return _digestionParams; }
        }

        public Dictionary<int, Modification> AllModsOneIsNterminus
        {
            get { return _allModsOneIsNterminus; }
        }

        public int NumMods
        {
            get { return AllModsOneIsNterminus.Count; }
        }

        public int NumVariableMods
        {
            get { return NumMods - NumFixedMods; }
        }

        public double MonoisotopicMass
        {
            get
            {
                if (!_monoisotopicMass.HasValue)
                {
                    double monoMass = WaterMonoisotopicMass;

                    foreach (var mod in AllModsOneIsNterminus.Values)
                    {
                        monoMass += mod.MonoisotopicMass.Value;
                    }
                    monoMass += BaseSequence.Select(b => Residue.ResidueMonoisotopicMass[b]).Sum();

                    _monoisotopicMass = monoMass;
                }
                return (double)ClassExtensions.RoundedDouble(_monoisotopicMass.Value);
            }
        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (!_hasChemicalFormulas.HasValue)
                {
                    _hasChemicalFormulas = true;
                    var subsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                    {
                        if (pep_n_term_variable_mod is Modification jj)
                        {
                            subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                        }
                        else
                        {
                            return null;
                        }
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        subsequence.Append(this[r]);
                        // variable modification on this residue
                        if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                        {
                            if (residue_variable_mod is Modification jj)
                            {
                                subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                            }
                            else
                            {
                                return null;
                            }
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                    {
                        if (pep_c_term_variable_mod is Modification jj)
                        {
                            subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                        }
                        else
                        {
                            return null;
                        }
                    }

                    _sequenceWithChemicalFormulas = subsequence.ToString();
                }
                return _sequenceWithChemicalFormulas;
            }
        }

        /// <summary>
        /// Generates theoretical fragments for given dissociation type for this peptide
        /// </summary>
        public IEnumerable<Product> Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus)
        {
            // molecular ion
            //yield return new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.None, this.MonoisotopicMass, Length, Length), 0);

            var productCollection = TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[fragmentationTerminus].Intersect(DissociationTypeCollection.ProductsFromDissociationType[dissociationType]);

            List<(ProductType, int)> skippers = new List<(ProductType, int)>();
            foreach (var product in productCollection.Where(f => f != ProductType.zPlusOne))
            {
                skippers.Add((product, BaseSequence.Length));
            }

            switch (dissociationType)
            {
                case DissociationType.CID:
                    skippers.Add((ProductType.b, 1));
                    break;

                case DissociationType.ETD:
                case DissociationType.ECD:
                case DissociationType.EThcD:
                    skippers.AddRange(GetProlineZIonIndicies());
                    break;
            }

            foreach (var productType in productCollection)
            {
                // we're separating the N and C terminal masses and computing a separate compact peptide for each one
                // this speeds calculations up without producing unnecessary terminus fragment info
                FragmentationTerminus temporaryFragmentationTerminus = TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType];
                NeutralTerminusFragment[] terminalMasses = CompactPeptide(temporaryFragmentationTerminus).TerminalMasses;

                for (int f = 0; f < terminalMasses.Length; f++)
                {
                    // fragments with neutral loss
                    if (AllModsOneIsNterminus.TryGetValue(terminalMasses[f].AminoAcidPosition + 1, out Modification mod) && mod.NeutralLosses != null
                        && mod.NeutralLosses.TryGetValue(dissociationType, out List<double> neutralLosses))
                    {
                        foreach (double neutralLoss in neutralLosses)
                        {
                            if (neutralLoss == 0)
                            {
                                continue;
                            }

                            for (int n = f; n < terminalMasses.Length; n++)
                            {
                                if (!skippers.Contains((productType, terminalMasses[n].FragmentNumber)))
                                {
                                    yield return new Product(productType, terminalMasses[n], neutralLoss);
                                }
                            }
                        }
                    }

                    // "normal" fragment without neutral loss
                    if (!skippers.Contains((productType, terminalMasses[f].FragmentNumber)))
                    {
                        yield return new Product(productType, terminalMasses[f], 0);
                    }
                }
            }

            if (AllModsOneIsNterminus != null)
            {
                foreach (Modification mod in AllModsOneIsNterminus.Values)
                {
                    // molecular ion minus neutral losses
                    if (mod.NeutralLosses != null && mod.NeutralLosses.TryGetValue(dissociationType, out List<double> losses))
                    {
                        foreach (double neutralLoss in losses)
                        {
                            if (neutralLoss != 0)
                            {
                                yield return new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.Both, MonoisotopicMass, 0, 0), neutralLoss);
                            }
                        }
                    }

                    // diagnostic ions
                    if (mod.DiagnosticIons != null && mod.DiagnosticIons.TryGetValue(dissociationType, out List<double> diagnosticIons))
                    {
                        foreach (double diagnosticIon in diagnosticIons)
                        {
                            // the diagnostic ion is assumed to be annotated in the mod info as the *neutral mass* of the diagnostic ion, not the ionized species
                            yield return new Product(ProductType.D, new NeutralTerminusFragment(FragmentationTerminus.Both, diagnosticIon, 0, 0), 0);
                        }
                    }
                }
            }
        }

        public virtual string EssentialSequence(IReadOnlyDictionary<string, int> modstoWritePruned)
        {
            string essentialSequence = BaseSequence;
            if (modstoWritePruned != null)
            {
                var sbsequence = new StringBuilder();

                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_n_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_n_term_variable_mod.ModificationType + ":" + pep_n_term_variable_mod.IdWithMotif + ']');
                    }
                }
                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // variable modification on this residue
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                    {
                        if (modstoWritePruned.ContainsKey(residue_variable_mod.ModificationType))
                        {
                            sbsequence.Append('[' + residue_variable_mod.ModificationType + ":" + residue_variable_mod.IdWithMotif + ']');
                        }
                    }
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_c_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_c_term_variable_mod.ModificationType + ":" + pep_c_term_variable_mod.IdWithMotif + ']');
                    }
                }

                essentialSequence = sbsequence.ToString();
            }
            return essentialSequence;
        }

        private IEnumerable<(ProductType, int)> GetProlineZIonIndicies()
        {
            for (int i = BaseSequence.IndexOf('P'); i > -1; i = BaseSequence.IndexOf('P', i + 1))
            {
                yield return (ProductType.zPlusOne, BaseSequence.Length - i);
            }
        }

        public CompactPeptide CompactPeptide(FragmentationTerminus fragmentationTerminus)
        {
            // need this for deserialization
            if (_compactPeptides == null)
            {
                _compactPeptides = new Dictionary<FragmentationTerminus, CompactPeptide>();
            }

            if (_compactPeptides.TryGetValue(fragmentationTerminus, out CompactPeptide compactPeptide))
            {
                return compactPeptide;
            }

            CompactPeptide cp = new CompactPeptide(this, fragmentationTerminus);

            lock (_compactPeptides)
            {
                if (!_compactPeptides.ContainsKey(fragmentationTerminus))
                {
                    _compactPeptides.Add(fragmentationTerminus, cp);
                }
            }

            return cp;
        }

        public PeptideWithSetModifications Localize(int j, double massToLocalize)
        {
            var dictWithLocalizedMass = new Dictionary<int, Modification>(AllModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (dictWithLocalizedMass.TryGetValue(j + 2, out Modification modToReplace))
            {
                massOfExistingMod = (double)modToReplace.MonoisotopicMass;
                dictWithLocalizedMass.Remove(j + 2);
            }

            dictWithLocalizedMass.Add(j + 2, new Modification(_locationRestriction: "Anywhere.", _monoisotopicMass: massToLocalize + massOfExistingMod));

            var peptideWithLocalizedMass = new PeptideWithSetModifications(Protein, _digestionParams, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                PeptideDescription, MissedCleavages, dictWithLocalizedMass, NumFixedMods);

            return peptideWithLocalizedMass;
        }

        public override string ToString()
        {
            return FullSequence + string.Join("\t", AllModsOneIsNterminus.Select(m => m.ToString()));
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;

            if (Protein == null && q.Protein == null)
            {
                return q.FullSequence.Equals(this.FullSequence);
            }

            return q != null
                && q.FullSequence.Equals(this.FullSequence)
                && q.OneBasedStartResidueInProtein == this.OneBasedStartResidueInProtein
                && (q.Protein.Accession == null && this.Protein.Accession == null || q.Protein.Accession.Equals(this.Protein.Accession))
                && q.DigestionParams.Protease == this.DigestionParams.Protease;
        }

        public override int GetHashCode()
        {
            return FullSequence.GetHashCode() + DigestionParams.Protease.GetHashCode();
        }

        /// <summary>
        /// This should be run after deserialization of a PeptideWithSetModifications, in order to set its Protein and Modification objects, which were not serialized
        /// </summary>
        public void SetNonSerializedPeptideInfo(Dictionary<string, Modification> idToMod, Dictionary<string, Protein> accessionToProtein)
        {
            GetModsAfterDeserialization(idToMod, out string baseSequence);
            GetProteinAfterDeserialization(accessionToProtein);
            GetDigestionParamsAfterDeserialization();
        }

        private void GetDigestionParamsAfterDeserialization()
        {
            if (DigestionParamString != null)
            {
                _digestionParams = DigestionParams.FromString(DigestionParamString);
            }
        }

        private void GetModsAfterDeserialization(Dictionary<string, Modification> idToMod, out string baseSequence)
        {
            _allModsOneIsNterminus = new Dictionary<int, Modification>();
            StringBuilder baseSequenceSb = new StringBuilder();
            StringBuilder currentModification = new StringBuilder();
            int currentModificationLocation = 1;
            bool currentlyReadingMod = false;
            int bracketCount = 0;

            for (int r = 0; r < FullSequence.Length; r++)
            {
                char c = FullSequence[r];

                switch (c)
                {
                    case '[':
                        currentlyReadingMod = true;

                        if (bracketCount > 0)
                        {
                            currentModification.Append(c);
                        }

                        bracketCount++;

                        break;

                    case ']':
                        string modId = null;
                        bracketCount--;

                        if (bracketCount == 0)
                        {
                            try
                            {
                                string modString = currentModification.ToString();
                                int splitIndex = modString.IndexOf(':');
                                string modType = modString.Substring(0, splitIndex);
                                modId = modString.Substring(splitIndex + 1, modString.Length - splitIndex - 1);
                            }
                            catch (Exception e)
                            {
                                throw new MzLibUtil.MzLibException(
                                    "Error while trying to parse string into peptide: " + e.Message);
                            }

                            if (!idToMod.TryGetValue(modId, out Modification mod))
                            {
                                throw new MzLibUtil.MzLibException(
                                    "Could not find modification while reading string: " + FullSequence);
                            }

                            AllModsOneIsNterminus.Add(currentModificationLocation, mod);
                            currentlyReadingMod = false;
                            currentModification = new StringBuilder();
                        }
                        else
                        {
                            currentModification.Append(c);
                        }

                        break;

                    default:
                        if (currentlyReadingMod)
                        {
                            currentModification.Append(c);
                        }
                        else
                        {
                            currentModificationLocation++;
                            baseSequenceSb.Append(c);
                        }
                        break;
                }
            }

            baseSequence = baseSequenceSb.ToString();
        }

        private void GetProteinAfterDeserialization(Dictionary<string, Protein> idToProtein)
        {
            Protein protein = null;

            if (ProteinAccession != null && !idToProtein.TryGetValue(ProteinAccession, out protein))
            {
                throw new MzLibUtil.MzLibException("Could not find protein accession after deserialization! " + ProteinAccession);
            }

            Protein = protein;
        }

        private void DetermineFullSequence()
        {
            var subsequence = new StringBuilder();

            // modification on peptide N-terminus
            if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
            {
                subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subsequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            FullSequence = subsequence.ToString();
        }
    }
}