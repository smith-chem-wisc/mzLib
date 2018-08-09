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
        [NonSerialized] public readonly DigestionParams DigestionParams;
        public string Sequence { get; private set; }
        public readonly int NumFixedMods;

        /// <summary>
        /// Dictionary of modifications on the peptide. The N terminus is index 1.
        /// The key indicates which residue modification is on (with 1 being N terminus).
        /// </summary>
        [NonSerialized] private Dictionary<int, Modification> _allModsOneIsNterminus; //we currently only allow one mod per position
        [NonSerialized] private bool? HasChemicalFormulas;
        [NonSerialized] private string _sequenceWithChemicalFormulas;
        [NonSerialized] private double? _monoisotopicMass;
        [NonSerialized] private Dictionary<FragmentationTerminus, CompactPeptide> _compactPeptides;
        private readonly string ProteinAccession; // used to get protein object after deserialization
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        
        
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
            DigestionParams = digestionParams;
            DetermineFullSequence();
            this.ProteinAccession = protein.Accession;
        }

        /// <summary>
        /// Creates a PeptideWithSetModifications object from a sequence string.
        /// Useful for reading in MetaMorpheus search engine output into mzLib objects
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

            Sequence = sequence;
            GetModsAfterDeserialization(allKnownMods, out _baseSequence);
            NumFixedMods = numFixedMods;
            DigestionParams = digestionParams;

            if (p != null)
            {
                ProteinAccession = p.Accession;
            }
        }

        public Dictionary<int, Modification> AllModsOneIsNterminus
        {
            get { return _allModsOneIsNterminus; }
        }

        public double MonoisotopicMass
        {
            get
            {
                if (!_monoisotopicMass.HasValue)
                {
                    _monoisotopicMass = WaterMonoisotopicMass;
                    foreach (var mod in AllModsOneIsNterminus.Values)
                    {
                        _monoisotopicMass += mod.MonoisotopicMass;
                    }
                    _monoisotopicMass += BaseSequence.Select(b => Residue.ResidueMonoisotopicMass[b]).Sum();
                }
                return (double)ClassExtensions.RoundedDouble(_monoisotopicMass.Value);
            }
        }

        public int NumMods
        {
            get
            {
                return AllModsOneIsNterminus.Count;
            }
        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (!HasChemicalFormulas.HasValue)
                {
                    HasChemicalFormulas = true;
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
                                yield return new Product(productType, terminalMasses[n], neutralLoss);
                            }
                        }
                    }

                    // "normal" fragment without neutral loss
                    yield return new Product(productType, terminalMasses[f], 0);
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

        public int NumVariableMods { get { return NumMods - NumFixedMods; } }

        public virtual string EssentialSequence(IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            string essentialSequence = BaseSequence;
            if (ModstoWritePruned != null)
            {
                var sbsequence = new StringBuilder();

                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                {
                    if (ModstoWritePruned.ContainsKey(pep_n_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_n_term_variable_mod.ModificationType + ":" + pep_n_term_variable_mod.Id + ']');
                    }
                }
                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // variable modification on this residue
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                    {
                        if (ModstoWritePruned.ContainsKey(residue_variable_mod.ModificationType))
                        {
                            sbsequence.Append('[' + residue_variable_mod.ModificationType + ":" + residue_variable_mod.Id + ']');
                        }
                    }
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                {
                    if (ModstoWritePruned.ContainsKey(pep_c_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_c_term_variable_mod.ModificationType + ":" + pep_c_term_variable_mod.Id + ']');
                    }
                }

                essentialSequence = sbsequence.ToString();
            }
            return essentialSequence;
        }

        public CompactPeptide CompactPeptide(FragmentationTerminus fragmentationTerminus)
        {
            if (_compactPeptides == null)
            {
                _compactPeptides = new Dictionary<FragmentationTerminus, CompactPeptide>();
            }
            if (_compactPeptides.TryGetValue(fragmentationTerminus, out CompactPeptide compactPeptide))
            {
                return compactPeptide;
            }
            else
            {
                CompactPeptide cp = new CompactPeptide(this, fragmentationTerminus);
                _compactPeptides.Add(fragmentationTerminus, cp);
                return cp;
            }
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

            var peptideWithLocalizedMass = new PeptideWithSetModifications(Protein, DigestionParams, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                PeptideDescription, MissedCleavages, dictWithLocalizedMass, NumFixedMods);

            return peptideWithLocalizedMass;
        }

        public override string ToString()
        {
            return Sequence + string.Join("\t", AllModsOneIsNterminus.Select(m => m.ToString()));
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;

            if (Protein == null && q.Protein == null)
            {
                return q.Sequence.Equals(Sequence);
            }

            return q != null
                && q.Sequence.Equals(Sequence)
                && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein
                && (q.Protein.Accession == null && Protein.Accession == null || q.Protein.Accession.Equals(Protein.Accession));
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }

        /// <summary>
        /// This should be run after deserialization of a PeptideWithSetModifications, in order to set its Protein and Modification objects, which were not serialized
        /// </summary>
        public static void SetNonSerializedPeptideInfo(Dictionary<string, Modification> idToMod, Dictionary<string, Protein> accessionToProtein, PeptideWithSetModifications peptide)
        {
            peptide.GetModsAfterDeserialization(idToMod, out string baseSequence);
            peptide.GetProteinAfterDeserialization(accessionToProtein);
        }

        private void GetModsAfterDeserialization(Dictionary<string, Modification> idToMod, out string baseSequence)
        {
            _allModsOneIsNterminus = new Dictionary<int, Modification>();
            StringBuilder baseSequenceSb = new StringBuilder();
            StringBuilder currentModification = new StringBuilder();
            int currentModificationLocation = 1;
            bool currentlyReadingMod = false;

            for (int r = 0; r < Sequence.Length; r++)
            {
                char c = Sequence[r];

                switch (c)
                {
                    case '[':
                        currentlyReadingMod = true;
                        break;

                    case ']':
                        string modId = null;

                        try
                        {
                            var split = currentModification.ToString().Split(new char[] { ':' });
                            string modType = split[0];
                            modId = split[1];
                        }
                        catch (Exception e)
                        {
                            throw new MzLibUtil.MzLibException("Error while trying to parse string into peptide: " + e.Message);
                        }

                        if (!idToMod.TryGetValue(modId, out Modification mod))
                        {
                            throw new MzLibUtil.MzLibException("Could not find modification while reading string: " + Sequence);
                        }

                        AllModsOneIsNterminus.Add(currentModificationLocation, mod);
                        currentlyReadingMod = false;
                        currentModification = new StringBuilder();
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
                subsequence.Append('[' + mod.ModificationType + ":" + mod.Id + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subsequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subsequence.Append('[' + mod.ModificationType + ":" + mod.Id + ']');
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                subsequence.Append('[' + mod.ModificationType + ":" + mod.Id + ']');
            }

            Sequence = subsequence.ToString();
        }
    }
}