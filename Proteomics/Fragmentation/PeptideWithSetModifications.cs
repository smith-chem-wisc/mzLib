using Chemistry;
using MassSpectrometry;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.Fragmentation
{
    public class PeptideWithSetModifications : ProteolyticPeptide
    {
        /// <summary>
        /// dictionary of modifications on a peptide the N terminus is index 1
        /// key indicates which residue modification is on (with 1 being N terminus)
        /// </summary>
        public readonly Dictionary<int, Modification> AllModsOneIsNterminus;

        public readonly DigestionParams DigestionParams;
        public readonly int NumFixedMods;

        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private bool? HasChemicalFormulas;
        private readonly Dictionary<FragmentationTerminus, CompactPeptide> _compactPeptides = new Dictionary<FragmentationTerminus, CompactPeptide>();
        private string _sequence;
        private string _sequenceWithChemicalFormulas;
        private double? _monoisotopicMass;

        public PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, int proteinOneBasedStart, int proteinOneBasedEnd)
            : base(modsFromThisOne.Protein, proteinOneBasedStart, proteinOneBasedEnd, proteinOneBasedEnd - proteinOneBasedStart,
            modsFromThisOne.PeptideDescription)
        {
            AllModsOneIsNterminus = modsFromThisOne.AllModsOneIsNterminus
                .Where(b => b.Key > (1 + proteinOneBasedStart - modsFromThisOne.OneBasedStartResidueInProtein)
                    && b.Key <= (2 + proteinOneBasedEnd - modsFromThisOne.OneBasedStartResidueInProtein))
                .ToDictionary(b => (b.Key + modsFromThisOne.OneBasedStartResidueInProtein - proteinOneBasedStart), b => b.Value);
        }

        public PeptideWithSetModifications(int numFixedMods, Protein protein, int proteinOneBasedStart, int proteinOneBasedEnd,
            Dictionary<int, Modification> allModsOneIsNterminus = null, int missedCleavages = 0)
            : base(protein, proteinOneBasedStart, proteinOneBasedEnd, missedCleavages, null)
        {
            NumFixedMods = numFixedMods;
            AllModsOneIsNterminus = allModsOneIsNterminus ?? new Dictionary<int, Modification>();
        }

        public PeptideWithSetModifications(Protein protein, DigestionParams digestionParams, int oneBasedStartResidueInProtein,
            int oneBasedEndResidueInProtein, string peptideDescription, int missedCleavages,
           Dictionary<int, Modification> allModsOneIsNterminus, int numFixedMods)
           : base(protein, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, peptideDescription)
        {
            AllModsOneIsNterminus = allModsOneIsNterminus;
            NumFixedMods = numFixedMods;
            DigestionParams = digestionParams;
        }

        /// <summary>
        /// Creates a PeptideWithSetModifications object from a sequence string. 
        /// Useful for reading in MetaMorpheus search engine output into mzLib objects
        /// </summary>
        public PeptideWithSetModifications(string sequence, IEnumerable<Modification> allKnownModifications, int numFixedMods = 0,
            DigestionParams digestionParams = null, Protein p = null, int oneBasedStartResidueInProtein = int.MinValue,
            int oneBasedEndResidueInProtein = int.MinValue, int missedCleavages = int.MinValue, string peptideDescription = null)
            : base(p, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, peptideDescription)
        {
            if (sequence.Contains("|"))
            {
                throw new MzLibUtil.MzLibException("Ambiguous peptide cannot be parsed from string: " + sequence);
            }

            AllModsOneIsNterminus = new Dictionary<int, Modification>();
            StringBuilder baseSequence = new StringBuilder();
            StringBuilder currentModification = new StringBuilder();
            int currentModificationLocation = 1;
            bool currentlyReadingMod = false;

            for (int r = 0; r < sequence.Length; r++)
            {
                char c = sequence[r];

                switch (c)
                {
                    case '[':
                        currentlyReadingMod = true;
                        break;
                    case ']':

                        Modification mod = null;

                        try
                        {
                            var split = currentModification.ToString().Split(new char[] { ':' });
                            string modType = split[0];
                            string id = split[1];
                            mod = allKnownModifications.Where(m => m.Id == id && m.ModificationType == modType).FirstOrDefault();
                        }
                        catch (Exception e)
                        {
                            throw new MzLibUtil.MzLibException("Error while trying to parse string into peptide: " + e.Message);
                        }

                        if (mod == null)
                        {
                            throw new MzLibUtil.MzLibException("Could not find modification while reading string: " + sequence);
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
                            baseSequence.Append(c);
                        }
                        break;
                }
            }

            _baseSequence = baseSequence.ToString();
            NumFixedMods = numFixedMods;
            DigestionParams = digestionParams;
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

        public virtual string Sequence
        {
            get
            {
                if (_sequence == null)
                {
                    var subsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                    {
                        subsequence.Append('[' + pep_n_term_variable_mod.ModificationType + ":" + pep_n_term_variable_mod.Id + ']');
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        subsequence.Append(this[r]);
                        // variable modification on this residue
                        if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                        {
                            subsequence.Append('[' + residue_variable_mod.ModificationType + ":" + residue_variable_mod.Id + ']');
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                    {
                        subsequence.Append('[' + pep_c_term_variable_mod.ModificationType + ":" + pep_c_term_variable_mod.Id + ']');
                    }

                    _sequence = subsequence.ToString();
                }
                return _sequence;
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
        /// Generates theoretical fragment ions for given product types for this peptide
        /// </summary>
        public IEnumerable<NeutralTheoreticalProduct> GetTheoreticalFragments(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus)
        {
            var productCollection = TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[fragmentationTerminus].Intersect(DissociationTypeCollection.ProductsFromDissociationType[dissociationType]);

            foreach (var productType in productCollection)
            {
                // we're separating the N and C terminal masses and computing a separate compact peptide for each one
                // this speeds calculations up without producing unnecessary terminus fragment info
                FragmentationTerminus temp = TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType];

                foreach (var neutralTerminusFragment in CompactPeptide(temp).TerminalMasses)
                {
                    //TODO: compute mass of fragment given neutral loss
                    yield return new NeutralTheoreticalProduct(productType, neutralTerminusFragment, null);
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
            var vvv = new Dictionary<int, Modification>(AllModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out Modification modToReplace))
            {
                massOfExistingMod = (double)modToReplace.MonoisotopicMass;
                vvv.Remove(j + 2);
            }

            vvv.Add(j + 2, new Modification(_locationRestriction: "Anywhere.", _monoisotopicMass: massToLocalize + massOfExistingMod));
            var hm = new PeptideWithSetModifications(NumFixedMods, Protein, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein, vvv, MissedCleavages);

            return hm;
        }

        public override string ToString()
        {
            return Sequence + string.Join("\t", AllModsOneIsNterminus.Select(m => m.ToString()));
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;
            return q != null
                && q.Sequence.Equals(Sequence)
                && q.OneBasedStartResidueInProtein == OneBasedStartResidueInProtein
                && (q.Protein.Accession == null && Protein.Accession == null || q.Protein.Accession.Equals(Protein.Accession));
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode();
        }
    }
}