using Chemistry;
using Proteomics.AminoAcidPolymer;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.ProteolyticDigestion
{
    public class PeptideWithSetModifications : Peptide
    {
        /// <summary>
        /// dictionary of modifications on a peptide the N terminus is index 1
        /// key indicates which residue modification is on (with 1 being N terminus)
        /// </summary>
        public readonly Dictionary<int, ModificationWithMass> AllModsOneIsNterminus;
        public readonly int NumFixedMods;

        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private bool? HasChemicalFormulas;
        private readonly Dictionary<TerminusType, CompactPeptide> _compactPeptides = new Dictionary<TerminusType, CompactPeptide>();
        private string _sequence;
        private string _sequenceWithChemicalFormulas;
        private double? _monoisotopicMass;

        public PeptideWithSetModifications(Protein protein, int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, string peptideDescription, int missedCleavages,
            Dictionary<int, ModificationWithMass> allModsOneIsNterminus, int numFixedMods)
            : base(protein, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, peptideDescription)
        {
            AllModsOneIsNterminus = allModsOneIsNterminus;
            NumFixedMods = numFixedMods;
        }

        public PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, PeptideWithSetModifications everythingElseFromThisOne)
            : base(everythingElseFromThisOne.Protein, everythingElseFromThisOne.OneBasedStartResidueInProtein, everythingElseFromThisOne.OneBasedEndResidueInProtein,
                  everythingElseFromThisOne.MissedCleavages, everythingElseFromThisOne.PeptideDescription)
        {
            AllModsOneIsNterminus = modsFromThisOne.AllModsOneIsNterminus;
            NumFixedMods = modsFromThisOne.NumFixedMods;
        }

        public PeptideWithSetModifications(PeptideWithSetModifications modsFromThisOne, int proteinOneBasedStart, int proteinOneBasedEnd)
            : base(modsFromThisOne.Protein, proteinOneBasedStart, proteinOneBasedEnd, proteinOneBasedEnd - proteinOneBasedStart, modsFromThisOne.PeptideDescription)
        {
            AllModsOneIsNterminus = modsFromThisOne.AllModsOneIsNterminus
                .Where(b => b.Key > (1 + proteinOneBasedStart - modsFromThisOne.OneBasedStartResidueInProtein) 
                    && b.Key <= (2 + proteinOneBasedEnd - modsFromThisOne.OneBasedStartResidueInProtein))
                .ToDictionary(b => (b.Key + modsFromThisOne.OneBasedStartResidueInProtein - proteinOneBasedStart), b => b.Value);
        }

        public PeptideWithSetModifications(int numFixedMods, Protein protein, int proteinOneBasedStart, int proteinOneBasedEnd,
            Dictionary<int, ModificationWithMass> allModsOneIsNterminus = null, int missedCleavages = 0)
            : base(protein, proteinOneBasedStart, proteinOneBasedEnd, missedCleavages, null)
        {
            NumFixedMods = numFixedMods;
            AllModsOneIsNterminus = allModsOneIsNterminus ?? new Dictionary<int, ModificationWithMass>();
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
                        _monoisotopicMass += mod.monoisotopicMass;
                    }
                    _monoisotopicMass += BaseSequence.Select(b => Residue.ResidueMonoisotopicMass[b]).Sum();
                }
                return _monoisotopicMass.Value;
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
                    if (AllModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                    {
                        subsequence.Append('[' + pep_n_term_variable_mod.modificationType + ":" + pep_n_term_variable_mod.id + ']');
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        subsequence.Append(this[r]);
                        // variable modification on this residue
                        if (AllModsOneIsNterminus.TryGetValue(r + 2, out ModificationWithMass residue_variable_mod))
                        {
                            subsequence.Append('[' + residue_variable_mod.modificationType + ":" + residue_variable_mod.id + ']');
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (AllModsOneIsNterminus.TryGetValue(Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                    {
                        subsequence.Append('[' + pep_c_term_variable_mod.modificationType + ":" + pep_c_term_variable_mod.id + ']');
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
                    if (AllModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                    {
                        if (pep_n_term_variable_mod is ModificationWithMassAndCf jj)
                        {
                            subsequence.Append('[' + jj.chemicalFormula.Formula + ']');
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
                        if (AllModsOneIsNterminus.TryGetValue(r + 2, out ModificationWithMass residue_variable_mod))
                        {
                            if (residue_variable_mod is ModificationWithMassAndCf jj)
                            {
                                subsequence.Append('[' + jj.chemicalFormula.Formula + ']');
                            }
                            else
                            {
                                return null;
                            }
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (AllModsOneIsNterminus.TryGetValue(Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                    {
                        if (pep_c_term_variable_mod is ModificationWithMassAndCf jj)
                        {
                            subsequence.Append('[' + jj.chemicalFormula.Formula + ']');
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

        public int NumVariableMods { get { return NumMods - NumFixedMods; } }

        public virtual string EssentialSequence(IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            string essentialSequence = BaseSequence;
            if (ModstoWritePruned != null)
            {
                var sbsequence = new StringBuilder();

                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                {
                    if (ModstoWritePruned.ContainsKey(pep_n_term_variable_mod.modificationType))
                    {
                        sbsequence.Append('[' + pep_n_term_variable_mod.modificationType + ":" + pep_n_term_variable_mod.id + ']');
                    }
                }
                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // variable modification on this residue
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out ModificationWithMass residue_variable_mod))
                    {
                        if (ModstoWritePruned.ContainsKey(residue_variable_mod.modificationType))
                        {
                            sbsequence.Append('[' + residue_variable_mod.modificationType + ":" + residue_variable_mod.id + ']');
                        }
                    }
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                {
                    if (ModstoWritePruned.ContainsKey(pep_c_term_variable_mod.modificationType))
                    {
                        sbsequence.Append('[' + pep_c_term_variable_mod.modificationType + ":" + pep_c_term_variable_mod.id + ']');
                    }
                }

                essentialSequence = sbsequence.ToString();
            }
            return essentialSequence;
        }

        public CompactPeptide CompactPeptide(TerminusType terminusType)
        {
            if (_compactPeptides.TryGetValue(terminusType, out CompactPeptide compactPeptide))
            {
                return compactPeptide;
            }
            else
            {
                CompactPeptide cp = new CompactPeptide(this, terminusType);
                _compactPeptides.Add(terminusType, cp);
                return cp;
            }
        }

        public PeptideWithSetModifications Localize(int j, double massToLocalize)
        {
            var vvv = new Dictionary<int, ModificationWithMass>(AllModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (vvv.TryGetValue(j + 2, out ModificationWithMass modToReplace))
            {
                massOfExistingMod = modToReplace.monoisotopicMass;
                vvv.Remove(j + 2);
            }

            vvv.Add(j + 2, new ModificationWithMass(null, null, null, TerminusLocalization.Any, massToLocalize + massOfExistingMod));
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