using Chemistry;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.Fragmentation
{
    [Serializable]
    public abstract class CompactPeptideBase : IEquatable<CompactPeptideBase>
    {
        public NeutralTerminusFragment[] TerminalMasses { get; protected set; }
        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        /// <summary>
        /// Sometimes says not equal when in reality should be equal, due to rounding errors. Small but annoying bug. Careful when fixing! Make sure Indexing runs at a reasonable speed.
        /// </summary>
        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptideBase;
            return cp != null && Equals(cp);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var result = 0;
                if (TerminalMasses != null)
                {
                    foreach (var mass in TerminalMasses)
                    {
                        result += (result * 31) ^ mass.GetHashCode();
                    }
                }
                return result;
            }
        }

        public bool Equals(CompactPeptideBase cp)
        {
            if (TerminalMasses != null && cp.TerminalMasses != null) //neither series is nulll
            {
                return TerminalMasses.SequenceEqual(cp.TerminalMasses);
            }
            else //Cannot compare
            {
                return false;
            }
        }

        protected static IEnumerable<NeutralTerminusFragment> ComputeNeutralTerminusFragments(PeptideWithSetModifications peptide, FragmentationTerminus fragmentationTerminus)
        {
            double mass = 0;

            if (fragmentationTerminus == FragmentationTerminus.N || fragmentationTerminus == FragmentationTerminus.Both)
            {
                for (int r = 0; r <= peptide.Length - 1; r++)//This is a zero based indexed for residues. The index of the first amino acid in the peptide is 0.
                {
                    mass += Residue.ResidueMonoisotopicMass[peptide[r]];//This is a zero based indexed for residues. The index of the first amino acid in the peptide is 0.

                    // side-chain mod
                    if (peptide.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification currentModification))//This is a one based index. The index of the fragment from the first amino acid is 1.
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }

                    // N-terminal mod
                    if (r == 0 && peptide.AllModsOneIsNterminus.TryGetValue(1, out currentModification))
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }
                    
                    if (r != peptide.Length - 1)
                    {
                        yield return new NeutralTerminusFragment(FragmentationTerminus.N, mass, r + 1, r + 1);//This is a one based index. The index of the fragment from the first amino acid is 1.
                    }
                }
            }

            if (fragmentationTerminus == FragmentationTerminus.C || fragmentationTerminus == FragmentationTerminus.Both)
            {
                mass = 0;

                for (int r = peptide.Length - 1; r >= 0; r--)
                {
                    mass += Residue.ResidueMonoisotopicMass[peptide[r]];

                    // side-chain mod
                    if (peptide.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification currentModification))
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }
                    
                    // C-terminal mod
                    if (r == peptide.Length - 1 && peptide.AllModsOneIsNterminus.TryGetValue(peptide.Length + 2, out currentModification))
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }

                    if (r != -1)
                    {
                        yield return new NeutralTerminusFragment(FragmentationTerminus.C, mass, peptide.Length - r, r + 1);
                    }
                }
            }
        }
    }
}
