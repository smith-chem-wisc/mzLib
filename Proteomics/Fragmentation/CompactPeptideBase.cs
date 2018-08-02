using Chemistry;
using MassSpectrometry;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.Fragmentation
{
    [Serializable]
    public abstract class CompactPeptideBase : IEquatable<CompactPeptideBase>
    {
        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        private const int digitsForRoundingMasses = 9;
        private const double massTolForPeptideEquality = 1e-9;

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
                if (TerminalMasses == null)
                {
                    foreach (var mass in TerminalMasses)
                    {
                        result = (result * 31) ^ mass.GetHashCode();
                    }
                }
                return result;
            }
        }

        public bool Equals(CompactPeptideBase cp)
        {
            if (TerminalMasses != null && cp.TerminalMasses != null) //neither series is nulll
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods)
                        && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods))
                        || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                    && TerminalMasses.SequenceEqual(cp.TerminalMasses)
                    );
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
                for (int r = 0; r < peptide.Length; r++)
                {
                    mass += Residue.ResidueMonoisotopicMass[peptide[r]];

                    if (peptide.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification currentModification))
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }

                    if (r != peptide.Length - 1)
                    {
                        yield return new NeutralTerminusFragment(FragmentationTerminus.N, Math.Round(mass, digitsForRoundingMasses), r + 1);
                    }
                }
            }

            if (fragmentationTerminus == FragmentationTerminus.C || fragmentationTerminus == FragmentationTerminus.Both)
            {
                mass = 0;

                for (int r = peptide.Length - 1; r >= 0; r--)
                {
                    mass += Residue.ResidueMonoisotopicMass[peptide[r]];

                    if (peptide.AllModsOneIsNterminus.TryGetValue(r + 2, out Modification currentModification))
                    {
                        mass += (double)currentModification.MonoisotopicMass;
                    }

                    if (r != 0)
                    {
                        yield return new NeutralTerminusFragment(FragmentationTerminus.C, Math.Round(mass, digitsForRoundingMasses), r);
                    }
                }
            }
        }
    }
}