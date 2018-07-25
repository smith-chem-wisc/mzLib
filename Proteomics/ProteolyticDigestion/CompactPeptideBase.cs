using Chemistry;
using MassSpectrometry;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics.ProteolyticDigestion
{
    [Serializable]
    public abstract class CompactPeptideBase : IEquatable<CompactPeptideBase>
    {
        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private const int digitsForRoundingMasses = 9;
        private const double massTolForPeptideEquality = 1e-9;

        public double[] CTerminalMasses { get; protected set; } //
        public double[] NTerminalMasses { get; protected set; }
        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        public double[] ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
            {
                throw new NotImplementedException();
            }
            if (containsBnoB1)
            {
                massLen += NTerminalMasses.Length - 1;
            }
            else if (containsB)
            {
                massLen += NTerminalMasses.Length;
            }
            if (containsC)
            {
                massLen += NTerminalMasses.Length;
            }
            if (containsX)
            {
                throw new NotImplementedException();
            }
            if (containsY)
            {
                massLen += CTerminalMasses.Length;
            }
            if (containsZdot)
            {
                massLen += CTerminalMasses.Length;
            }

            if (massLen < 0)
                return new double[0];

            double[] massesToReturn = new double[massLen];

            int i = 0;
            if (NTerminalMasses != null)
            {
                for (int j = 0; j < NTerminalMasses.Length; j++)
                {
                    var hm = NTerminalMasses[j];
                    if (containsB || (containsBnoB1 && j > 0))
                    {
                        massesToReturn[i] = ClassExtensions.RoundedDouble(hm).Value;
                        i++;
                    }
                    if (containsC)
                    {
                        massesToReturn[i] = ClassExtensions.RoundedDouble(hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass).Value;
                        i++;
                    }
                }
            }
            if (CTerminalMasses != null)
            {
                foreach (double hm in CTerminalMasses)
                {
                    if (containsY)
                    {
                        massesToReturn[i] = ClassExtensions.RoundedDouble(hm + waterMonoisotopicMass).Value;
                        i++;
                    }
                    if (containsZdot)
                    {
                        massesToReturn[i] = ClassExtensions.RoundedDouble(hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass).Value;
                        i++;
                    }
                }
            }
            return massesToReturn;
        }

        /// <summary>
        /// Sometimes says not equal when in reality should be equal, due to rounding errors. Small but annoying bug. Careful when fixing! Make sure Indexing runs at a reasonable speed.
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
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
                if (CTerminalMasses == null)
                {
                    foreach (double b in NTerminalMasses)
                    {
                        result = (result * 31) ^ b.GetHashCode();
                    }
                }
                else
                {
                    foreach (double b in CTerminalMasses)
                    {
                        result = (result * 31) ^ b.GetHashCode();
                    }
                }
                return result;
            }
        }

        public bool Equals(CompactPeptideBase cp)
        {
            if (CTerminalMasses != null && cp.CTerminalMasses != null)
            {
                if (NTerminalMasses != null && cp.NTerminalMasses != null) //neither series is nulll
                {
                    return (
                        ((double.IsNaN(MonoisotopicMassIncludingFixedMods)
                            && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods))
                            || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                        && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                        && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                        );
                }
                else //No N-terminal ions
                {
                    return (
                        ((double.IsNaN(MonoisotopicMassIncludingFixedMods)
                            && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods))
                            || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                        && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                        );
                }
            }
            else if (NTerminalMasses != null && cp.NTerminalMasses != null) //No C-terminal ions
            {
                return (
                    ((double.IsNaN(MonoisotopicMassIncludingFixedMods)
                        && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods))
                        || (Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < massTolForPeptideEquality))
                    && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                    );
            }
            else //Cannot compare
            {
                return false;
            }
        }

        protected static IEnumerable<double> ComputeFollowingFragmentMasses(PeptideWithSetModifications peptide, double prevMass, int residue, int direction, DissociationType dissociationType)//we're going to have to pass fragmentation type
        {
            do
            {
                if(residue != 0 && ((residue > 1 && direction == -1) || (residue != peptide.Length && direction == 1)))//This equates to true if you're at the beginning or end of the peptide and about to jump off. This should be true at the final step of the do loop.
                {
                    prevMass += Residue.ResidueMonoisotopicMass[peptide[residue - 1]];
                }

                // If modification exists on the particular residue that is being evaluated in the do loop.
                if (peptide.AllModsOneIsNterminus.TryGetValue(residue + 1, out Modification currentModification))
                {
                    if ((currentModification.NeutralLosses == null || currentModification.NeutralLosses.Count == 0) && residue != 0 && residue != peptide.Length + 1)
                    {
                        // no neutral losses - just add the modification's mass
                        prevMass += (double)currentModification.MonoisotopicMass;
                        yield return Math.Round(prevMass, digitsForRoundingMasses);
                    }
                    else if (currentModification.NeutralLosses != null && ((residue > 1 && direction == -1) || (residue != peptide.Length && direction == 1))) // we don't want to consider neutral losses on the complete peptide
                    {
                        // neutral losses
                        if (currentModification.NeutralLosses.TryGetValue(dissociationType, out var neutralLosses))
                        {
                            // return mass without neutral loss
                            prevMass += (double)currentModification.MonoisotopicMass;
                            yield return Math.Round(prevMass, digitsForRoundingMasses);

                            foreach (double loss in neutralLosses)
                            {
                                if (loss == 0)
                                {
                                    continue;
                                }

                                // return the current fragment minus neutral loss
                                yield return Math.Round(prevMass - loss, digitsForRoundingMasses);

                                // generate the remainder of the series with the neutral loss
                                if ((residue > 1 && direction == -1) || (residue < (peptide.Length-1) && direction == 1))
                                {
                                    foreach (var followingMass in ComputeFollowingFragmentMasses(peptide, prevMass, residue + direction, direction, dissociationType))
                                    {
                                        yield return Math.Round(followingMass - loss, digitsForRoundingMasses);
                                    }
                                }
                            }
                        }
                    }
                }
                else if (residue != 0) // No modification exists
                {
                    yield return Math.Round(prevMass, digitsForRoundingMasses);
                }
                residue += direction;
            } while ((residue > 1 && direction == -1) || (residue < peptide.Length && direction == 1));
        }
    }
}