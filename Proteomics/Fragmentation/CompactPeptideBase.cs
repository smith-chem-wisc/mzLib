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
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private const int digitsForRoundingMasses = 9;
        private const double massTolForPeptideEquality = 1e-9;

        public NeutralTerminusFragment[] TerminalMasses { get; protected set; } //

        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        public TheoreticalFragmentIon[] ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType dissociationType, int charge = 0)
        {

            List<ProductType> nTerminalProductTypes = DissociationTypeCollection.ProductsFromDissociationType[dissociationType].Except(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.C]).ToList();
            List<ProductType> cTerminalProductTypes = DissociationTypeCollection.ProductsFromDissociationType[dissociationType].Except(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[FragmentationTerminus.N]).ToList();

            List<double> calculated_N_terminalMasses = new List<double>();
            List<double> calculated_C_terminalMasses = new List<double>();

            if(NTerminalMasses != null && nTerminalProductTypes.Count > 0)
            {
                foreach (var NTerminalMass in NTerminalMasses)
                {
                    foreach (ProductType p in nTerminalProductTypes)
                    {
                        calculated_N_terminalMasses.Add(DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(NTerminalMass, p));
                    }                   
                }
            }

            if (CTerminalMasses != null && cTerminalProductTypes.Count > 0)
            {
                foreach (var CTerminalMass in CTerminalMasses)
                {
                    foreach (ProductType p in cTerminalProductTypes)
                    {
                        calculated_N_terminalMasses.Add(DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(CTerminalMass, p));
                    }
                }
            }

            return calculated_N_terminalMasses.Concat(calculated_C_terminalMasses).ToArray();
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

        protected static IEnumerable<NeutralTerminusFragment> ComputeFollowingFragmentMasses(PeptideWithSetModifications peptide, double prevMass, int residue, int direction, DissociationType dissociationType)//we're going to have to pass fragmentation type
        {
            do
            {
                if (residue != 0 && ((residue > 1 && direction == -1) || (residue != peptide.Length && direction == 1)))//This equates to true if you're at the beginning or end of the peptide and about to jump off. This should be true at the final step of the do loop.
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
                        yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(prevMass, digitsForRoundingMasses), "", residue) ;
                    }
                    else if (currentModification.NeutralLosses != null && ((residue > 1 && direction == -1) || (residue != peptide.Length && direction == 1))) // we don't want to consider neutral losses on the complete peptide
                    {
                        // neutral losses

                        List<double> theseNeutralLosses = new List<double>();
                        if (currentModification.NeutralLosses.TryGetValue(dissociationType, out var specificNeutralLosses))
                            theseNeutralLosses.AddRange(specificNeutralLosses);
                        if(dissociationType != DissociationType.AnyActivationType)
                            if (currentModification.NeutralLosses.TryGetValue(DissociationType.AnyActivationType, out var anyNeutralLosses))
                                theseNeutralLosses.AddRange(anyNeutralLosses);
                        if (theseNeutralLosses.Count > 0)
                        {
                            // return mass without neutral loss
                            prevMass += (double)currentModification.MonoisotopicMass;
                            yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(prevMass, digitsForRoundingMasses), "", residue) ;

                            foreach (double loss in theseNeutralLosses)
                            {
                                if (loss == 0)
                                {
                                    continue;
                                }

                                // return the current fragment minus neutral loss
                                yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(prevMass - loss, digitsForRoundingMasses), "-" + loss.ToString(), residue) ;

                                // generate the remainder of the series with the neutral loss
                                if ((residue > 1 && direction == -1) || (residue < (peptide.Length - 1) && direction == 1))
                                {
                                    foreach (var followingMass in ComputeFollowingFragmentMasses(peptide, prevMass, residue + direction, direction, dissociationType))
                                    {
                                        yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(followingMass.NeutralMass - loss, digitsForRoundingMasses), "-" + loss.ToString(), residue) ;
                                    }
                                }
                            }
                        }
                        else//there were neutral losses but they were the wrong kind so just add the mass of the mod.
                        {
                            prevMass += (double)currentModification.MonoisotopicMass;
                            yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(prevMass, digitsForRoundingMasses), "", ) ;
                        }
                    }
                }
                else if (residue != 0) // No modification exists
                {
                    yield return new NeutralTerminusFragment(GetTermus(direction), Math.Round(prevMass, digitsForRoundingMasses), "", ) ;
                }
                residue += direction;
            } while ((residue > 1 && direction == -1) || (residue < peptide.Length && direction == 1));
        }

        private static FragmentationTerminus GetTermus(int direction)
        {
            if (direction == 1)
                return FragmentationTerminus.N;
            else
                return FragmentationTerminus.C;
        }

    }
}