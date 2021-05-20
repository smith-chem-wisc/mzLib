// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (Peptide.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using System.Collections.Generic;
using System.Linq;

namespace Proteomics.AminoAcidPolymer
{
    public class Peptide : AminoAcidPolymer
    {
        public Peptide()
        {
        }

        public Peptide(string sequence) : base(sequence)
        {
        }

        public Peptide(AminoAcidPolymer aminoAcidPolymer)
                            : this(aminoAcidPolymer, true)
        {
        }

        /// <summary>
        /// Create a new peptide based on another amino acid polymer
        /// </summary>
        /// <param name="aminoAcidPolymer">The other amino acid polymer to copy</param>
        /// <param name="includeModifications">Whether to copy the modifications to the new peptide</param>
        public Peptide(AminoAcidPolymer aminoAcidPolymer, bool includeModifications)
            : base(aminoAcidPolymer, includeModifications)
        {
            Parent = aminoAcidPolymer;
            StartResidue = 0;
            EndResidue = Length - 1;
        }

        public Peptide(AminoAcidPolymer aminoAcidPolymer, int firstResidue, int length)
                    : this(aminoAcidPolymer, firstResidue, length, true)
        {
        }

        public Peptide(AminoAcidPolymer aminoAcidPolymer, int firstResidue, int length, bool includeModifications)
                    : base(aminoAcidPolymer, firstResidue, length, includeModifications)
        {
            Parent = aminoAcidPolymer;
            StartResidue = firstResidue;
            EndResidue = firstResidue + length - 1;
            PreviousResidue = aminoAcidPolymer.GetResidue(StartResidue - 1);
            NextResidue = aminoAcidPolymer.GetResidue(EndResidue + 1);
        }

        /// <summary>
        /// The amino acid number this peptide is located in its parent
        /// </summary>
        public int StartResidue { get; set; }

        /// <summary>
        /// The amino acid number this peptide is located in its parent
        /// </summary>
        public int EndResidue { get; set; }

        /// <summary>
        /// The amino acid polymer this peptide came from. Could be null
        /// </summary>
        public AminoAcidPolymer Parent { get; set; }

        /// <summary>
        /// The preceding amino acid in its parent
        /// </summary>
        public Residue PreviousResidue { get; set; }

        /// <summary>
        /// The next amino acid in its parent
        /// </summary>
        public Residue NextResidue { get; set; }

        public IEnumerable<Peptide> GenerateAllModificationCombinations()
        {
            // Get all the modifications that are isotopologues
            var isotopologues = GetUniqueModifications<ModificationWithMultiplePossibilitiesCollection>().ToArray();

            // Base condition, no more isotopologues to make, so just return
            if (isotopologues.Length < 1)
            {
                yield break;
            }

            // Grab the the first isotopologue
            ModificationWithMultiplePossibilitiesCollection isotopologue = isotopologues[0];

            // Loop over each modification in the isotopologue
            foreach (OldSchoolModification mod in isotopologue)
            {
                // Create a clone of the peptide, cloning modifications as well.
                Peptide peptide = new Peptide(this);

                // Replace the base isotopologue mod with the specific version
                peptide.ReplaceModification(isotopologue, mod);

                // There were more than one isotopologue, so we must go deeper
                if (isotopologues.Length > 1)
                {
                    // Call the same rotuine on the newly generate peptide that has one less isotopologue
                    foreach (var subpeptide in peptide.GenerateAllModificationCombinations())
                    {
                        yield return subpeptide;
                    }
                }
                else
                {
                    // Return this peptide
                    yield return peptide;
                }
            }
        }

        public Peptide GetSubPeptide(int firstResidue, int length)
        {
            return new Peptide(this, firstResidue, length);
        }
    }
}