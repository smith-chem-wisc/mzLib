// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (AminoAcidPolymerExtensions.cs) is part of Proteomics.
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

using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;

namespace Proteomics.AminoAcidPolymer
{
    public static class AminoAcidPolymerExtensions
    {
        public static double GetSequenceCoverageFraction(this AminoAcidPolymer baseSequence, IEnumerable<AminoAcidPolymer> sequences)
        {
            return GetSequenceCoverageFraction(baseSequence, sequences, true);
        }

        public static double GetSequenceCoverageFraction(this AminoAcidPolymer baseSequence, IEnumerable<AminoAcidPolymer> sequences, bool useLeucineSequence)
        {
            int[] counts = baseSequence.GetSequenceCoverage(sequences, useLeucineSequence);
            return ((double)counts.Count(x => x > 0)) / baseSequence.Length;
        }

        public static int[] GetSequenceCoverage(this AminoAcidPolymer baseSequence, IEnumerable<AminoAcidPolymer> sequences)
        {
            return GetSequenceCoverage(baseSequence, sequences, true);
        }

        public static int[] GetSequenceCoverage(this AminoAcidPolymer baseSequence, IEnumerable<AminoAcidPolymer> allPolymers, bool useLeucineSequence)
        {
            int[] bits = new int[baseSequence.Length];

            string masterSequence = useLeucineSequence ? baseSequence.BaseLeucineSequence : baseSequence.BaseSequence;

            foreach (AminoAcidPolymer polymer in allPolymers)
            {
                string seq = useLeucineSequence ? polymer.BaseLeucineSequence : polymer.BaseSequence;

                int startIndex = 0;
                while (true)
                {
                    int index = masterSequence.IndexOf(seq, startIndex, StringComparison.Ordinal);

                    if (index < 0)
                    {
                        break;
                    }

                    for (int aa = index; aa < index + polymer.Length; aa++)
                    {
                        bits[aa]++;
                    }

                    startIndex = index + 1;
                }
            }
            return bits;
        }

        /// <summary>
        /// Checks whether a given string represents a valid peptide sequence without modifications. 
        /// Valid sequences contain only residues in the ResidueDictionary.
        /// Note: the residue dictionary can be externally modified. Unusual amino acid letters can and do appear in peptide sequences.
        /// </summary>
        /// <param name="baseSequence"> Sequence to be checked </param>
        /// <returns> True if the sequence is valid. False otherwise. </returns>
        public static bool AllSequenceResiduesAreValid(this string baseSequence)
        {
            if (baseSequence.IsNullOrEmpty()) return false;

            return !baseSequence.ToCharArray() //the input is the base sequence.
                .Distinct() //this line eliminates any duplicated amino acids to save search time.
                .Except(Residue.ResiduesDictionary.Values.Select(l => l.Letter)) //This line removes from the base sequence array any residues that are known in the dictionary
                .Any(); //If there are any leftovers, then that means that there are base sequence characters not in the dictionary. therefore the sequence is not valide.
                // please note the "!" at the start of the whole linq statement.
        }
    }
}