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

using System;
using System.Collections.Generic;
using System.Linq;

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
    }
}