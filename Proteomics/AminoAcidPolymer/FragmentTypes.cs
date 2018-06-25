// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (FragmentTypes.cs) is part of Proteomics.
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

using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;

namespace Proteomics.AminoAcidPolymer
{
    [Flags]
    public enum FragmentTypes
    {
        None = 0,
        a = 1 << 0,
        adot = 1 << 1,
        b = 1 << 2,
        bdot = 1 << 3,
        c = 1 << 4,
        cdot = 1 << 5,
        x = 1 << 6,
        xdot = 1 << 7,
        y = 1 << 8,
        ydot = 1 << 9,
        z = 1 << 10,
        zdot = 1 << 11,
        Internal = 1 << 12,
        All = (1 << 12) - 1, // Handy way of setting all below the 12th bit
    }

    public static class FragmentTypesExtension
    {
        private static readonly Dictionary<FragmentTypes, ChemicalFormula> FragmentIonCaps = new Dictionary<FragmentTypes, ChemicalFormula>
        {
            {FragmentTypes.a, ChemicalFormula.ParseFormula("C-1H-1O-1")},
            {FragmentTypes.adot, ChemicalFormula.ParseFormula("C-1O-1")},
            {FragmentTypes.b, ChemicalFormula.ParseFormula("H-1")},
            {FragmentTypes.bdot, new ChemicalFormula()},
            {FragmentTypes.c, ChemicalFormula.ParseFormula("NH2")},
            {FragmentTypes.cdot, ChemicalFormula.ParseFormula("NH3")},
            {FragmentTypes.x, ChemicalFormula.ParseFormula("COH-1")},
            {FragmentTypes.xdot, ChemicalFormula.ParseFormula("CO")},
            {FragmentTypes.y, ChemicalFormula.ParseFormula("H")},
            {FragmentTypes.ydot, ChemicalFormula.ParseFormula("H2")},
            {FragmentTypes.z, ChemicalFormula.ParseFormula("N-1H-2")},
            {FragmentTypes.zdot, ChemicalFormula.ParseFormula("N-1H-1")}
        };

        public static IEnumerable<FragmentTypes> GetIndividualFragmentTypes(this FragmentTypes fragmentTypes)
        {
            foreach (FragmentTypes site in Enum.GetValues(typeof(FragmentTypes)))
            {
                if (site == FragmentTypes.None || site == FragmentTypes.All || site == FragmentTypes.Internal)
                {
                    continue;
                }
                if ((fragmentTypes & site) == site)
                {
                    yield return site;
                }
            }
        }

        public static Terminus GetTerminus(this FragmentTypes fragmentType)
        {
            // Super handy: http://stackoverflow.com/questions/4624248/c-logical-riddle-with-bit-operations-only-one-bit-is-set
            if (fragmentType == FragmentTypes.None || (fragmentType & (fragmentType - 1)) != FragmentTypes.None)
            {
                throw new MzLibException("Fragment Type must be a single value to determine the terminus");
            }
            return fragmentType >= FragmentTypes.x ? Terminus.C : Terminus.N;
        }

        public static ChemicalFormula GetIonCap(this FragmentTypes fragmentType)
        {
            if (fragmentType == FragmentTypes.None || (fragmentType & (fragmentType - 1)) != FragmentTypes.None)
            {
                throw new MzLibException("Fragment Type must be a single value to determine the ion cap");
            }
            return FragmentIonCaps[fragmentType];
        }
    }
}