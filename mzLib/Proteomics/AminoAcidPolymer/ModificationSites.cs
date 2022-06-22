// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ModificationSites.cs) is part of Proteomics.
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

namespace Proteomics.AminoAcidPolymer
{
    [Flags]
    public enum ModificationSites
    {
        None = 0,
        A = 1 << 0,
        R = 1 << 1,
        N = 1 << 2,
        D = 1 << 3,
        C = 1 << 4,
        E = 1 << 5,
        Q = 1 << 6,
        G = 1 << 7,
        H = 1 << 8,
        I = 1 << 9,
        L = 1 << 10,
        K = 1 << 11,
        M = 1 << 12,
        F = 1 << 13,
        P = 1 << 14,
        S = 1 << 15,
        T = 1 << 16,
        U = 1 << 17,
        W = 1 << 18,
        Y = 1 << 19,
        V = 1 << 20,
        NPep = 1 << 21,
        PepC = 1 << 22,
        NProt = 1 << 23,
        ProtC = 1 << 24,
        All = (1 << 25) - 1, // Handy way of setting all below the 24th bit
        NTerminus = NPep | NProt,
        TerminusC = PepC | ProtC,
        Any = 1 << 31 // Acts like none, but is equal to all
    }

    public static class ModificationSiteExtensions
    {
        public static IEnumerable<ModificationSites> EnumerateActiveSites(this ModificationSites sites)
        {
            foreach (ModificationSites site in Enum.GetValues(typeof(ModificationSites)))
            {
                if (site == ModificationSites.None)
                {
                    continue;
                }
                if ((sites & site) == site)
                {
                    yield return site;
                }
            }
        }

        public static bool ContainsSites(this ModificationSites sites, ModificationSites otherSites)
        {
            // By convention, if the other site is 'Any', they are always equal
            if (otherSites == ModificationSites.Any)
            {
                return true;
            }

            if (otherSites == ModificationSites.None)
            {
                return sites == ModificationSites.None;
            }

            return (~sites & otherSites) == ModificationSites.None;
        }
    }
}