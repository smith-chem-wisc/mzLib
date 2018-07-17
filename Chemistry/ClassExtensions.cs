// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (MassExtensions.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>.

using System;

namespace Chemistry
{
    public static class ClassExtensions
    {
        /// <summary>
        /// Calculates m/z value for a given mass assuming charge comes from losing or gaining protons
        /// </summary>
        public static double ToMz(this IHasMass objectWithMass, int charge)
        {
            return ToMz(objectWithMass.MonoisotopicMass, charge);
        }

        /// <summary>
        /// Calculates m/z value for a given mass assuming charge comes from losing or gaining protons
        /// </summary>
        public static double ToMz(this double mass, int charge)
        {
            return mass / Math.Abs(charge) + Math.Sign(charge) * Constants.ProtonMass;
        }

        /// <summary>
        /// Determines the original mass from an m/z value, assuming charge comes from a proton
        /// </summary>
        public static double ToMass(this double massToChargeRatio, int charge)
        {
            return Math.Abs(charge) * massToChargeRatio - charge * Constants.ProtonMass;
        }

        public static double? RoundedDouble(this double? myNumber, int places = 9)
        {
            if (myNumber != null)
            {
                myNumber = Math.Round((double)myNumber, places, MidpointRounding.AwayFromZero);
            }
            return myNumber;
        }

    }
}