// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
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
        /// The mass difference tolerance for having identical masses
        /// </summary>
        private const double MassEqualityEpsilon = 1e-10;

        /// <summary>
        /// Calculates m/z value for a given mass assuming charge comes from losing or gaining protons
        /// </summary>
        public static double ToMZ(this IHasMass objectWithMass, int charge)
        {
            if (objectWithMass == null)
                throw new ArgumentNullException("objectWithMass", "Cannot compute an MZ value for a null object");
            return ToMassToChargeRatio(objectWithMass.MonoisotopicMass, charge);
        }

        /// <summary>
        /// Calculates m/z value for a given mass assuming charge comes from losing or gaining protons
        /// </summary>
        public static double ToMassToChargeRatio(this double mass, int charge)
        {
            if (charge == 0)
                throw new DivideByZeroException("Charge cannot be zero");
            return mass / Math.Abs(charge) + Math.Sign(charge) * Constants.ProtonMass;
        }

        /// <summary>
        /// Determines the original mass from an m/z value, assuming charge comes from a proton
        /// </summary>
        public static double ToMass(this double massToChargeRatio, int charge)
        {
            if (charge == 0)
                throw new DivideByZeroException("Charge cannot be zero");
            return Math.Abs(charge) * massToChargeRatio - charge * Constants.ProtonMass;
        }

        /// <summary>
        /// Approximate equality check, error due to rounding
        /// </summary>
        public static bool MassEquals(this double mass1, double mass2, double epsilon)
        {
            return Math.Abs(mass1 - mass2) < epsilon;
        }

        /// <summary>
        /// Approximate equality check, error due to rounding
        /// </summary>
        public static bool MassEquals(this double mass1, double mass2)
        {
            return Math.Abs(mass1 - mass2) < MassEqualityEpsilon;
        }
    }
}