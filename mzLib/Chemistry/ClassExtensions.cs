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
using System.Collections.Generic;
using System.Linq;
using static System.Net.Mime.MediaTypeNames;

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

        public class TupleList<T1, T2> : List<Tuple<T1, T2>>
        {
            public void Add(T1 item, T2 item2)
            {
                Add(new Tuple<T1, T2>(item, item2));
            }
        }


        /// <summary>
        /// The mass difference tolerance for having identical masses
        /// </summary>
        public const double AbsoluteMassDifferenceAcceptor = 1e-10;

        public static bool MassEquals(this double mass1, IHasMass mass2, double epsilon = AbsoluteMassDifferenceAcceptor)
        {
            if (mass2 == null)
                return false;
            return mass1.MassEquals(mass2.MonoisotopicMass);
        }

        public static bool MassEquals(this double mass1, double mass2, double epsilon = AbsoluteMassDifferenceAcceptor)
        {
            return Math.Abs(mass1 - mass2) < epsilon;
        }

        public static bool MassEquals(this IHasMass mass1, double mass2, double epsilon = AbsoluteMassDifferenceAcceptor)
        {
            if (mass1 == null)
                return false;
            return mass1.MonoisotopicMass.MassEquals(mass2);
        }

        public static bool MassEquals(this IHasMass mass1, IHasMass mass2, double epsilon = AbsoluteMassDifferenceAcceptor)
        {
            if (mass1 == null || mass2 == null)
                return false;
            return mass1.MonoisotopicMass.MassEquals(mass2.MonoisotopicMass);
        }
    }
}