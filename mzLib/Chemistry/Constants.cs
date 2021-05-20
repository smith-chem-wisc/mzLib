// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (Constants.cs) is part of Chemistry Library.
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

namespace Chemistry
{
    /// <summary>
    /// A collection of immutable constants and physical properties.
    /// Masses are given for the most abundant isotope unless otherwise stated
    ///
    /// Sources include:
    /// http://physics.nist.gov/cuu/Constants/index.html
    /// </summary>
    public static class Constants
    {
        /// <summary>
        /// The mass of the subatomic particle with a single elementary charge in atomic
        /// units (u)
        /// </summary>
        public const double ProtonMass = 1.007276466879;

        /// <summary>
        /// The electron mass in atomic units (u)
        /// </summary>
        public const double ElectronMass = 5.48579909070e-4;

        /// <summary>
        /// The difference between the carbon 13 isotope and carbon 12 isotope in
        /// atomic units (u)
        /// </summary>
        public const double C13MinusC12 = 1.00335483810;

        /// <summary>
        /// The largest number of elements to consider
        /// </summary>
        internal const int MaximumNumberOfElementsAllowed = 128;

        /// <summary>
        /// The largest mass number
        /// </summary>
        internal const int MaximumMassNumberPossible = 295;

        internal const int CarbonAtomicNumber = 6;

        internal const int HydrogenAtomicNumber = 1;
    }
}