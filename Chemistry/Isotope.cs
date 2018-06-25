// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (Isotope.cs) is part of Chemistry Library.
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

using System.Globalization;

namespace Chemistry
{
    /// <summary>
    /// Represents a single isotope of a chemical element. Contains a unique number
    /// of protons and neutrons compared to every other isotope.
    /// </summary>
    public sealed class Isotope
    {
        /// <summary>
        /// Create a new isotope
        /// </summary>
        /// <param name="parentElement">The parent element of the isotope</param>
        /// <param name="massNumber">The mass number of the isotope</param>
        /// <param name="atomicMass">The atomic mass of the isotope</param>
        /// <param name="abundance">The natural relative abundance of the isotope</param>
        internal Isotope(Element parentElement, int massNumber, double atomicMass, double abundance)
        {
            Element = parentElement;
            MassNumber = massNumber;
            AtomicMass = atomicMass;
            RelativeAbundance = abundance;
        }

        /// <summary>
        /// The atomic number of the isotope's parent element (also the number of protons)
        /// </summary>
        public int AtomicNumber
        {
            get { return Element.AtomicNumber; }
        }

        /// <summary>
        /// The number of protons in this isotope
        /// </summary>
        public int Protons
        {
            get { return Element.AtomicNumber; }
        }

        /// <summary>
        /// The number of neutrons in this isotope
        /// </summary>
        public int Neutrons
        {
            get { return MassNumber - Element.AtomicNumber; }
        }

        /// <summary>
        /// The element this isotope is apart of (based on atomic number)
        /// </summary>
        public Element Element { get; }

        /// <summary>
        /// The atomic mass of this isotope (in unified atomic mass units)
        /// </summary>
        public double AtomicMass { get; }

        /// <summary>
        /// The total number of nucleons (protons and neutrons) in this isotope
        /// </summary>
        public int MassNumber { get; }

        /// <summary>
        /// The relative natural abundance of this isotope in nature (on Earth)
        /// </summary>
        public double RelativeAbundance { get; }

        /// <summary>
        /// Returns a textual representation of this isotope in the following format: H{1} He{4} O{16}
        /// </summary>
        /// <returns>The atomic symbol and mass number combined</returns>
        public override string ToString()
        {
            return string.Format(CultureInfo.InvariantCulture, "{0}{{{1}}}", Element.AtomicSymbol, MassNumber);
        }
    }
}