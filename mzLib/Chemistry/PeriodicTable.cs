// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (PeriodicTable.cs) is part of Chemistry Library.
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

namespace Chemistry
{
    /// <summary>
    /// A static store of elements accessible by anyone
    /// </summary>
    public static class PeriodicTable
    {
        // Two datastores storing same elements! Need both for efficient access by both symbol and atomic number

        /// <summary>
        /// The internal dictionary housing elements, keyed by their unique atomic symbol
        /// </summary>
        private static readonly Dictionary<string, Element> _elements = new Dictionary<string, Element>();

        /// <summary>
        /// The internal dictionary housing elements, keyed by their unique atomic number
        /// </summary>
        private static readonly Element[] _elementsArray = new Element[Constants.MaximumNumberOfElementsAllowed];

        /// <summary>
        /// Populate the periodic table by calling this method
        /// </summary>
        public static void Add(Element element)
        {
            if (!_elements.ContainsKey(element.AtomicSymbol))
            {
                _elements.Add(element.AtomicSymbol, element);
                _elementsArray[element.AtomicNumber] = element;
            }
        }

        /// <summary>
        /// Fast method for getting an element by its atomic symbol
        /// </summary>
        public static Element GetElement(string atomicSymbol)
        {
            return _elements[atomicSymbol];
        }

        /// <summary>
        /// Fast method for getting an element by its atomic number
        /// </summary>
        public static Element GetElement(int atomicNumber)
        {
            return _elementsArray[atomicNumber];
        }

        /// <summary>
        /// Validates the abundances in the periodic table
        /// </summary>
        public static bool ValidateAbundances(double epsilon)
        {
            foreach (var e in _elements)
            {
                double totalAbundance = e.Value.Isotopes.Sum(b => b.RelativeAbundance);
                if (Math.Abs(totalAbundance - 1) > epsilon)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Validates the average masses in the periodic table
        /// </summary>
        public static bool ValidateAverageMasses(double epsilon)
        {
            foreach (var e in _elements)
            {
                double averageMass = e.Value.Isotopes.Sum(b => b.RelativeAbundance * b.AtomicMass);
                if (Math.Abs(averageMass - e.Value.AverageMass) / e.Value.AverageMass > epsilon)
                {
                    return false;
                }
            }
            return true;
        }
    }
}