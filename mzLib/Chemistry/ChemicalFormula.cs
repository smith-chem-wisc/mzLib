// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (ChemicalFormula.cs) is part of Chemistry Library.
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

using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Chemistry
{
    /// <summary>
    /// A chemical formula class. This does NOT correspond to a physical object. A physical object can have a chemical formula.
    /// Formula can change!!! If isotopes or elements are changed.
    /// </summary>
    [Serializable]
    public sealed class ChemicalFormula : IEquatable<ChemicalFormula>
    {
        // Main data stores, the isotopes and elements

        /// <summary>
        /// A regular expression for matching chemical formulas such as: C2C{13}3H5NO5
        /// \s* (at end as well) allows for optional spacing among the elements, i.e. C2 C{13}3 H5 N O5
        /// The first group is the only non-optional group and that handles the chemical symbol: H, He, etc..
        /// The second group is optional, which handles isotopes of elements: C{13} means carbon-13, while C is the carbon element with unspecified mass number
        /// The third group is optional and indicates if we are adding or subtracting the elements form the formula, C-2C{13}5 would mean first subtract 2 carbons and then add 5 carbon-13
        /// The fourth group is optional and represents the number of isotopes or elements to add, if not present it assumes 1: H2O means 2 Hydrogen and 1 Oxygen
        /// Modified from: http://stackoverflow.com/questions/4116786/parsing-a-chemical-formula-from-a-string-in-c
        /// </summary>
        private static readonly Regex FormulaRegex = new Regex(@"\s*([A-Z][a-z]*)(?:\{([0-9]+)\})?(-)?([0-9]+)?\s*", RegexOptions.Compiled);

        /// <summary>
        /// A wrapper for the formula regex that validates if a string is in the correct chemical formula format or not
        /// </summary>
        private static readonly Regex ValidateFormulaRegex = new Regex("^(" + FormulaRegex + ")+$", RegexOptions.Compiled);

        private string formulaString;

        public ChemicalFormula()
        {
            Isotopes = new Dictionary<Isotope, int>();
            Elements = new Dictionary<Element, int>();
        }

        public ChemicalFormula(ChemicalFormula capFormula)
        {
            Isotopes = new Dictionary<Isotope, int>(capFormula.Isotopes);
            Elements = new Dictionary<Element, int>(capFormula.Elements);
        }

        /// <summary>
        /// Gets the average mass of this chemical formula
        /// </summary>
        public double AverageMass
        {
            get
            {
                return Isotopes.Sum(b => b.Key.AtomicMass * b.Value) + Elements.Sum(b => b.Key.AverageMass * b.Value);
            }
        }

        /// <summary>
        /// Gets the monoisotopic mass of this chemical formula: for elements use the principle isotope mass, not average mass
        /// </summary>
        public double MonoisotopicMass
        {
            get
            {
                return Isotopes.Sum(b => b.Key.AtomicMass * b.Value) + Elements.Sum(b => b.Key.PrincipalIsotope.AtomicMass * b.Value);
            }
        }

        /// <summary>
        /// Gets the number of atoms in this chemical formula
        /// </summary>
        public int AtomCount
        {
            get
            {
                return Isotopes.Sum(b => b.Value) + Elements.Sum(b => b.Value);
            }
        }

        /// <summary>
        /// Gets the number of unique chemical elements in this chemical formula
        /// </summary>
        public int NumberOfUniqueElementsByAtomicNumber
        {
            get
            {
                HashSet<int> ok = new HashSet<int>();
                foreach (var i in Isotopes)
                    ok.Add(i.Key.AtomicNumber);
                foreach (var i in Elements)
                    ok.Add(i.Key.AtomicNumber);
                return ok.Count;
            }
        }

        /// <summary>
        /// Gets the number of unique chemical isotopes in this chemical formula
        /// </summary>
        public int NumberOfUniqueIsotopes
        {
            get
            {
                return Isotopes.Count;
            }
        }

        /// <summary>
        /// Gets the string representation (Hill Notation) of this chemical formula
        /// </summary>
        public string Formula
        {
            get
            {
                if (formulaString == null)
                    formulaString = GetHillNotation();
                return formulaString;
            }
        }

        public int ProtonCount
        {
            get
            {
                return Isotopes.Sum(b => b.Key.Protons * b.Value) + Elements.Sum(b => b.Key.Protons * b.Value);
            }
        }

        /// <summary>
        /// The ratio of the number of Carbon to Hydrogen in this chemical formula
        /// </summary>
        /// <returns></returns>
        public double HydrogenCarbonRatio
        {
            get
            {
                int carbonCount = CountWithIsotopes("C");
                int hydrogenCount = CountWithIsotopes("H");
                return hydrogenCount / (double)carbonCount;
            }
        }

        internal Dictionary<Isotope, int> Isotopes { get; private set; }
        internal Dictionary<Element, int> Elements { get; private set; }

        public static ChemicalFormula Combine(IEnumerable<IHasChemicalFormula> formulas)
        {
            ChemicalFormula returnFormula = new ChemicalFormula();
            foreach (IHasChemicalFormula iformula in formulas)
                returnFormula.Add(iformula);
            return returnFormula;
        }

        /// <summary>
        /// Parses a string representation of chemical formula and adds the elements
        /// to this chemical formula.
        /// Use brackets for isotopes (C6H12N2O -> C{13}6H12N2O)
        /// </summary>
        /// <param name="formula">the Chemical Formula to parse</param>
        public static ChemicalFormula ParseFormula(string formula)
        {
            ChemicalFormula f = new ChemicalFormula();

            if (!ValidateFormulaRegex.IsMatch(formula))
                throw new MzLibException("Input string for chemical formula was in an incorrect format: " + formula);

            foreach (Match match in FormulaRegex.Matches(formula))
            {
                string chemsym = match.Groups[1].Value; // Group 1: Chemical Symbol

                Element element = PeriodicTable.GetElement(chemsym);

                int sign = match.Groups[3].Success ? // Group 3 (optional): Negative Sign
                    -1 :
                    1;

                int numofelem = match.Groups[4].Success ? // Group 4 (optional): Number of Elements
                    int.Parse(match.Groups[4].Value, CultureInfo.InvariantCulture) :
                    1;

                if (match.Groups[2].Success) // Group 2 (optional): Isotope Mass Number
                {
                    // Adding isotope!
                    f.Add(element[int.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture)], sign * numofelem);
                }
                else
                {
                    // Adding element!
                    f.Add(element, numofelem * sign);
                }
            }
            return f;
        }

        public int NeutronCount()
        {
            if (Elements.Count > 0)
                throw new MzLibException("Cannot know for sure what the number of neutrons is!");
            return Isotopes.Sum(b => b.Key.Neutrons * b.Value);
        }

        /// <summary>
        /// Replaces one isotope with another.
        /// Replacement happens on a 1 to 1 basis, i.e., if you remove 5 you will add 5
        /// </summary>
        /// <param name="isotopeToRemove">The isotope to remove</param>
        /// <param name="isotopeToAdd">The isotope to add</param>
        public void Replace(Isotope isotopeToRemove, Isotope isotopeToAdd)
        {
            int numberRemoved = Remove(isotopeToRemove);
            Add(isotopeToAdd, numberRemoved);
        }

        /// <summary>
        /// Add a chemical formula containing object to this chemical formula
        /// </summary>
        /// <param name="item">The object that contains a chemical formula</param>
        public void Add(IHasChemicalFormula item)
        {
            Add(item.ThisChemicalFormula);
        }

        public void Multiply(int multiplier)
        {
            List<Element> keys = new List<Element>(Elements.Keys);
            foreach (var key in keys)
                Elements[key] *= multiplier;
            List<Isotope> keysIsotope = new List<Isotope>(Isotopes.Keys);
            foreach (var key in keysIsotope)
                Isotopes[key] *= multiplier;
        }

        /// <summary>
        /// Add a chemical formula to this chemical formula.
        /// </summary>
        /// <param name="formula">The chemical formula to add to this</param>
        public void Add(ChemicalFormula formula)
        {
            foreach (var e in formula.Elements)
            {
                Add(e.Key, e.Value);
            }
            foreach (var i in formula.Isotopes)
            {
                Add(i.Key, i.Value);
            }
        }

        /// <summary>
        /// Add the principal isotope of the element to this chemical formula
        /// given its chemical symbol
        /// </summary>
        /// <param name="element">The chemical element to add</param>
        /// <param name="count">The number of the element to add</param>
        public void AddPrincipalIsotopesOf(Element element, int count)
        {
            Isotope isotope = element.PrincipalIsotope;
            Add(isotope, count);
        }

        public void Add(Element element, int count)
        {
            if (count == 0)
                return;
            if (!Elements.ContainsKey(element))
                Elements.Add(element, count);
            else
            {
                Elements[element] += count;
                if (Elements[element] <= 0)
                    Elements.Remove(element);
            }
            formulaString = null;
        }

        /// <summary>
        /// Add an isotope to this chemical formula
        /// </summary>
        /// <param name="isotope">The isotope to add</param>
        /// <param name="count">The number of the isotope to add</param>
        public void Add(Isotope isotope, int count)
        {
            if (count == 0)
                return;
            if (!Isotopes.ContainsKey(isotope))
                Isotopes.Add(isotope, count);
            else
            {
                Isotopes[isotope] += count;
                if (Isotopes[isotope] <= 0)
                    Isotopes.Remove(isotope);
            }
            formulaString = null;
        }

        /// <summary>
        /// Remove a chemical formula containing object from this chemical formula
        /// </summary>
        /// <param name="item">The object that contains a chemical formula</param>
        public void Remove(IHasChemicalFormula item)
        {
            Remove(item.ThisChemicalFormula);
        }

        /// <summary>
        /// Remove a chemical formula from this chemical formula
        /// </summary>
        /// <param name="formula">The chemical formula to remove</param>
        public void Remove(ChemicalFormula formula)
        {
            foreach (var e in formula.Elements)
                Remove(e.Key, e.Value);
            foreach (var i in formula.Isotopes)
                Remove(i.Key, i.Value);
        }

        /// <summary>
        /// Remove the provided number of elements (not isotopes!) from formula
        /// </summary>
        /// <param name="element">The chemical element to remove</param>
        /// <param name="count">The number of elements to remove</param>
        public void Remove(Element element, int count)
        {
            Add(element, -count);
        }

        /// <summary>
        /// Remove a isotope from this chemical formula
        /// </summary>
        /// <param name="isotope">The isotope to remove</param>
        /// <param name="count">The number of isotopes to remove</param>
        public void Remove(Isotope isotope, int count)
        {
            Add(isotope, -count);
        }

        /// <summary>
        /// Completely removes a particular isotope from this chemical formula.
        /// </summary>
        /// <param name="isotope">The isotope to remove</param>
        /// <returns>Number of removed isotopes</returns>
        public int Remove(Isotope isotope)
        {
            int count = Isotopes[isotope];
            Add(isotope, -count);
            return count;
        }

        /// <summary>
        /// Remove all the isotopes of an chemical element from this
        /// chemical formula
        /// </summary>
        /// <param name="element">The chemical element to remove</param>
        /// <returns>Number of removed isotopes</returns>
        public int RemoveIsotopesOf(Element element)
        {
            int count = Elements[element];
            Add(element, -count);
            foreach (var k in Isotopes.Where(b => b.Key.Element == element).ToList())
                Isotopes.Remove(k.Key);
            return count;
        }

        /// <summary>
        /// Remove all isotopes and elements
        /// </summary>
        public void Clear()
        {
            Isotopes = new Dictionary<Isotope, int>();
            Elements = new Dictionary<Element, int>();
            formulaString = null;
        }

        /// <summary>
        /// Checks if any isotope of the specified element is present in this chemical formula
        /// </summary>
        /// <param name="element">The element to look for</param>
        /// <returns>True if there is a non-zero number of the element in this formula</returns>
        public bool ContainsIsotopesOf(Element element)
        {
            return CountWithIsotopes(element) != 0;
        }

        public bool IsSubsetOf(ChemicalFormula formula)
        {
            return formula.IsSupersetOf(this);
        }

        /// <summary>
        /// Checks whether this formula contains all the isotopes of the specified formula
        /// MIGHT CONSIDER ELEMENTS TO BE SUPERSET OF ISOTOPES IF NEEDED!!!
        /// Right now they are distinct
        /// </summary>
        /// <param name="formula"></param>
        /// <returns></returns>
        public bool IsSupersetOf(ChemicalFormula formula)
        {
            foreach (var aa in formula.Elements)
                if (!Elements.ContainsKey(aa.Key) || aa.Value > Elements[aa.Key])
                    return false;
            foreach (var aa in formula.Isotopes)
                if (!Isotopes.ContainsKey(aa.Key) || aa.Value > Isotopes[aa.Key])
                    return false;
            return true;
        }

        public bool ContainsSpecificIsotope(Element element, int atomicNumber)
        {
            return CountSpecificIsotopes(element, atomicNumber) != 0;
        }

        /// <summary>
        /// Checks if the isotope is present in this chemical formula
        /// </summary>
        /// <param name="isotope">The isotope to look for</param>
        /// <returns>True if there is a non-negative number of the isotope in this formula</returns>
        public bool ContainsSpecificIsotope(Isotope isotope)
        {
            return CountSpecificIsotopes(isotope) != 0;
        }

        /// <summary>
        /// Return the number of given isotopes in this chemical fomrula
        /// </summary>
        /// <param name="isotope"></param>
        /// <returns></returns>
        public int CountSpecificIsotopes(Isotope isotope)
        {
            return (Isotopes.TryGetValue(isotope, out int isotopeCount) ? isotopeCount : 0);
        }

        /// <summary>
        /// Count the number of isotopes and elements from this element that are
        /// present in this chemical formula
        /// </summary>
        /// <param name="element">The element to search for</param>
        /// <returns>The total number of all the element isotopes in this chemical formula</returns>
        public int CountWithIsotopes(Element element)
        {
            var isotopeCount = element.Isotopes.Sum(isotope => CountSpecificIsotopes(isotope));
            return isotopeCount + (Elements.TryGetValue(element, out int ElementCount) ? ElementCount : 0);
        }

        public int CountSpecificIsotopes(Element element, int massNumber)
        {
            Isotope isotope = element[massNumber];
            return CountSpecificIsotopes(isotope);
        }

        public override int GetHashCode()
        {
            return Tuple.Create(Isotopes.Sum(b => b.Key.AtomicMass * b.Value), Elements.Sum(b => b.Key.AverageMass * b.Value)).GetHashCode();
        }

        public bool Equals(ChemicalFormula other)
        {
            if (other == null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (Math.Abs(MonoisotopicMass - other.MonoisotopicMass) > 1e-9)
                return false;
            if (Math.Abs(AverageMass - other.AverageMass) > 1e-9)
                return false;
            return true;
        }

        /// <summary>
        /// Produces the Hill Notation of the chemical formula
        /// </summary>
        private string GetHillNotation()
        {
            StringBuilder s = new StringBuilder();

            // Find carbons
            if (Elements.ContainsKey(PeriodicTable.GetElement(Constants.CarbonAtomicNumber)))
                s.Append("C" + (Elements[PeriodicTable.GetElement(Constants.CarbonAtomicNumber)] == 1 ? "" : "" + Elements[PeriodicTable.GetElement(Constants.CarbonAtomicNumber)]));

            // Find carbon isotopes
            foreach (var i in PeriodicTable.GetElement(Constants.CarbonAtomicNumber).Isotopes)
                if (Isotopes.ContainsKey(i))
                    s.Append("C{" + i.MassNumber + "}" + (Isotopes[i] == 1 ? "" : "" + Isotopes[i]));

            // Find hydrogens
            if (Elements.ContainsKey(PeriodicTable.GetElement(Constants.HydrogenAtomicNumber)))
                s.Append("H" + (Elements[PeriodicTable.GetElement(Constants.HydrogenAtomicNumber)] == 1 ? "" : "" + Elements[PeriodicTable.GetElement(Constants.HydrogenAtomicNumber)]));

            // Find hydrogen isotopes
            foreach (var i in PeriodicTable.GetElement(Constants.HydrogenAtomicNumber).Isotopes)
                if (Isotopes.ContainsKey(i))
                    s.Append("H{" + i.MassNumber + "}" + (Isotopes[i] == 1 ? "" : "" + Isotopes[i]));

            List<string> otherParts = new List<string>();

            // Find other elements
            foreach (var i in Elements)
                if (i.Key != PeriodicTable.GetElement(Constants.CarbonAtomicNumber) && i.Key != PeriodicTable.GetElement(Constants.HydrogenAtomicNumber))
                    otherParts.Add(i.Key.AtomicSymbol + (i.Value == 1 ? "" : "" + i.Value));

            // Find other isotopes
            foreach (var i in Isotopes)
                if (i.Key.Element != PeriodicTable.GetElement(Constants.CarbonAtomicNumber) && i.Key.Element != PeriodicTable.GetElement(Constants.HydrogenAtomicNumber))
                    otherParts.Add(i.Key.Element.AtomicSymbol + "{" + i.Key.MassNumber + "}" + (i.Value == 1 ? "" : "" + i.Value));

            otherParts.Sort();
            return s + string.Join("", otherParts);
        }
    }
}