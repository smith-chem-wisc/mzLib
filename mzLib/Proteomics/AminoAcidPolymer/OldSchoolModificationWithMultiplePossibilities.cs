// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (Isotopologue.cs) is part of Proteomics.
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

using MzLibUtil;
using System.Collections;
using System.Collections.Generic;

namespace Proteomics.AminoAcidPolymer
{
    public class ModificationWithMultiplePossibilitiesCollection : OldSchoolModification, IEnumerable<OldSchoolModification>
    {
        private readonly SortedList<double, OldSchoolModification> _modifications;

        public ModificationWithMultiplePossibilitiesCollection(string name, ModificationSites sites)
            : base(0, name, sites)
        {
            _modifications = new SortedList<double, OldSchoolModification>();
        }

        public int Count
        {
            get { return _modifications.Count; }
        }

        public OldSchoolModification this[int index]
        {
            get { return _modifications.Values[index]; }
        }

        public void AddModification(OldSchoolModification modification)
        {
            if (!Sites.ContainsSites(modification.Sites))
                throw new MzLibException("Unable to add a modification with sites other than " + Sites);

            _modifications.Add(modification.MonoisotopicMass, modification);
        }

        public bool Contains(OldSchoolModification modification)
        {
            return _modifications.ContainsValue(modification);
        }

        public IEnumerator<OldSchoolModification> GetEnumerator()
        {
            return _modifications.Values.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return _modifications.Values.GetEnumerator();
        }
    }
}