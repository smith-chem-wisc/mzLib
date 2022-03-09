// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (Fragment.cs) is part of Proteomics.
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
using System;
using System.Collections.Generic;
using System.Globalization;

namespace Proteomics.AminoAcidPolymer
{
    public class Fragment : IHasMass, IEquatable<Fragment>
    {
        public Fragment(FragmentTypes type, int number, double monoisotopicMass, AminoAcidPolymer parent)
        {
            FragmentType = type;
            Number = number;
            Parent = parent;
            MonoisotopicMass = monoisotopicMass;
        }

        public double MonoisotopicMass { get; private set; }

        public int Number { get; private set; }

        public AminoAcidPolymer Parent { get; private set; }

        public FragmentTypes FragmentType { get; private set; }

        public IEnumerable<IHasMass> Modifications
        {
            get
            {
                var mods = Parent.Modifications;
                if (FragmentType.GetTerminus() == Terminus.N)
                {
                    for (int i = 0; i <= Number; i++)
                    {
                        if (mods[i] != null)
                            yield return mods[i];
                    }
                }
                else
                {
                    int length = Parent.Length + 1;
                    for (int i = length - Number; i <= length; i++)
                    {
                        if (mods[i] != null)
                            yield return mods[i];
                    }
                }
            }
        }

        public string Sequence
        {
            get
            {
                string parentSeq = Parent.BaseSequence;
                if (FragmentType.GetTerminus() == Terminus.N)
                {
                    return parentSeq.Substring(0, Number);
                }

                return parentSeq.Substring(parentSeq.Length - Number, Number);
            }
        }

        public override string ToString()
        {
            return string.Format(CultureInfo.InvariantCulture, "{0}{1}", Enum.GetName(typeof(FragmentTypes), FragmentType), Number);
        }

        public override int GetHashCode()
        {
            return MonoisotopicMass.GetHashCode();
        }

        public bool Equals(Fragment other)
        {
            return FragmentType.Equals(other.FragmentType) && Number.Equals(other.Number) && Math.Abs(MonoisotopicMass - other.MonoisotopicMass) < 1e-9;
        }
    }
}