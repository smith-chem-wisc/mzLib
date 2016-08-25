// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
// 
// This file (Modification.cs) is part of Proteomics.
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
using System.Globalization;

namespace Proteomics
{
    /// <summary>
    /// Represents a modification with a mass and name and default amino acid sites of modification
    /// </summary>
    public class Modification : IHasMass, IEquatable<Modification>
    {
        /// <summary>
        /// The name of the modification
        /// </summary>
        public string Name { get; protected set; }

        /// <summary>
        /// The monoisotopic mass of the modification, commoningly known as the delta mass
        /// </summary>
        public double MonoisotopicMass { get; protected set; }

        /// <summary>
        /// The potentially modified sites of this modification
        /// </summary>
        public ModificationSites Sites { get; set; }

        /// <summary>
        /// Displays the name of the mod and the sites it modified in a formated string
        /// </summary>
        public string NameAndSites
        {
            get { return string.Format(CultureInfo.InvariantCulture, "{0} ({1})", Name, Sites); }
        }

        public Modification(Modification modification)
            : this(modification.MonoisotopicMass, modification.Name, modification.Sites)
        {
        }

        public Modification()
            : this(0.0, "", ModificationSites.Any)
        {

        }
        public Modification(double monoMass)
            : this(monoMass, "", ModificationSites.Any)
        {

        }

        public Modification(double monoMass, string name)
            : this(monoMass, name, ModificationSites.Any)
        {

        }
        public Modification(double monoMass, string name, ModificationSites sites)
        {
            MonoisotopicMass = monoMass;
            Name = name;
            Sites = sites;
        }

        public override string ToString()
        {
            return Name;
        }

        public override int GetHashCode()
        {
            return MonoisotopicMass.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Modification modObj = obj as Modification;
            return modObj != null && Equals(modObj);
        }


        public bool Equals(Modification other)
        {
            // TODO: Is this enough? Does it guarantee equality?
            if (ReferenceEquals(this, other))
                return true;

            if (!MonoisotopicMass.MassEquals(other.MonoisotopicMass))
                return false;

            if (!Name.Equals(other.Name))
                return false;

            if (!Sites.Equals(other.Sites))
                return false;

            return true;
        }
    }
}