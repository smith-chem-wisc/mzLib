using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace Transcriptomics
{
    public class Modification : IHasMass, IEquatable<Modification>
    {
        /// <summary>
        /// The default empty modification
        /// </summary>
        public static readonly Modification Empty = new Modification();

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
        public ModificationSite Site { get; set; }

        /// <summary>
        /// Displays the name of the mod and the sites it modified in a formated string
        /// </summary>
        public string NameAndSites => $"{Name} ({Site})";

        public Modification(Modification modification)
            : this(modification.MonoisotopicMass, modification.Name, modification.Site)
        {
        }

        public Modification(double monoMass = 0.0, string name = "", ModificationSite sites = ModificationSite.Any)
        {
            MonoisotopicMass = monoMass;
            Name = name;
            Site = sites;
        }

        #region Interface Implemenations and Overrides

        public override int GetHashCode()
        {
            return MonoisotopicMass.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
                return false;
            Modification modObj = obj as Modification;
            return modObj != null && Equals(modObj);
        }

        public bool Equals(Modification other)
        {
            if (other == null)
                return false;

            if (ReferenceEquals(this, other))
                return true;

            if (!this.MassEquals(other))
                return false;

            if (!Name.Equals(other.Name))
                return false;

            if (!Site.Equals(other.Site))
                return false;

            return true;
        }

        public override string ToString()
        {
            return Name;
        }

        #endregion

    }
}
