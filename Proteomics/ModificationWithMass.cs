using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationWithMass : ModificationWithLocation
    {

        #region Public Fields

        public readonly double monoisotopicMass;
        public readonly IEnumerable<double> diagnosticIons;
        public readonly IEnumerable<double> neutralLosses;

        #endregion Public Fields

        #region Public Constructors

        public ModificationWithMass(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites modificationSites, double monoisotopicMass, IDictionary<string, IList<string>> externalDatabaseReferences, IEnumerable<double> neutralLosses, IEnumerable<double> diagnosticIons, string modificationType)
            : base(id, accession, motif, modificationSites, externalDatabaseReferences, modificationType)
        {
            this.monoisotopicMass = monoisotopicMass;
            this.neutralLosses = neutralLosses;
            this.diagnosticIons = diagnosticIons;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            if (neutralLosses.Count() != 1 || neutralLosses.First() != 0)
                sb.AppendLine("NL   " + string.Join(" or ", neutralLosses.Select(b => b.ToString(CultureInfo.InvariantCulture))));
            if (diagnosticIons != null)
                sb.AppendLine("DI   " + string.Join(" or ", diagnosticIons.Select(b => b.ToString(CultureInfo.InvariantCulture))));
            sb.Append("MM   " + monoisotopicMass.ToString(CultureInfo.InvariantCulture));
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithMass m = o as ModificationWithMass;
            return m == null ? false :

                base.Equals(m)

                && (this.diagnosticIons == null && m.diagnosticIons == null
                || this.diagnosticIons != null && m.diagnosticIons != null
                && this.diagnosticIons.OrderBy(x => x).SequenceEqual(m.diagnosticIons.OrderBy(x => x)))

                && (this.neutralLosses == null && m.neutralLosses == null
                || this.neutralLosses != null && m.neutralLosses != null
                && this.neutralLosses.OrderBy(x => x).SequenceEqual(m.neutralLosses.OrderBy(x => x)))

                && this.monoisotopicMass == m.monoisotopicMass; ;
        }

        public override int GetHashCode()
        {
            int hash = base.GetHashCode() ^ monoisotopicMass.GetHashCode();
            if (neutralLosses != null) foreach (double x in neutralLosses) hash = hash ^ x.GetHashCode();
            if (diagnosticIons != null) foreach (double x in diagnosticIons) hash = hash ^ x.GetHashCode();
            return hash;
        }

        #endregion Public Methods

    }
}