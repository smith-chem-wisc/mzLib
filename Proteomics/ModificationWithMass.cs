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

        public ModificationWithMass(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites modificationSites, double monoisotopicMass, IDictionary<string, IList<string>> externalDatabaseReferences, List<double> neutralLosses, List<double> diagnosticIons, string modificationType)
            : base(id, accession, motif, modificationSites, externalDatabaseReferences, modificationType)
        {
            this.monoisotopicMass = monoisotopicMass;
            if (neutralLosses == null || neutralLosses.Count() == 0)
                neutralLosses = new List<double> { 0 };
            this.neutralLosses = neutralLosses ?? new List<double> { 0 };
            this.diagnosticIons = diagnosticIons ?? new List<double>();
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            if (neutralLosses.Count() != 1 || neutralLosses.First() != 0)
                sb.AppendLine("NL   " + string.Join(" or ", neutralLosses.Select(b => b.ToString(CultureInfo.InvariantCulture))));
            if (diagnosticIons.Count() > 0)
                sb.AppendLine("DI   " + string.Join(" or ", diagnosticIons.Select(b => b.ToString(CultureInfo.InvariantCulture))));
            sb.Append("MM   " + monoisotopicMass.ToString(CultureInfo.InvariantCulture));
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithMass m = o as ModificationWithMass;
            return m == null ? false :
                base.Equals(m)
                && (this.diagnosticIons.OrderBy(x => x).SequenceEqual(m.diagnosticIons.OrderBy(x => x)))
                && (this.neutralLosses.OrderBy(x => x).SequenceEqual(m.neutralLosses.OrderBy(x => x)))
                && Math.Abs(this.monoisotopicMass - m.monoisotopicMass) < 1e-9;
        }

        public override int GetHashCode()
        {
            int hash = base.GetHashCode();
            foreach (double x in neutralLosses) hash = hash ^ x.GetHashCode();
            foreach (double x in diagnosticIons) hash = hash ^ x.GetHashCode();
            return hash;
        }

        #endregion Public Methods

    }
}