using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class ModificationWithMass : ModificationWithLocation
    {

        #region Public Fields

        public readonly IEnumerable<double> massesObserved;
        public readonly double monoisotopicMass;
        public readonly IEnumerable<double> diagnosticIons;
        public readonly double neutralLoss;

        #endregion Public Fields

        #region Public Constructors

        public ModificationWithMass(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites modificationSites, double monoisotopicMass, IDictionary<string, IList<string>> externalDatabaseReferences, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons, string modificationType)
            : base(id, accession, motif, modificationSites, externalDatabaseReferences, modificationType)
        {
            this.monoisotopicMass = monoisotopicMass;
            this.neutralLoss = neutralLoss;
            this.massesObserved = massesObserved;
            this.diagnosticIons = diagnosticIons;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            if (neutralLoss != 0)
                sb.AppendLine("NL   " + neutralLoss);
            if (massesObserved.Count() != 1 || massesObserved.First() != monoisotopicMass)
                sb.AppendLine("MO   " + string.Join(" or ", massesObserved));
            if (diagnosticIons != null)
                sb.AppendLine("DI   " + string.Join(" or ", diagnosticIons));
            sb.Append("MM   " + monoisotopicMass);
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithMass m = o as ModificationWithMass;
            return m == null ? false :

                base.Equals(m)
                && (this.massesObserved == null && m.massesObserved == null
                || this.massesObserved != null && m.massesObserved != null
                && this.massesObserved.OrderBy(x => x).SequenceEqual(m.massesObserved.OrderBy(x => x)))

                && (this.diagnosticIons == null && m.diagnosticIons == null
                || this.diagnosticIons != null && m.diagnosticIons != null
                && this.diagnosticIons.OrderBy(x => x).SequenceEqual(m.diagnosticIons.OrderBy(x => x)))

                && this.monoisotopicMass == m.monoisotopicMass
                && this.neutralLoss == m.neutralLoss;
        }

        public override int GetHashCode()
        {
            int hash = base.GetHashCode() ^ monoisotopicMass.GetHashCode() ^ neutralLoss.GetHashCode();
            if (massesObserved != null) foreach (double x in massesObserved) hash = hash ^ x.GetHashCode();
            if (diagnosticIons != null) foreach (double x in diagnosticIons) hash = hash ^ x.GetHashCode();
            return hash;
        }

        #endregion Public Methods

    }
}