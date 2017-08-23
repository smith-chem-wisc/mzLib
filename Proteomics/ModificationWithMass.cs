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
        public readonly List<double> diagnosticIons;
        public readonly List<double> neutralLosses;

        #endregion Public Fields

        #region Protected Fields

        protected const double tolForEquality = 1e-9;

        #endregion Protected Fields

        #region Public Constructors

        public ModificationWithMass(string id, string modificationType, ModificationMotif motif, TerminusLocalization terminusLocalization, double monoisotopicMass, IDictionary<string, IList<string>> externalDatabaseReferences = null, List<string> keywords = null, List<double> neutralLosses = null, List<double> diagnosticIons = null)
            : base(id, modificationType, motif, terminusLocalization, externalDatabaseReferences, keywords)
        {
            this.monoisotopicMass = monoisotopicMass;

            // Optional
            this.neutralLosses = neutralLosses ?? new List<double> { 0 };
            this.diagnosticIons = diagnosticIons ?? new List<double>();

            this.neutralLosses = this.neutralLosses.OrderBy(b => b).ToList();
            this.diagnosticIons = this.diagnosticIons.OrderBy(b => b).ToList();
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
            return m != null
                && base.Equals(m)
                && ApproxSequenceEqual(diagnosticIons, m.diagnosticIons, tolForEquality)
                && ApproxSequenceEqual(neutralLosses, m.neutralLosses, tolForEquality)
                && Math.Abs(monoisotopicMass - m.monoisotopicMass) < tolForEquality;
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() ^ diagnosticIons.Count ^ diagnosticIons.Count;
        }

        #endregion Public Methods

        #region Private Methods

        private bool ApproxSequenceEqual(List<double> a, List<double> b, double tolForEquality)
        {
            for (int i = 0; i < a.Count; i++)
                if (Math.Abs(a[i] - b[i]) >= tolForEquality)
                    return false;
            return true;
        }

        #endregion Private Methods
    }
}