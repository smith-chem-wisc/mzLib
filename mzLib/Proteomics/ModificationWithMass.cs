using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

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

        public ModificationWithMass(string uniprotID, Tuple<string, string> uniprotAC, ModificationMotif uniprotTG, ModificationSites uniprotPP, double uniprotMM, IDictionary<string, IList<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons, string modificationType)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotDR, modificationType)
        {
            this.monoisotopicMass = uniprotMM;
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
            sb.AppendLine("NL   " + neutralLoss);
            sb.AppendLine("MO   " + string.Join(" or ", massesObserved));
            sb.AppendLine("DI   " + (diagnosticIons != null ? string.Join(" or ", diagnosticIons) : ""));
            sb.Append("MM   " + monoisotopicMass);
            return sb.ToString();
        }

        public override bool Equals(Modification m)
        {
            return base.Equals(m)
                && m as ModificationWithMass != null

                && (this.massesObserved == null && ((ModificationWithMass)m).massesObserved == null
                || this.massesObserved != null && ((ModificationWithMass)m).massesObserved != null
                && this.massesObserved.All(x => ((ModificationWithMass)m).massesObserved.Contains(x))
                && ((ModificationWithMass)m).massesObserved.All(x => this.massesObserved.Contains(x)))

                && (this.diagnosticIons == null && ((ModificationWithMass)m).diagnosticIons == null
                || this.diagnosticIons != null && ((ModificationWithMass)m).diagnosticIons != null
                && this.diagnosticIons.All(x => ((ModificationWithMass)m).diagnosticIons.Contains(x))
                && ((ModificationWithMass)m).diagnosticIons.All(x => this.diagnosticIons.Contains(x)))

                && this.monoisotopicMass == ((ModificationWithMass)m).monoisotopicMass
                && this.neutralLoss == ((ModificationWithMass)m).neutralLoss;
        }

        public override int GetCustomHashCode()
        {
            int hash = base.GetCustomHashCode()
                + sum_string_chars(monoisotopicMass.ToString())
                + sum_string_chars(monoisotopicMass.ToString());
            if (massesObserved != null) hash += massesObserved.Sum(x => sum_string_chars(x.ToString()));
            if (diagnosticIons != null) hash += diagnosticIons.Sum(x => sum_string_chars(x.ToString()));
            return hash;
        }

        #endregion Public Methods

        #region Private Methods

        private int sum_string_chars(string s)
        {
            if (s == null) return 0;
            int sum = int.MinValue + s.ToCharArray().Sum(c => 37 * c);
            return sum != 0 ? sum : sum + 1;
        }

        #endregion Private Methods

    }
}