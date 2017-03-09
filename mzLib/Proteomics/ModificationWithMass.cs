using System;
using System.Collections.Generic;
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

        #endregion Public Methods

    }
}