using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class ModificationWithMass : Modification
    {

        #region Public Fields

        public readonly IEnumerable<double> massesObserved;
        public readonly double monoisotopicMass;
        public readonly IEnumerable<double> diagnosticIons;
        public readonly double neutralLoss;

        #endregion Public Fields

        #region Public Constructors

        public ModificationWithMass(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, string uniprotPP, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons, string database)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotDR, database)
        {
            this.monoisotopicMass = uniprotMM;
            this.neutralLoss = neutralLoss;
            this.massesObserved = massesObserved;
            this.diagnosticIons = diagnosticIons;
        }

        public ModificationWithMass(string id, Tuple<string, string> ac, string tg, ModificationSites pos, double mm, double neutralLoss, string database)
            : base(id, ac, tg, pos, null, database)
        {
            this.neutralLoss = neutralLoss;
            this.monoisotopicMass = mm;
            massesObserved = new HashSet<double> { mm };
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