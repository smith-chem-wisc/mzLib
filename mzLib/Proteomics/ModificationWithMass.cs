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

        #endregion Public Fields

        #region Private Fields

        private readonly double neutralLoss;

        #endregion Private Fields

        #region Public Constructors

        public ModificationWithMass(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotDR)
        {
            this.monoisotopicMass = uniprotMM;
            this.neutralLoss = neutralLoss;
            this.massesObserved = massesObserved;
            this.diagnosticIons = diagnosticIons;
        }

        public ModificationWithMass(string id, Tuple<string, string> ac, string tg, ModificationSites pos, double mm, double neutralLoss)
            : base(id, ac, tg, pos)
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