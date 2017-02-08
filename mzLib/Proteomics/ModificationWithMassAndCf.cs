using Chemistry;
using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class ModificationWithMassAndCf : ModificationWithMass
    {
        #region Public Fields

        public readonly ChemicalFormula chemicalFormula;

        #endregion Public Fields

        #region Public Constructors

        public ModificationWithMassAndCf(string id, Tuple<string, string> ac, string tg, ModificationSites pos, ChemicalFormula cf, double mm, double nl, string database)
             : base(id, ac, tg, pos, mm, nl, database)
        {
            this.chemicalFormula = cf;
        }

        public ModificationWithMassAndCf(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, ChemicalFormula uniprotCF, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons, string database)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotMM, uniprotDR, neutralLoss, massesObserved, diagnosticIons, database)
        {
            this.chemicalFormula = uniprotCF;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("CF   " + chemicalFormula.Formula);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}