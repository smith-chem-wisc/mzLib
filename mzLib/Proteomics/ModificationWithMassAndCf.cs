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

        public ModificationWithMassAndCf(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites site, ChemicalFormula chemicalFormula, double mm, IDictionary<string, IList<string>> linksToOtherDbs, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons, string database)
            : base(id, accession, motif, site, mm, linksToOtherDbs, neutralLoss, massesObserved, diagnosticIons, database)
        {
            this.chemicalFormula = chemicalFormula;
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