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

        public ModificationWithMassAndCf(string id, Tuple<string, string> accession, ModificationMotif motif, ModificationSites site, ChemicalFormula chemicalFormula, double mm, IDictionary<string, IList<string>> linksToOtherDbs, List<double> neutralLosses, List<double> diagnosticIons, string modificationType)
            : base(id, accession, motif, site, mm, linksToOtherDbs, neutralLosses, diagnosticIons, modificationType)
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

        public override bool Equals(object o)
        {
            ModificationWithMassAndCf m = o as ModificationWithMassAndCf;
            return m == null ? false :
                base.Equals(m) && this.chemicalFormula.Equals(m.chemicalFormula);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() ^ chemicalFormula.GetHashCode();
        }

        #endregion Public Methods

    }
}