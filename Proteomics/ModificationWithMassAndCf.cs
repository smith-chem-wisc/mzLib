using Chemistry;
using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;

namespace Proteomics
{
    public class ModificationWithMassAndCf : ModificationWithMass
    {
        #region Public Fields

        public readonly ChemicalFormula chemicalFormula;

        #endregion Public Fields

        #region Public Constructors

        public ModificationWithMassAndCf(string id, Tuple<string, string> accession, ModificationMotif motif, TerminusLocalization site, ChemicalFormula chemicalFormula, double mm, IDictionary<string, IList<string>> linksToOtherDbs, IEnumerable<double> neutralLosses, IEnumerable<double> diagnosticIons, string modificationType)
            : base(id, accession, motif, site, mm, linksToOtherDbs, neutralLosses, diagnosticIons, modificationType)
        {
            this.chemicalFormula = chemicalFormula;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            var baseString = base.ToString();
            if (Math.Abs(monoisotopicMass - chemicalFormula.MonoisotopicMass) < 1e-9)
            {
                var a = Regex.Match(baseString, @"MM.*$");
                baseString = Regex.Replace(baseString, @"MM.*$", "CF   " + chemicalFormula.Formula);
                sb.Append(baseString);
            }
            else
            {
                sb.AppendLine(baseString);
                sb.Append("CF   " + chemicalFormula.Formula);
            }
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            ModificationWithMassAndCf m = o as ModificationWithMassAndCf;
            return m == null ? false :
                base.Equals(m) && chemicalFormula.Equals(m.chemicalFormula);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() ^ chemicalFormula.GetHashCode();
        }

        #endregion Public Methods
    }
}