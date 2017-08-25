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

        public ModificationWithMassAndCf(string id, string modificationType, ModificationMotif motif, TerminusLocalization terminusLocalization, ChemicalFormula chemicalFormula, double? mm = null, IDictionary<string, IList<string>> linksToOtherDbs = null, List<string> keywords = null, List<double> neutralLosses = null, List<double> diagnosticIons = null)
            : base(id, modificationType, motif, terminusLocalization, mm ?? chemicalFormula.MonoisotopicMass, linksToOtherDbs, keywords, neutralLosses, diagnosticIons)
        {
            this.chemicalFormula = chemicalFormula;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            var baseString = base.ToString();
            if (Math.Abs(monoisotopicMass - chemicalFormula.MonoisotopicMass) < tolForEquality)
            {
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
            return m != null
               && base.Equals(m)
               && chemicalFormula.Equals(m.chemicalFormula);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() ^ chemicalFormula.GetHashCode();
        }

        #endregion Public Methods
    }
}