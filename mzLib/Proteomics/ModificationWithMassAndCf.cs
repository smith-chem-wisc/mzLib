using System.Collections.Generic;
using Chemistry;
using Proteomics;
using System;
using System.Text;

namespace Proteomics
{
    public class ModificationWithMassAndCf : ModificationWithMass
    {
        public readonly ChemicalFormula chemicalFormula;

        public ModificationWithMassAndCf(string id, Tuple<string, string> ac, string tg, ModificationSites pos, ChemicalFormula cf, double mm, double nl)
             : base(id, ac, tg, pos, mm, nl)
        {
            this.chemicalFormula = cf;
        }

        public ModificationWithMassAndCf(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, ChemicalFormula uniprotCF, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotMM, uniprotDR, neutralLoss, massesObserved, diagnosticIons)
        {
            this.chemicalFormula = uniprotCF;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append("CF   " + chemicalFormula.Formula);
            return sb.ToString();
        }
    }
}