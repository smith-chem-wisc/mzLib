using System.Collections.Generic;
using Chemistry;
using Proteomics;
using UsefulProteomicsDatabases.Generated;
using System;

namespace UsefulProteomicsDatabases
{
    internal class ModificationWithMassAndCf : ModificationWithMass
    {
        private readonly ChemicalFormula chemicalFormula;

        public ModificationWithMassAndCf(string id, Tuple<string, string> ac, string tg, position_t pos, ChemicalFormula cf, double mm, double nl)
             : base(id, ac, tg, pos, mm, nl)
        {
            this.chemicalFormula = cf;
        }

        public ModificationWithMassAndCf(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, ChemicalFormula uniprotCF, double uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, double neutralLoss, IEnumerable<double> massesObserved, IEnumerable<double> diagnosticIons)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotMM, uniprotDR, neutralLoss, massesObserved, diagnosticIons)
        {
            this.chemicalFormula = uniprotCF;
        }
    }
}