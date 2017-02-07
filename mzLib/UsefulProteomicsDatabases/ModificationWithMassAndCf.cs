using System.Collections.Generic;
using Chemistry;
using Proteomics;

namespace UsefulProteomicsDatabases
{
    internal class ModificationWithMassAndCf : ModificationWithMass
    {
        private ChemicalFormula uniprotCF;

        public ModificationWithMassAndCf(string uniprotID, string uniprotAC, ModificationSites uniprotTG, ModificationSites uniprotPP, ChemicalFormula uniprotCF, double? uniprotMM, Dictionary<string, HashSet<string>> uniprotDR, IEnumerable<double> massesOnFragments, IEnumerable<double> massesObserved, IEnumerable<string> motifs)
            : base(uniprotID, uniprotAC, uniprotTG, uniprotPP, uniprotMM, uniprotDR, massesOnFragments, massesObserved, motifs)
        {
            this.uniprotCF = uniprotCF;
        }
    }
}